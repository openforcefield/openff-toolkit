# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Property calculator 'server' side API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import json
import logging
import os
import struct
import uuid
from os import path
from typing import Dict, List, Any

from pydantic import BaseModel
from simtk import unit
from tornado.ioloop import IOLoop, PeriodicCallback
from tornado.iostream import IOStream, StreamClosedError
from tornado.tcpserver import TCPServer

from openforcefield.properties import PhysicalProperty
from openforcefield.properties.estimator.client import PropertyEstimatorDataModel, PropertyEstimatorOptions
from openforcefield.properties.estimator.workflow.protocols import ProtocolPath, PropertyCalculatorException
from openforcefield.properties.estimator.layers import available_layers
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils.serialization import serialize_quantity
from .message_types import PropertyEstimatorMessageTypes

# Needed for server-client communication.
int_struct = struct.Struct("<i")

unpack_int = int_struct.unpack
pack_int = int_struct.pack


# =============================================================================================
# Property Runner
# =============================================================================================

class PropertyRunnerDataModel(BaseModel):
    """Represents a data packet to be calculated by the runner, along with
    the options which should be used when running the calculations.

    Attributes
    ----------
    id: str
        A unique id assigned to this calculation packet.
    queued_properties: list of PhysicalProperty
        A list of physical properties waiting to be calculated.
    calculated_properties: list of PhysicalProperty
        A list of physical properties which have been calculated.
    options: PropertyEstimatorOptions
        The options used to calculate the properties.
    parameter_set_id: str
        The unique server side id of the force field parameters used to calculate the properties.
    """

    id: str

    queued_properties: List[PhysicalProperty] = []
    calculated_properties: List[PhysicalProperty] = []

    options: PropertyEstimatorOptions = None

    parameter_set_id: str = None

    stored_data: Dict[str, Any] = {}

    class Config:

        arbitrary_types_allowed = True

        json_encoders = {
            unit.Quantity: lambda v: serialize_quantity(v),
            ProtocolPath: lambda v: v.full_path
        }


class PropertyCalculationRunner(TCPServer):
    """The class responsible for coordinating all calculations to be ran using
    the property estimator, in addition to deciding at which fidelity a property
    will be calculated.

    It acts as a server, which receives submitted jobs from clients
    launched via the property estimator.

    Notes
    -----
    Methods to handle the TCP messages are based on the StackOverflow response from
    A. Jesse Jiryu Davis: https://stackoverflow.com/a/40257248

    Examples
    --------
    Setting up a general server instance using a dask LocalCluster backend:

    >>> # Create the backend which will be responsible for distributing the calculations
    >>> from openforcefield.properties.estimator.backends import DaskLocalClusterBackend
    >>> calculation_backend = DaskLocalClusterBackend(1, 1)
    >>>
    >>> # Calculate the backend which will be responsible for storing and retrieving
    >>> # the data from previous calculations
    >>> from openforcefield.properties.estimator.storage import LocalFileStorage
    >>> storage_backend = LocalFileStorage()
    >>>
    >>> # Create the server to which all calculations will be submitted
    >>> from openforcefield.properties.estimator.runner import PropertyCalculationRunner
    >>> property_server = PropertyCalculationRunner(calculation_backend, storage_backend)
    >>>
    >>> # Instruct the server to listen for incoming submissions
    >>> property_server.run_until_killed()
    """

    @property
    def queued_calculations(self):
        return self._queued_calculations

    @property
    def finished_calculations(self):
        return self._finished_calculations

    def __init__(self, calculation_backend, storage_backend,
                 port=8000, working_directory='working-data'):
        """Constructs a new PropertyCalculationRunner object.

        Parameters
        ----------
        calculation_backend: PropertyEstimatorBackend
            The backend to use for executing calculations.
        storage_backend: PropertyEstimatorStorage
            The backend to use for storing information from any calculations.
        port: int
            The port one which to listen for incoming client
            requests.
        working_directory: str
            The local directory in which to store all local, temporary calculation data.
        """
        self._calculation_backend = calculation_backend
        self._storage_backend = storage_backend

        self.working_directory = working_directory

        if not path.isdir(self.working_directory):
            os.makedirs(self.working_directory)

        self._port = port

        self._queued_calculations = {}
        self._finished_calculations = {}

        self._periodic_loops = []

        super().__init__()

        # self.listen(self._port)
        self.bind(self._port)
        self.start(1)

        calculation_backend.start()

    async def handle_job_submission(self, stream, address, message_length):
        """An asynchronous routine for handling the receiving and processing
        of job submissions from a client.

        Parameters
        ----------
        stream: IOStream
            An IO stream used to pass messages between the
            server and client.
        address: str
            The address from which the request came.
        message_length: int
            The length of the message being recieved.
        """

        logging.info('Received job request from {}'.format(address))

        # Read the incoming request from the server. The first four bytes
        # of the response should be the length of the message being sent.

        # Decode the client submission json.
        encoded_json = await stream.read_bytes(message_length)
        json_model = encoded_json.decode()

        # TODO: Add exeception handling so the server can gracefully reject bad json.
        client_data_model = PropertyEstimatorDataModel.parse_raw(json_model)

        runner_data_model = self._prepare_data_model(client_data_model)
        should_launch = True

        # Make sure this job is not already in the queue.
        calculation_id = None
        cached_data_id = runner_data_model.id

        for existing_id in self._queued_calculations:

            runner_data_model.id = existing_id

            if runner_data_model.json() == self._queued_calculations[existing_id].json():

                calculation_id = existing_id
                should_launch = False

                break

        runner_data_model.id = cached_data_id

        if calculation_id is None:
            # Instruct the server to setup and queue the calculations.
            calculation_id = cached_data_id
            self._queued_calculations[calculation_id] = runner_data_model

        # Pass the ids of the submitted calculations back to the
        # client.
        encoded_job_ids = json.dumps(calculation_id).encode()
        length = pack_int(len(encoded_job_ids))

        await stream.write(length + encoded_job_ids)

        logging.info('Jobs ids sent to the client ({}): {}'.format(address, calculation_id))

        if should_launch:
            self.schedule_calculation(runner_data_model)

    async def handle_job_query(self, stream, address, message_length):
        """An asynchronous routine for handling the receiving and processing
        of job queries from a client

        Parameters
        ----------
        stream: IOStream
            An IO stream used to pass messages between the
            server and client.
        address: str
            The address from which the request came.
        message_length: int
            The length of the message being recieved.
        """

        logging.info('Received job query from {}'.format(address))

        encoded_ticket_id = await stream.read_bytes(message_length)
        ticket_id = encoded_ticket_id.decode()

        logging.info('Looking up ticket id {}'.format(ticket_id))

        response = None

        if (ticket_id not in self._queued_calculations and
            ticket_id not in self._finished_calculations):

            response = PropertyCalculatorException(directory='',
                                                   message='The {} ticket id was not found '
                                                           'on the server.'.format(ticket_id)).json()

        elif ticket_id in self._finished_calculations:
            response = self._finished_calculations[ticket_id].json()

        else:
            response = ''

        encoded_response = response.encode()
        length = pack_int(len(encoded_response))

        await stream.write(length + encoded_response)

        logging.info('Job results sent to the client {}: {}'.format(address, response))

    async def handle_stream(self, stream, address):
        """A routine to handle incoming requests from
        a property estimator TCP client.

        Notes
        -----

        This method is based on the StackOverflow response from
        A. Jesse Jiryu Davis: https://stackoverflow.com/a/40257248

        Parameters
        ----------
        stream: IOStream
            An IO stream used to pass messages between the
            server and client.
        address: str
            The address from which the request came.
        """
        logging.info("Incoming connection from {}".format(address))

        try:
            while True:

                # Receive a hello message with the message type.

                packed_message_type = await stream.read_bytes(4)
                message_type_int = unpack_int(packed_message_type)[0]

                packed_message_length = await stream.read_bytes(4)
                message_length = unpack_int(packed_message_length)[0]

                logging.info('Introductory packet recieved: {} {}'.format(message_type_int, message_length))

                message_type = None

                try:
                    message_type = PropertyEstimatorMessageTypes(message_type_int)
                    logging.info('Message type: {}'.format(message_type))

                except Exception as e:

                    logging.info('Bad message type recieved: {}'.format(e))

                    # Discard the unrecognised message.
                    if message_length > 0:
                        await stream.read_bytes(message_length)

                    continue

                if message_type is PropertyEstimatorMessageTypes.Submission:
                    await self.handle_job_submission(stream, address, message_length)
                elif message_type is PropertyEstimatorMessageTypes.Query:
                    await self.handle_job_query(stream, address, message_length)

        except StreamClosedError as e:

            # Handle client disconnections gracefully.
            logging.info("Lost connection to {}:{} : {}.".format(address, self._port, e))

    def _prepare_data_model(self, client_data_model):
        """Turns a client data model into a form more useful to
        the server side runner.

        Parameters
        ----------
        client_data_model: openforcefield.properties.estimator.PropertyEstimatorDataModel
            The client data model.

        Returns
        -------
        PropertyRunnerDataModel
            The server side data model.
        """
        parameter_set = ForceField([])
        parameter_set.__setstate__(client_data_model.parameter_set)

        parameter_set_id = self._storage_backend.has_force_field(parameter_set)

        if parameter_set_id is None:

            parameter_set_id = str(uuid.uuid4())
            self._storage_backend.store_force_field(parameter_set_id, parameter_set)

        calculation_id = str(uuid.uuid4())

        # Make sure we don't somehow generate the same uuid
        # twice (although this is very unlikely to ever happen).
        while calculation_id in self._queued_calculations:
            calculation_id = str(uuid.uuid4())

        runner_data = PropertyRunnerDataModel(id=calculation_id,
                                              queued_properties=client_data_model.properties,
                                              options=client_data_model.options,
                                              parameter_set_id=parameter_set_id)

        return runner_data

    def schedule_calculation(self, data_model):
        """Schedules the calculation of the given properties using the passed parameters.

        This method will recursively cascade through all allowed calculation
        layers or until all properties have been calculated.

        Parameters
        ----------
        data_model : PropertyRunnerDataModel
            The object containing instructions about which calculations
            should be performed.
        """

        if len(data_model.options.allowed_calculation_layers) == 0 or \
           len(data_model.queued_properties) == 0:

            self._queued_calculations.pop(data_model.id)
            self._finished_calculations[data_model.id] = data_model

            logging.info('Finished calculation {}'.format(data_model.id))
            return

        current_layer_type = data_model.options.allowed_calculation_layers.pop(0)

        if current_layer_type not in available_layers:

            # Kill all remaining properties if we reach an unsupported calculation layer.
            error_object = PropertyCalculatorException(directory='',
                                                       message='The {} calculation layer is not supported by '
                                                               'the server.'.format(current_layer_type))

            for queued_calculation in data_model.queued_properties:

                from openforcefield.properties import CalculationSource

                queued_calculation.source = CalculationSource(fidelity=current_layer_type,
                                                              provenance=error_object.json())

                data_model.calculated_properties.append(queued_calculation)

            data_model.options.allowed_calculation_layers.append(current_layer_type)
            data_model.queued_properties = []

            self.schedule_calculation(data_model)
            return

        logging.info('Launching calculation {} using the {} layer'.format(data_model.id,
                                                                          current_layer_type))

        layer_directory = path.join(self.working_directory, current_layer_type)

        if not path.isdir(layer_directory):
            os.makedirs(layer_directory)

        current_layer = available_layers[current_layer_type]

        current_layer.schedule_calculation(self._calculation_backend,
                                           self._storage_backend,
                                           layer_directory,
                                           data_model,
                                           self.schedule_calculation)

    def run_until_killed(self):
        """Starts the main (blocking) server IOLoop which will run until
        the user kills the process.
        """
        logging.info('Server listening at port {}'.format(self._port))

        try:
            IOLoop.current().start()
        except KeyboardInterrupt:
            self.stop()

    def run_until_complete(self):
        """Run the server loop until no more jobs are left in the queue.

        This is mainly made available for debug purposes for now.
        """
        def check_if_queue_empty():

            if len(self._queued_calculations) > 0:
                return

            logging.info('The queue is now empty - closing the server.')
            self.stop()

        check_callback = PeriodicCallback(check_if_queue_empty, 5000)
        self._periodic_loops.append(check_callback)

        check_callback.start()
        self.run_until_killed()

    def stop(self):
        """Stops the property calculation server and it's
        provided backend.
        """
        self._calculation_backend.stop()

        for periodic_loop in self._periodic_loops:
            periodic_loop.stop()

        IOLoop.current().stop()
