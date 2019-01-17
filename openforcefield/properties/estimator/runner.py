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

import logging
import uuid
import json
import pickle
import struct
import os
import hashlib

from os import path

from pydantic import BaseModel
from typing import Dict, List, Any

from simtk import unit

from tornado.ioloop import IOLoop, PeriodicCallback
from tornado.iostream import IOStream, StreamClosedError
from tornado.tcpserver import TCPServer

from openforcefield.utils.serialization import serialize_quantity

from openforcefield.typing.engines.smirnoff import ForceField

from openforcefield.properties import PhysicalProperty
from openforcefield.properties.estimator.client import PropertyEstimatorDataModel, PropertyEstimatorOptions
from openforcefield.properties.estimator.layers import available_layers
from openforcefield.properties.estimator.components.protocols import ProtocolPath

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
    """

    id: str

    queued_properties: List[PhysicalProperty] = []
    calculated_properties: List[PhysicalProperty] = []

    options: PropertyEstimatorOptions = None

    parameter_set_path: str = None

    stored_data: Dict[str, Any] = {}

    class Config:

        # A dirty hack to allow simtk.unit.Quantities...
        # TODO: Should really investigate QCElemental as an alternative.
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

    Examples
    --------
    Setting up a general server instance using a dask LocalCluster backend:

    >>> # Create the backend which will be responsible for distributing the calculations
    >>> from openforcefield.properties.estimator.backends import DaskLocalClusterBackend
    >>> backend = DaskLocalClusterBackend(1, 1)
    >>>
    >>> # Create the server to which all calculations will be submitted
    >>> from openforcefield.properties.estimator.runner import PropertyCalculationRunner
    >>> property_server = PropertyCalculationRunner(backend)
    >>>
    >>> # Instruct the server to listen for incoming submissions
    >>> property_server.run_until_complete()
    """

    def __init__(self, backend, port=8000):
        """Constructs a new PropertyCalculationRunner object.

        Parameters
        ----------
        backend: PropertyEstimatorBackend
            The backend to use for executing calculations.
        port: int
            The port one which to listen for incoming client
            requests.
        """
        self._backend = backend
        self._port = port

        self._queued_calculations = {}
        self._finished_calculations = {}

        # Keep a track of the parameters used for each of the
        # calculations and assign each set a server unique id.
        self._force_field_hashes = {}

        self._force_field_directory = 'cached_force_fields'
        self._force_field_hash_name = 'hash_keys'

        if not path.isdir(self._force_field_directory):
            os.makedirs(self._force_field_directory)
        else:
            self._import_cached_force_fields()

        self._periodic_loops = []

        super().__init__()

        # self.listen(self._port)
        self.bind(self._port)
        self.start(1)

        backend.start()

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

                # Read the incoming request from the server. The first four bytes
                # of the response should be the length of the message being sent.
                header = await stream.read_bytes(4)
                length = unpack_int(header)[0]

                # Decode the client submission json.
                encoded_json = await stream.read_bytes(length)
                json_model = encoded_json.decode()

                # TODO: Add exeception handling so the server can gracefully reject bad json.
                client_data_model = PropertyEstimatorDataModel.parse_raw(json_model)

                logging.info('Received job request from {}'.format(address))

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

        parameter_set_id = self._is_force_field_in_cache(parameter_set)

        if parameter_set_id is None:

            parameter_set_id = str(uuid.uuid4())
            self._cache_force_field(parameter_set_id, parameter_set)

        parameter_set_path = path.join(self._force_field_directory, parameter_set_id)

        calculation_id = str(uuid.uuid4())

        # Make sure we don't somehow generate the same uuid
        # twice (although this is very unlikely to ever happen).
        while calculation_id in self._queued_calculations:
            calculation_id = str(uuid.uuid4())

        runner_data = PropertyRunnerDataModel(id=calculation_id,
                                              queued_properties=client_data_model.properties,
                                              options=client_data_model.options,
                                              parameter_set_path=parameter_set_path)

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
            # TODO: Implement graceful error handling.
            return

        logging.info('Launching calculation {} using the {} layer'.format(data_model.id,
                                                                          current_layer_type))

        current_layer = available_layers[current_layer_type]
        current_layer.schedule_calculation(self._backend, data_model, {}, self.schedule_calculation)

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
        self._backend.stop()

        for periodic_loop in self._periodic_loops:
            periodic_loop.stop()

        IOLoop.current().stop()

    def _is_force_field_in_cache(self, force_field):
        """Checks whether the force field has been used
        as part of a previous calculation.

        TODO: Load in the hashed file

        Parameters
        ----------
        force_field: ForceField
            The force field to check for.

        Returns
        -------
        str, None
            None if the force field has not been cached, otherwise
            the unique id of the cached force field.
        """

        force_field_pickle = pickle.dumps(force_field.__getstate__())
        hash_string = hashlib.sha256(force_field_pickle).hexdigest()

        for unique_id in self._force_field_hashes:

            existing_hash = self._force_field_hashes[unique_id].replace('\n', '')

            if hash_string != existing_hash:
                continue

            existing_path = path.join(self._force_field_directory, unique_id)

            if not path.isfile(existing_path):
                # For some reason the force field got deleted..
                continue

            return unique_id

        return None

    def _cache_force_field(self, unique_id, force_field):
        """Store the force field in the cached force field
        directory.

        TODO: In future will most likely need to also serialize the
              SMIRNOFF engine version for correct comparision

        Parameters
        ----------
        unique_id: str
            The unique id assigned to the force field.
        force_field: ForceField
            The force field to cache.
        """

        force_field_pickle = pickle.dumps(force_field.__getstate__())
        hash_string = hashlib.sha256(force_field_pickle).hexdigest()

        with open(path.join(self._force_field_directory, unique_id), 'wb') as file:
            pickle.dump(force_field.__getstate__(), file)

        hash_file_path = path.join(self._force_field_directory,
                                   self._force_field_hash_name)

        with open(hash_file_path, 'a+') as file:
            file.write('{} {}\n'.format(unique_id, hash_string))

        self._force_field_hashes[unique_id] = hash_string

    def _import_cached_force_fields(self):
        """Loads all of the the force fields which are cached in the
        force field directory.
        """

        hash_file_path = path.join(self._force_field_directory,
                                   self._force_field_hash_name)

        if not path.isfile(hash_file_path):
            return

        with open(hash_file_path, 'r') as file:

            hash_lines = file.readlines()

            for hash_line in hash_lines:

                unique_id = hash_line.split(' ')[0]
                hash_key = hash_line.split(' ')[1]

                self._force_field_hashes[unique_id] = hash_key
