#!/usr/bin/env python

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
import struct

from pydantic import BaseModel
from typing import Dict, List, Any

from simtk import unit

from tornado.ioloop import IOLoop
from tornado.iostream import IOStream, StreamClosedError
from tornado.tcpserver import TCPServer

from openforcefield.utils.serialization import serialize_quantity

from openforcefield.properties.estimator.client import PropertyEstimatorDataModel, PropertyEstimatorOptions
from openforcefield.properties import PhysicalProperty
from openforcefield.properties.estimator.layers import available_layers

# Needed for server-client communication.
int_struct = struct.Struct("<i")

unpack_int = int_struct.unpack
pack_int = int_struct.pack


# =============================================================================================
# Property Runner
# =============================================================================================

class PropertyRunnerDataModel(BaseModel):

    id: str

    queued_properties: List[PhysicalProperty] = []
    calculated_properties: List[PhysicalProperty] = []

    options: PropertyEstimatorOptions = None

    parameter_set: Dict[int, str] = None

    stored_data: Dict[str, Any] = {}

    class Config:

        # A dirty hack to allow simtk.unit.Quantities...
        # TODO: Should really investigate QCElemental as an alternative.
        arbitrary_types_allowed = True

        json_encoders = {
            unit.Quantity: lambda v: serialize_quantity(v),
        }


class PropertyCalculationRunner(TCPServer):
    """The class responsible for coordinating all calculations to be ran using
    the property estimator, in addition to deciding at which fidelity a property
    will be calculated.

    It acts as a server, which receives submitted jobs from clients
    launched via the property estimator.
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
                data_model = PropertyEstimatorDataModel.parse_raw(json_model)

                logging.info('Received job request from {}'.format(address))

                calculation_id = None
                should_launch = True

                # Make sure this job is not already in the queue.
                # for existing_id in self._queued_calculations:
                #
                #     if data_model.json() == self._queued_calculations[existing_id].json():
                #
                #         calculation_id = existing_id
                #         should_launch = False
                #
                #         break

                if calculation_id is None:
                    # Instruct the server to setup and queue the calculations.
                    calculation_id = str(uuid.uuid4())

                    # Make sure we don't somehow generate the same uuid
                    # twice (although this is very unlikely to ever happen).
                    while calculation_id in self._queued_calculations:
                        calculation_id = str(uuid.uuid4())

                    self._queued_calculations[calculation_id] = data_model

                # Pass the ids of the submitted calculations back to the
                # client.
                encoded_job_ids = json.dumps(calculation_id).encode()
                length = pack_int(len(encoded_job_ids))

                await stream.write(length + encoded_job_ids)

                logging.info('Jobs ids sent to the client ({}): {}'.format(address, calculation_id))

                if should_launch:

                    runner_data = PropertyRunnerDataModel(id=calculation_id,
                                                          queued_properties=data_model.properties,
                                                          options=data_model.options,
                                                          parameter_set=data_model.parameter_set)

                    self.launch_calculation(runner_data)

        except StreamClosedError as e:

            # Handle client disconnections gracefully.
            logging.info("Lost connection to {}:{} : {}.".format(address, self._port, e))

    def launch_calculation(self, data_model):
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

            # All properties have been calculated, or have they...?
            return

        initial_layer_type = data_model.options.allowed_calculation_layers.pop(0)

        if initial_layer_type not in available_layers:
            # TODO Graceful error handling.
            return

        logging.info('Launching calculation {} at {} fidelity'.format(data_model.id, initial_layer_type))

        initial_layer = available_layers[initial_layer_type]
        initial_layer.perform_calculation(self._backend, data_model, {}, self.launch_calculation)

    def start_listening_loop(self):
        """Starts the main (blocking) server IOLoop.
        """
        logging.info('Server listening at port {}'.format(self._port))
        IOLoop.current().start()
