# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Property estimator client side API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import json
import logging
import struct
from time import sleep
from typing import Dict, List

from pydantic import BaseModel
from simtk import unit
from tornado.ioloop import IOLoop
from tornado.iostream import StreamClosedError
from tornado.tcpclient import TCPClient

from openforcefield.properties import PhysicalProperty
from openforcefield.properties.estimator import CalculationSchema
from openforcefield.properties.estimator.layers import SurrogateLayer, ReweightingLayer, SimulationLayer
from openforcefield.properties.estimator.workflow.protocols import ProtocolPath
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils.serialization import serialize_quantity, PolymorphicDataType
from .utils import PropertyEstimatorMessageTypes, PropertyEstimatorException

int_struct = struct.Struct("<i")

unpack_int = int_struct.unpack
pack_int = int_struct.pack


# =============================================================================================
# Registration Decorators
# =============================================================================================

def register_estimable_property():
    """A decorator which registers a property as being estimable
    by the property estimator.

    Notes
    -----
    The property must implement a static get_calculation_template method
    which returns the calculation schema to follow.
    """

    def decorator(cls):

        if cls.__name__ in PropertyEstimator.registered_properties:
            raise ValueError('The {} property is already registered.'.format(cls.__name__))

        PropertyEstimator.registered_properties[cls.__name__] = cls
        return cls

    return decorator


# =============================================================================================
# Property Estimator
# =============================================================================================

class PropertyEstimatorOptions(BaseModel):
    """Represents additional options that can be passed to the
    property estimator backend.

    Attributes
    ----------
    allowed_calculation_layers: list of str
        A list of allowed calculation layers. The order of the layers in the list is the order
        that the calculator will attempt to execute the layers in.
    calculation_schemas: Dict[str, CalculationSchema]
        A dictionary of the CalculationSchema which will be used to calculate any properties.
        The dictionary key represents the type of property the schema will calculate. The
        dictionary will be automatically populated with defaults if no entries are added.
    relative_uncertainty: float, default = 1.0
        Control the desired uncertainty of any calculated properties, relative to measured ones.
    allow_protocol_merging: bool, default = True
        If true, allows individual, identical steps in a property calculation to be merged.
    """
    allowed_calculation_layers: List[str] = [
        SurrogateLayer.__name__,
        ReweightingLayer.__name__,
        SimulationLayer.__name__
    ]

    calculation_schemas: Dict[str, CalculationSchema] = {}

    relative_uncertainty: float = 1.0
    allow_protocol_merging: bool = True

    class Config:

        arbitrary_types_allowed = True

        json_encoders = {
            unit.Quantity: lambda v: serialize_quantity(v),
            ProtocolPath: lambda v: v.full_path,
            PolymorphicDataType: lambda value: PolymorphicDataType.serialize(value)
        }


class PropertyEstimatorSubmission(BaseModel):
    """Represents a set of properties to be calculated by the estimator,
    the parameters which will be used to calculate them, and options about
    how the properties will be calculated.

    Attributes
    ----------
    properties: list of PhysicalProperty
        The list of physical properties to calculate.
    options: PropertyEstimatorOptions
        The options used to calculate the properties.
    parameter_set: dict of str and int
        The force field parameters used during the calculations. These should be
        obtained by calling .__getstate__() on a `ForceField` object.
    """
    properties: List[PhysicalProperty] = []
    options: PropertyEstimatorOptions = None

    parameter_set: Dict[int, str] = None

    class Config:

        arbitrary_types_allowed = True

        json_encoders = {
            unit.Quantity: lambda v: serialize_quantity(v),
            ProtocolPath: lambda v: v.full_path,
            PolymorphicDataType: lambda value: PolymorphicDataType.serialize(value)
        }


class PropertyEstimatorResult(BaseModel):
    """Represents a set of properties to be calculated by the estimator,
    the parameters which will be used to calculate them, and options about
    how the properties will be calculated.

    Attributes
    ----------
    properties: list of PhysicalProperty
        The list of physical properties to calculate.
    options: PropertyEstimatorOptions
        The options used to calculate the properties.
    parameter_set: dict of str and int
        The force field parameters used during the calculations. These should be
        obtained by calling .__getstate__() on a `ForceField` object.
    """
    id: str

    calculated_properties: Dict[str, PhysicalProperty] = {}
    unsuccessful_properties: Dict[str, PropertyEstimatorException] = {}

    parameter_set_id: str = None

    class Config:
        arbitrary_types_allowed = True

        json_encoders = {
            unit.Quantity: lambda v: serialize_quantity(v),
            ProtocolPath: lambda v: v.full_path,
            PolymorphicDataType: lambda value: PolymorphicDataType.serialize(value)
        }


class PropertyEstimator(object):
    """The main object responsible for requesting a set of properties
    be calculated by the low-level property calculation backend.

    The property estimator acts as a client object, which is able to connect
    and submit calculations to a remote (or even a local) property estimator server.

    Examples
    --------

    Setting up the server instance:

    >>> # Create the backend which will be responsible for distributing the calculations
    >>> from openforcefield.properties.estimator.backends import DaskLocalClusterBackend, BackendResources
    >>> calculation_backend = DaskLocalClusterBackend(1, 1, BackendResources(1, 0))
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

    If the server and client initialisation occur in the same function call (e.g. the same
    main function), then `property_server.run_until_complete()` must be called after
    all computations have been submitted.

    Setting up the client instance:

    >>> # Load in the data set of properties which will be used for comparisons
    >>> from openforcefield.properties.datasets import ThermoMLDataSet
    >>> data_set = ThermoMLDataSet.from_doi_list('10.1016/j.jct.2016.10.001')
    >>> # Filter the dataset to only include densities measured between 130-260 K
    >>> from openforcefield.properties import Density
    >>>
    >>> data_set.filter_by_properties(types=[Density.__name__])
    >>> data_set.filter_by_temperature(min_temperature=130*unit.kelvin, max_temperature=260*unit.kelvin)
    >>>
    >>> # Load initial parameters
    >>> from openforcefield.typing.engines.smirnoff import ForceField
    >>> parameters = ForceField('smirnoff99Frosst.offxml')
    >>>
    >>> # Create a property estimator
    >>> from openforcefield.properties.estimator import PropertyEstimator
    >>> property_estimator = PropertyEstimator()

    Submit the calculations using the default submission options:

    >>> options = PropertyEstimatorOptions()
    >>> ticket_ids = property_estimator.submit_computations(data_set, parameters, options)

    The layers which will be used to calculate the properties can be controlled
    using the submission options:

    >>> options = PropertyEstimatorOptions(allowed_calculation_layers = [ReweightingLayer.__name__,
    >>>                                                                  SimulationLayer.__name__])

    As can the uncertainty tolerance:

    >>> options = PropertyEstimatorOptions(relative_uncertainty = 0.1)

    Once a job has been submitted to the server, it's status can be tracked by calling

    >>> property_estimator.query_computation_state(ticket_ids)
    """

    registered_properties = {}

    def __init__(self, server_address='localhost', port=8000):
        """Constructs a new PropertyEstimator object.

        Parameters
        ----------
        server_address: str
            The address of the calculation server.
        port: int
            The port that the server is listening on.
        """

        self._server_address = server_address

        if server_address is None:

            raise ValueError('The address of the server which will run'
                             'these calculations must be given.')

        self._port = port
        self._tcp_client = TCPClient()

    def submit_computations(self, data_set, parameter_set, options=None):
        """
        Submit the property and parameter set for calculation.

        Parameters
        ----------
        data_set : PropertyDataSet
            The set of properties to attempt to compute.
        parameter_set : ForceField
            The OpenFF parameter set to use for the calculations.
        options : PropertyEstimatorOptions, optional
            A set of additional calculation options.

        Returns
        -------
        list of str:
            A list unique ids which can be used to retrieve the submitted calculations
            when they have finished running.
        """

        if data_set is None or parameter_set is None:

            raise ValueError('Both a data set and parameter set must be '
                             'present to compute physical properties.')

        if options is None:
            options = PropertyEstimatorOptions()

        properties_list = []

        for substance_tag in data_set.properties:

            for physical_property in data_set.properties[substance_tag]:

                properties_list.append(physical_property)

                type_name = type(physical_property).__name__

                if type_name not in PropertyEstimator.registered_properties:

                    raise ValueError('The property estimator does not support {} '
                                     'properties.'.format(type_name))

                if type_name in options.calculation_schemas:
                    continue

                property_type = PropertyEstimator.registered_properties[type_name]()

                options.calculation_schemas[type_name] = \
                    property_type.get_default_calculation_schema()

        for property_schema_name in options.calculation_schemas:

            options.calculation_schemas[property_schema_name].validate_interfaces()

            for protocol_schema_name in options.calculation_schemas[property_schema_name].protocols:

                protocol_schema = options.calculation_schemas[
                    property_schema_name].protocols[protocol_schema_name]

                protocol_schema.inputs['.allow_merging'] = PolymorphicDataType(options.allow_protocol_merging)

        submission = PropertyEstimatorSubmission(properties=properties_list,
                                                 parameter_set=parameter_set.__getstate__(),
                                                 options=options)

        # For now just do a blocking submit to the server.
        ticket_ids = IOLoop.current().run_sync(lambda: self._send_calculations_to_server(submission))

        return ticket_ids

    def query_computation_state(self, ticket_id):
        """
        Submit the property and parameter set for calculation.

        Parameters
        ----------
        ticket_id: str
            The id of the calculation to query.

        Returns
        -------
        PropertyEstimatorResult or PropertyCalculatorException, optional:
           The result of the submitted job. Returns None if the calculation has
           not yet completed.
        """

        # For now just do a blocking submit to the server.
        result = IOLoop.current().run_sync(lambda: self._send_query_server(ticket_id))

        return result

    def wait_for_result(self, ticket_id, interval=10):
        """
        Synchronously wait for the result of a calculation

        Parameters
        ----------
        ticket_id: str
            The id of the calculation to wait for.
        interval: int
            The time interval (seconds) between checking if the calculation has finished.

        Returns
        -------
        PropertyEstimatorResult or PropertyCalculatorException, optional:
           The result of the submitted job. Returns None if the calculation has
           not yet completed.
        """
        assert interval >= 1

        response = None
        should_run = True

        while should_run:

            sleep(interval)

            response = IOLoop.current().run_sync(lambda: self._send_query_server(ticket_id))

            if response is None:
                continue

            logging.info('The server has returned a response.')
            should_run = False

        return response

    async def _send_calculations_to_server(self, submission):
        """Attempts to connect to the calculation server, and
        submit the requested calculations.

        Notes
        -----

        This method is based on the StackOverflow response from
        A. Jesse Jiryu Davis: https://stackoverflow.com/a/40257248

        Parameters
        ----------
        submission: PropertyEstimatorSubmission
            The jobs to submit.

        Returns
        -------
        str, optional:
           The id which the server has assigned the submitted calculations.
           This can be used to query the server for when the calculation
           has completed.

           Returns None if the calculation could not be submitted.
        """
        ticket_id = None

        try:

            # Attempt to establish a connection to the server.
            logging.info("Attempting Connection to {}:{}".format(self._server_address, self._port))
            stream = await self._tcp_client.connect(self._server_address, self._port)
            logging.info("Connected to {}:{}".format(self._server_address, self._port))

            stream.set_nodelay(True)

            # Encode the submission json into an encoded
            # packet ready to submit to the server. The
            # Length of the packet is encoded in the first
            # four bytes.
            message_type = pack_int(PropertyEstimatorMessageTypes.Submission)

            encoded_json = submission.json().encode()
            length = pack_int(len(encoded_json))

            await stream.write(message_type + length + encoded_json)

            logging.info("Sent calculations to {}:{}. Waiting for a response from"
                         " the server...".format(self._server_address, self._port))

            # Wait for confirmation that the server has submitted
            # the jobs. The first four bytes of the response should
            # be the length of the message being sent.
            header = await stream.read_bytes(4)
            length = unpack_int(header)[0]

            # Decode the response from the server. If everything
            # went well, this should be a list of ids of the submitted
            # calculations.
            encoded_json = await stream.read_bytes(length)
            ticket_id = json.loads(encoded_json.decode())

            logging.info('Received job id from server: {}'.format(ticket_id))
            stream.close()
            self._tcp_client.close()

        except StreamClosedError as e:

            # Handle no connections to the server gracefully.
            logging.info("Error connecting to {}:{} : {}. Please ensure the server is running and"
                         "that the server address / port is correct.".format(self._server_address, self._port, e))

        # Return the ids of the submitted jobs.
        return ticket_id

    async def _send_query_server(self, ticket_id):
        """Attempts to connect to the calculation server, and
        submit the requested calculations.

        Notes
        -----

        This method is based on the StackOverflow response from
        A. Jesse Jiryu Davis: https://stackoverflow.com/a/40257248

        Parameters
        ----------
        ticket_id: str
            The id of the job to query.

        Returns
        -------
        str, optional:
           The status of the submitted job.
           Returns None if the calculation has not yet completed.
        """
        server_response = None

        try:

            # Attempt to establish a connection to the server.
            logging.info("Attempting Connection to {}:{}".format(self._server_address, self._port))
            stream = await self._tcp_client.connect(self._server_address, self._port)
            logging.info("Connected to {}:{}".format(self._server_address, self._port))

            stream.set_nodelay(True)

            # Encode the ticket id into the message.
            message_type = pack_int(PropertyEstimatorMessageTypes.Query)

            encoded_ticket_id = ticket_id.encode()
            length = pack_int(len(encoded_ticket_id))

            await stream.write(message_type + length + encoded_ticket_id)

            logging.info("Querying the server {}:{}...".format(self._server_address, self._port))

            # Wait for the server response.
            header = await stream.read_bytes(4)
            length = unpack_int(header)[0]

            # Decode the response from the server. If everything
            # went well, this should be the finished calculation.
            if length > 0:

                encoded_json = await stream.read_bytes(length)
                server_response = encoded_json.decode()

            logging.info('Received response from server of length {}'.format(length))

            stream.close()
            self._tcp_client.close()

        except StreamClosedError as e:

            # Handle no connections to the server gracefully.
            logging.info("Error connecting to {}:{} : {}. Please ensure the server is running and"
                         "that the server address / port is correct.".format(self._server_address, self._port, e))

        if server_response is not None:

            try:
                server_response = PropertyEstimatorResult.parse_raw(server_response)
            except:
                server_response = PropertyEstimatorException.parse_raw(server_response)

        # Return the ids of the submitted jobs.
        return server_response
