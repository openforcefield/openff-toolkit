# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Direct simulation layer.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>
"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import copy
import json
import logging
import os
import pickle
import traceback
import uuid
from os import path

import mdtraj

import openforcefield.properties.estimator.utils
from openforcefield.properties import PhysicalProperty, CalculationSource
from openforcefield.properties.estimator import CalculationSchema
from openforcefield.properties.estimator.layers.base import register_calculation_layer, PropertyCalculationLayer
from openforcefield.properties.estimator.storage import StoredSimulationData
from openforcefield.properties.estimator.utils import PropertyEstimatorException
from openforcefield.properties.estimator.workflow import protocols, groups
from openforcefield.utils import graph
from openforcefield.utils.serialization import PolymorphicDataType
from .base import CalculationLayerResult


# =============================================================================================
# Direct Calculation Classes
# =============================================================================================

class DirectCalculation:
    """Encapsulates and prepares the workflow needed to calculate a physical
    property by direct simulation methods.
    """

    def __init__(self, physical_property, force_field_path, schema, options):
        """
        Constructs a new DirectCalculation object.

        Parameters
        ----------
        physical_property: PhysicalProperty
            The protocol this node will execute.
        force_field_path: str
            The force field to use for this calculation.
        schema: CalculationSchema
            The schema to use to calculate this property.
        options: PropertyEstimatorOptions
            The options to run the calculation with.
        """
        self.physical_property = physical_property

        self.uuid = str(uuid.uuid4())

        self.protocols = {}

        self.starting_protocols = []
        self.dependants_graph = {}

        self.final_value_source = None
        self.final_uncertainty_source = None
        self.final_trajectory_source = None
        self.final_coordinate_source = None

        self.schema = CalculationSchema.parse_raw(schema.json())

        # Define a dictionary of accessible 'global' properties.
        self.global_properties = {
            "thermodynamic_state": physical_property.thermodynamic_state,
            "substance": physical_property.substance,
            "target_uncertainty": physical_property.uncertainty * options.relative_uncertainty,
            "force_field_path": force_field_path
        }

        # TODO: Nasty hack to turn a unitless quantity back into a unitless quantity after
        #       scalar multiplication.
        from simtk import unit

        if (isinstance(physical_property.uncertainty, unit.Quantity) and not
            isinstance(self.global_properties['target_uncertainty'], unit.Quantity)):

            self.global_properties['target_uncertainty'] = unit.Quantity(self.global_properties['target_uncertainty'],
                                                                         physical_property.uncertainty.unit)

        for protocol_name in self.schema.protocols:

            protocol_schema = self.schema.protocols[protocol_name]

            protocol = protocols.available_protocols[protocol_schema.type](protocol_schema.id)
            protocol.schema = protocol_schema

            # Try to set global properties on each of the protocols
            for input_path in protocol.required_inputs:

                input_value = protocol.get_value(input_path)

                if not isinstance(input_value, protocols.ProtocolPath):
                    continue

                if not input_value.is_global:
                    continue

                protocol.set_value(input_path, self.global_properties[input_value.property_name])

            if isinstance(protocol, groups.ConditionalGroup):

                for condition in protocol.conditions:

                    left_value = condition.left_hand_value

                    if isinstance(left_value, protocols.ProtocolPath) and left_value.is_global:
                        condition.left_hand_value = self.global_properties[left_value.property_name]

                    right_value = condition.right_hand_value
        
                    if isinstance(right_value, protocols.ProtocolPath) and right_value.is_global:
                        condition.right_hand_value = self.global_properties[right_value.property_name]

            protocol.set_uuid(self.uuid)
            self.protocols[protocol.id] = protocol

        self._build_dependants_graph()
        self.update_schema()

        self.final_value_source = copy.deepcopy(self.schema.final_value_source)
        self.final_value_source.append_uuid(self.uuid)

        self.final_uncertainty_source = copy.deepcopy(self.schema.final_uncertainty_source)
        self.final_uncertainty_source.append_uuid(self.uuid)

        self.final_trajectory_source = copy.deepcopy(self.schema.final_trajectory_source)
        self.final_trajectory_source.append_uuid(self.uuid)

        self.final_coordinate_source = copy.deepcopy(self.schema.final_coordinate_source)
        self.final_coordinate_source.append_uuid(self.uuid)

    def _build_dependants_graph(self):
        """Builds a dictionary of key value pairs where each key represents the id of a
        protocol to be executed in this calculation, and each value a list ids of protocols
        which must be ran after the protocol identified by the key.
        """

        for protocol_name in self.protocols:
            self.dependants_graph[protocol_name] = []

        for dependant_protocol_name in self.protocols:

            dependant_protocol = self.protocols[dependant_protocol_name]

            for dependency in dependant_protocol.dependencies:

                if dependency.is_global:
                    # Global inputs are outside the scope of the
                    # schema dependency graph.
                    continue

                if dependency.start_protocol == dependant_protocol_name and dependency.start_protocol:
                    # Don't add self to the dependency list.
                    continue

                # Only add a dependency on the protocol at the head of the path,
                # dependencies on the rest of protocols in the path is then implied.
                if dependant_protocol.id in self.dependants_graph[dependency.start_protocol]:
                    continue

                self.dependants_graph[dependency.start_protocol].append(dependant_protocol.id)

        self.starting_protocols = graph.find_root_nodes(self.dependants_graph)

    def replace_protocol(self, old_protocol, new_protocol):
        """Replaces an existing protocol with a new one, while
        updating all input and local references to point to the
        new protocol.

        The main use of this method is when merging multiple protocols
        into one.

        Parameters
        ----------
        old_protocol : protocols.BaseProtocol or str
            The protocol (or its id) to replace.
        new_protocol : protocols.BaseProtocol or str
            The new protocol (or its id) to use.
        """

        old_protocol_id = old_protocol
        new_protocol_id = new_protocol
        
        if isinstance(old_protocol, protocols.BaseProtocol):
            old_protocol_id = old_protocol.id
        if isinstance(new_protocol, protocols.BaseProtocol):
            new_protocol_id = new_protocol.id

        if new_protocol_id in self.protocols:
            raise ValueError('A protocol with the same id already exists in this calculation.')

        for protocol_id in self.protocols:

            protocol = self.protocols[protocol_id]
            protocol.replace_protocol(old_protocol_id, new_protocol_id)

        if old_protocol_id in self.protocols and isinstance(new_protocol, protocols.BaseProtocol):

            self.protocols.pop(old_protocol_id)
            self.protocols[new_protocol_id] = new_protocol

        for index, starting_id in enumerate(self.starting_protocols):

            if starting_id == old_protocol_id:
                starting_id = new_protocol_id

            self.starting_protocols[index] = starting_id

        for protocol_id in self.dependants_graph:

            for index, dependant_id in enumerate(self.dependants_graph[protocol_id]):

                if dependant_id == old_protocol_id:
                    dependant_id = new_protocol_id

                self.dependants_graph[protocol_id][index] = dependant_id

        if old_protocol_id in self.dependants_graph:
            self.dependants_graph[new_protocol_id] = self.dependants_graph.pop(old_protocol_id)

        self.final_value_source.replace_protocol(old_protocol_id, new_protocol_id)
        self.final_uncertainty_source.replace_protocol(old_protocol_id, new_protocol_id)
        self.final_trajectory_source.replace_protocol(old_protocol_id, new_protocol_id)
        self.final_coordinate_source.replace_protocol(old_protocol_id, new_protocol_id)

    def update_schema(self):

        for protocol_id in self.protocols:

            protocol = self.protocols[protocol_id]
            self.schema.protocols[protocol_id] = protocol.schema


class DirectCalculationGraph:
    """A hierarchical structure for storing and submitting the workflows
    which will calculate a set of physical properties..
    """

    def __init__(self, root_directory=None):
        """Constructs a new DirectCalculationGraph

        Parameters
        ----------
        root_directory: str, optional
            The root directory in which to store all outputs from
            this graph.
        """
        self._nodes_by_id = {}

        self._root_nodes = []
        self._root_directory = root_directory

        self._dependants_graph = {}

        self._calculations_to_run = {}

    def _insert_node(self, protocol_name, calculation, parent_node_name=None):
        """Recursively inserts a protocol node into the tree.

        Parameters
        ----------
        protocol_name : str
            The name of the protocol to insert.
        calculation : DirectCalculation
            The calculation being inserted.
        parent_node_name : str, optional
            The name of the new parent of the node to be inserted. If None,
            the protocol will be added as a new parent node.
        """

        if protocol_name in self._nodes_by_id:

            raise RuntimeError('A protocol with id ' + protocol_name + ' has already been inserted'
                                                                       ' into the graph.')

        nodes = self._root_nodes if parent_node_name is None else self._dependants_graph[parent_node_name]

        protocol_to_insert = calculation.protocols[protocol_name]
        existing_node = None

        # Start by checking to see if the starting node of the calculation graph is
        # already present in the full graph.
        for node_id in nodes:

            node = self._nodes_by_id[node_id]

            if not node.can_merge(protocol_to_insert):
                continue

            existing_node = node
            break

        if existing_node is not None:
            # Make a note that the existing node should be used in place
            # of this calculations version.

            merged_ids = existing_node.merge(protocol_to_insert)
            calculation.replace_protocol(protocol_to_insert, existing_node)

            for old_id, new_id in merged_ids.items():
                calculation.replace_protocol(old_id, new_id)

        else:

            parent_node = None if parent_node_name is None else self._nodes_by_id[parent_node_name]
            root_directory = self._root_directory if parent_node_name is None else None

            protocol_to_insert.directory = protocol_to_insert.id if parent_node is None \
                else path.join(parent_node.directory, protocol_to_insert.id)

            if root_directory is not None:

                protocol_to_insert.directory = path.join(root_directory,
                                                         protocol_to_insert.directory)

            # Add the protocol as a new node in the graph.
            self._nodes_by_id[protocol_name] = protocol_to_insert

            existing_node = self._nodes_by_id[protocol_name]
            self._dependants_graph[protocol_name] = []

            if parent_node_name is None:
                self._root_nodes.append(protocol_name)
            else:

                for node_id in calculation.dependants_graph:

                    if (protocol_name not in calculation.dependants_graph[node_id] or
                       node_id in self._dependants_graph[protocol_name]):

                        continue

                    self._dependants_graph[node_id].append(protocol_name)

        reduced_graph = copy.deepcopy(calculation.dependants_graph)
        graph.apply_transitive_reduction(reduced_graph)

        # Add all of the dependants to the existing node
        for dependant_name in reduced_graph[existing_node.id]:
            self._insert_node(dependant_name, calculation, existing_node.id)

    def add_calculation(self, calculation):
        """Insert a calculation into the calculation graph.

        Parameters
        ----------
        calculation : DirectCalculation
            The calculation to insert.
        """

        if calculation.uuid in self._calculations_to_run:

            # Quick sanity check.
            raise ValueError('A calculation with the same uuid ({}) is '
                             'trying to run twice.'.format(calculation.uuid))

        self._calculations_to_run[calculation.uuid] = calculation

        for starting_protocol_name in calculation.starting_protocols:
            self._insert_node(starting_protocol_name, calculation)

        self._calculations_to_run[calculation.uuid].update_schema()

    def submit(self, backend):
        """Submits the protocol graph to the backend of choice.

        Parameters
        ----------
        backend: PropertyEstimatorBackend
            The backend to launch the graph on.

        Returns
        -------
        list of Future:
            The futures of the submitted protocols.
        """
        submitted_futures = {}
        value_futures = []

        # Determine the ideal order in which to submit the
        # protocols.
        submission_order = graph.topological_sort(self._dependants_graph)
        # Build a dependency graph from the dependants graph so that
        # futures can be passed in the correct place.
        dependencies = graph.dependants_to_dependencies(self._dependants_graph)

        for node_id in submission_order:

            node = self._nodes_by_id[node_id]
            dependency_futures = []

            for dependency in dependencies[node_id]:
                dependency_futures.append(submitted_futures[dependency])

            # Pull out any 'global' properties.
            global_properties = {}

            for dependency in node.dependencies:

                if not dependency.is_global:
                    continue

                global_properties[dependency.property_name] = node.get_value(dependency)

            # Do a quick sanity check to make sure they've all been set already.
            if len(global_properties) > 0:
                raise ValueError('The global_properties array should be empty by this point.')

            submitted_futures[node_id] = backend.submit_task(DirectCalculationGraph._execute_protocol,
                                                             node.directory,
                                                             node.schema,
                                                             *dependency_futures)

        for calculation_id in self._calculations_to_run:

            calculation = self._calculations_to_run[calculation_id]
            calculation.update_schema()

            # TODO: Fill in any extra required provenance.
            provenance = {}

            for protocol_id in calculation.protocols:

                protocol = calculation.protocols[protocol_id]
                provenance[protocol_id] = PolymorphicDataType(protocol.schema)

            source = CalculationSource(fidelity=SimulationLayer.__name__,
                                       provenance=provenance)

            calculation.physical_property.source = source

            value_node_id = calculation.final_value_source.start_protocol
            uncertainty_node_id = calculation.final_uncertainty_source.start_protocol
            trajectory_node_id = calculation.final_trajectory_source.start_protocol
            coordinate_node_id = calculation.final_coordinate_source.start_protocol

            final_futures = [
                submitted_futures[value_node_id],
                submitted_futures[uncertainty_node_id],
                submitted_futures[trajectory_node_id],
                submitted_futures[coordinate_node_id]
            ]

            # Gather the values and uncertainties of each property being calculated.
            value_futures.append(backend.submit_task(DirectCalculationGraph._gather_results,
                                                     calculation.physical_property,
                                                     calculation.final_value_source,
                                                     calculation.final_uncertainty_source,
                                                     calculation.final_coordinate_source,
                                                     calculation.final_trajectory_source,
                                                     *final_futures))

        return value_futures

    @staticmethod
    def _execute_protocol(directory, protocol_schema, *parent_outputs, available_resources, **kwargs):
        """Executes a protocol defined by ``protocol_schema``, and with
        inputs sets via the global scope and from previously executed protocols.

        Parameters
        ----------
        protocol_schema: protocols.ProtocolSchema
            The schema defining the protocol to execute.
        parent_outputs: tuple of Any
            The results of previous protocol executions.

        Returns
        -------
        str, dict of str and Any
            Returns a tuple of the id of the executed protocol, and a dictionary
            which contains the outputs of the executed protocol.
        """

        # Store the results of the relevant previous protocols in a handy dictionary.
        # If one of the results is a failure, propagate it up the chain!
        parent_outputs_by_path = {}

        for parent_id, parent_output in parent_outputs:

            if isinstance(parent_output, openforcefield.properties.estimator.utils.PropertyEstimatorException):
                return protocol_schema.id, parent_output

            for output_path, output_value in parent_output.items():

                property_name, protocol_ids = protocols.ProtocolPath.to_components(output_path)

                if len(protocol_ids) == 0 or (len(protocol_ids) > 0 and protocol_ids[0] != parent_id):
                    protocol_ids.insert(0, parent_id)

                final_path = protocols.ProtocolPath(property_name, *protocol_ids)
                parent_outputs_by_path[final_path] = output_value

        # Recreate the protocol on the backend to bypass the need for static methods
        # and awkward args and kwargs syntax.
        protocol = protocols.available_protocols[protocol_schema.type](protocol_schema.id)
        protocol.schema = protocol_schema

        if not path.isdir(directory):
            os.makedirs(directory)

        for input_path in protocol.required_inputs:

            target_path = protocol.get_value(input_path)

            if not isinstance(target_path, protocols.ProtocolPath):
                continue

            if (target_path.start_protocol == input_path.start_protocol or
                target_path.start_protocol == protocol.id):

                continue

            input_value = parent_outputs_by_path[target_path]
            protocol.set_value(input_path, input_value)

        if isinstance(protocol, groups.ConditionalGroup):

            for condition in protocol.conditions:

                left_value = condition.left_hand_value

                if (isinstance(left_value, protocols.ProtocolPath) and 
                    left_value.start_protocol is not None and
                    left_value.start_protocol != protocol.id):
                    
                    condition.left_hand_value = parent_outputs_by_path[left_value]

                right_value = condition.right_hand_value

                if (isinstance(right_value, protocols.ProtocolPath) and
                        right_value.start_protocol is not None and
                        right_value.start_protocol != protocol.id):

                    condition.right_hand_value = parent_outputs_by_path[right_value]

        try:
            output_dictionary = protocol.execute(directory, available_resources)
        except Exception as e:
            # Except the unexpected...
            formatted_exception = traceback.format_exception(None, e, e.__traceback__)

            return protocol.id, PropertyEstimatorException(directory=directory,
                                                           message='An unhandled exception occurred: '
                                                                   '{}'.format(formatted_exception))

        return protocol.id, output_dictionary

    @staticmethod
    def _gather_results(property_to_return, value_reference, uncertainty_reference,
                        coordinate_reference, trajectory_reference, *protocol_results, **kwargs):
        """Gather the value and uncertainty calculated from the submission graph
        and store them in the property to return.

        Parameters
        ----------
        value_result: dict of string and Any
            The result dictionary of the protocol which calculated the value of the property.
        value_reference: ProtocolPath
            A reference to which property in the output dictionary is the actual value.
        uncertainty_result: dict of string and Any
            The result dictionary of the protocol which calculated the uncertainty in the value.
        uncertainty_reference: ProtocolPath
            A reference to which property in the output dictionary is the actual uncertainty.
        property_to_return: PhysicalProperty
            The property to which the value and uncertainty belong.

        Returns
        -------
        CalculationLayerResult
            The result of attempting to estimate this property by direct simulation.
        """

        return_object = CalculationLayerResult()
        return_object.property_id = property_to_return.id

        results_by_id = {}

        for protocol_id, protocol_results in protocol_results:

            # Make sure none of the protocols failed and we actually have a value
            # and uncertainty.
            if isinstance(protocol_results, openforcefield.properties.estimator.utils.PropertyEstimatorException):

                return_object.calculation_error = protocol_results
                return return_object

            for output_path, output_value in protocol_results.items():

                property_name, protocol_ids = protocols.ProtocolPath.to_components(output_path)

                if len(protocol_ids) == 0 or (len(protocol_ids) > 0 and protocol_ids[0] != protocol_id):
                    protocol_ids.insert(0, protocol_id)

                final_path = protocols.ProtocolPath(property_name, *protocol_ids)
                results_by_id[final_path] = output_value

        property_to_return.value = results_by_id[value_reference]
        property_to_return.uncertainty = results_by_id[uncertainty_reference]

        return_object.calculated_property = property_to_return

        data_to_store = StoredSimulationData()

        data_to_store.substance = property_to_return.substance
        data_to_store.thermodynamic_state = property_to_return.thermodynamic_state

        data_to_store.provenance = property_to_return.source

        data_to_store.source_calculation_id = property_to_return.id

        coordinate_path = results_by_id[coordinate_reference]
        trajectory_path = results_by_id[trajectory_reference]

        data_to_store.trajectory_data = mdtraj.load_dcd(trajectory_path, top=coordinate_path)

        return_object.data_to_store = data_to_store

        return return_object


# =============================================================================================
# Simulation Layer
# =============================================================================================

@register_calculation_layer()
class SimulationLayer(PropertyCalculationLayer):
    """A calculation layer which aims to calculate physical properties
    directly from molecular simulation.

    .. warning :: This class is experimental and should not be used in a production environment.
    """

    @staticmethod
    def _build_calculation_graph(working_directory, properties, force_field_path, options):
        """ Construct a graph of the protocols needed to calculate a set of properties.

        Parameters
        ----------
        working_directory: str
            The local directory in which to store all local,
            temporary calculation data from this graph.
        properties : list of PhysicalProperty
            The properties to attempt to compute.
        force_field_path : str
            The path to the force field parameters to use in the calculation.
        options: PropertyEstimatorOptions
            The options to run the calculations with.
        """
        calculation_graph = DirectCalculationGraph(working_directory)

        for property_to_calculate in properties:

            property_type = type(property_to_calculate).__name__

            if property_type not in options.calculation_schemas:

                logging.warning('The property calculator does not support {} '
                                'calculations.'.format(property_type))

                continue

            schema = options.calculation_schemas[property_type]

            calculation = DirectCalculation(property_to_calculate,
                                            force_field_path,
                                            schema,
                                            options)

            calculation_graph.add_calculation(calculation)

        return calculation_graph

    @staticmethod
    def schedule_calculation(calculation_backend, storage_backend, layer_directory,
                             data_model, callback, synchronous=False):

        # Store a temporary copy of the force field for protocols to easily access.
        force_field = storage_backend.retrieve_force_field(data_model.parameter_set_id)
        force_field_path = path.join(layer_directory, 'force_field_{}'.format(data_model.parameter_set_id))

        with open(force_field_path, 'wb') as file:
            pickle.dump(force_field, file)

        calculation_graph = SimulationLayer._build_calculation_graph(layer_directory,
                                                                     data_model.queued_properties,
                                                                     force_field_path,
                                                                     data_model.options)

        simulation_futures = calculation_graph.submit(calculation_backend)

        PropertyCalculationLayer._await_results(calculation_backend, storage_backend, layer_directory,
                                                data_model, callback, simulation_futures, synchronous)
