# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Direct simulation layer.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

.. todo::

    * Make all protocol execute methods static.
    * Create a pydantic data model to pass input / outputs between protocols.
    * Gather up the final results of the calculation.
"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import os
import copy
import logging
import uuid

from os import path

from openforcefield.utils import graph

from openforcefield.properties import PhysicalProperty
from openforcefield.properties.estimator.components import protocols
from openforcefield.properties.estimator.layers.base import register_calculation_layer, PropertyCalculationLayer

from openforcefield.typing.engines.smirnoff import ForceField

from .base import return_args


# =============================================================================================
# Direct Calculation Classes
# =============================================================================================

class DirectCalculation:
    """Defines the property to calculate and the calculation
    workflow needed to calculate it.
    """
    def __init__(self, physical_property, force_field, schema):
        """
        Constructs a new DirectCalculation object.

        Parameters
        ----------
        physical_property: PhysicalProperty
            The protocol this node will execute.
        force_field: openforcefield.typing.engines.smirnoff.ForceField
            The force field to use for this calculation.
        schema: CalculationSchema
            The schema to use to calculate this property.
        """
        self.physical_property = physical_property
        self.uuid = str(uuid.uuid4())

        self.protocols = {}

        self.starting_protocols = []
        self.dependants_graph = {}

        self.final_value_reference = None
        self.final_uncertainty_reference = None

        # Define a dictionary of accessible 'global' properties.
        global_properties = {
            "thermodynamic_state": physical_property.thermodynamic_state,
            "substance": physical_property.substance,
            "uncertainty": physical_property.uncertainty,
            "force_field": force_field
        }

        # A helper to map the template protocol ids to
        # the new, uuid appended ones.
        id_map = {}

        for protocol_name in schema.protocols:

            protocol_schema = schema.protocols[protocol_name]

            protocol = protocols.available_protocols[protocol_schema.type]()
            protocol.schema = protocol_schema

            # Try to set global properties on each of the protocols
            protocol.set_global_properties(global_properties)

            protocol.set_uuid(self.uuid)
            id_map[protocol_name] = protocol.id

            self.protocols[protocol.id] = protocol

        self._build_dependants_graph()
        self._apply_groups(schema)

        self.final_value_reference = copy.deepcopy(schema.final_value_reference)
        self.final_value_reference.set_uuid(self.uuid)

        self.final_uncertainty_reference = copy.deepcopy(schema.final_uncertainty_reference)
        self.final_uncertainty_reference.set_uuid(self.uuid)

    def _build_dependants_graph(self):
        """Builds a dictionary of key value pairs where each key represents the id of a
        protocol to be executed in this calculation, and each value a list ids of protocols
        which must be ran after the protocol identified by the key.
        """

        for protocol_name in self.protocols:
            self.dependants_graph[protocol_name] = []

        for dependant_protocol_name in self.protocols:

            dependant_protocol = self.protocols[dependant_protocol_name]

            for input_reference in dependant_protocol.input_references:

                if input_reference.output_protocol_id == 'global':
                    # Global inputs are outside the scope of the
                    # schema dependency graph.
                    continue

                if dependant_protocol.id in self.dependants_graph[input_reference.output_protocol_id]:
                    continue

                self.dependants_graph[input_reference.output_protocol_id].append(dependant_protocol.id)

        self.starting_protocols = graph.find_root_nodes(self.dependants_graph)

        # Remove all implicit dependencies from the dependants graph.
        graph.apply_transitive_reduction(self.dependants_graph)

    def _apply_groups(self, schema):
        """Groups protocols together into a set of user defined groups."""

        if len(schema.groups) == 0:
            # Nothing to do here.
            return

        # TODO: Implement groups.
        return
        #
        # for group_id in self.groups:
        #
        #     group = self.groups[group_id]
        #
        #     # Remove all grouped protocols from the protocol list.
        #     for grouped_protocol_id in group.protocols:
        #         self.protocols.pop(grouped_protocol_id)
        #
        #     # Point the other protocols to the groups rather than
        #     # the removed protocols
        #     for grouped_protocol_id in group.protocols:
        #
        #         for protocol_id in self.protocols:
        #
        #             protocol = self.protocols[protocol_id]
        #             protocol.replace_protocol(grouped_protocol_id, group.id)
        #
        #             for input_reference in protocol.input_references:
        #
        #                 if input_reference.output_protocol_id != group.id:
        #                     continue
        #
        #                 input_reference.grouped_protocol_id = grouped_protocol_id
        #
        #         if self.final_value_reference.output_protocol_id == grouped_protocol_id:
        #             self.final_value_reference.output_protocol_id = group.id
        #             self.final_value_reference.grouped_protocol_id = grouped_protocol_id
        #
        #         if self.final_uncertainty_reference.output_protocol_id == grouped_protocol_id:
        #             self.final_uncertainty_reference.output_protocol_id = group.id
        #             self.final_uncertainty_reference.grouped_protocol_id = grouped_protocol_id
        #
        #     # Add the group in their place.
        #     self.protocols[group.id] = group

    def replace_protocol(self, old_protocol, new_protocol):
        """Replaces an existing protocol with a new one, while
        updating all input and local references to point to the
        new protocol.

        The main use of this method is when merging multiple protocols
        into one.

        Parameters
        ----------
        old_protocol : protocols.BaseProtocol
            The protocol to replace.
        new_protocol : protocols.BaseProtocol
            The new protocol to use.
        """

        if new_protocol.id in self.protocols:
            raise ValueError('A protocol with the same id already exists in this calculation.')

        for protocol_id in self.protocols:

            protocol = self.protocols[protocol_id]
            protocol.replace_protocol(old_protocol.id, new_protocol.id)

        if old_protocol.id in self.protocols:

            self.protocols.pop(old_protocol.id)
            self.protocols[new_protocol] = new_protocol

        for index, starting_id in enumerate(self.starting_protocols):

            if starting_id == old_protocol.id:
                starting_id = new_protocol.id

            self.starting_protocols[index] = starting_id

        for protocol_id in self.dependants_graph:

            for index, dependant_id in enumerate(self.dependants_graph[protocol_id]):

                if dependant_id == old_protocol.id:
                    dependant_id = new_protocol.id

                self.dependants_graph[protocol_id][index] = dependant_id

        if old_protocol.id in self.dependants_graph:
            self.dependants_graph[new_protocol.id] = self.dependants_graph.pop(old_protocol.id)

        self.final_value_reference.replace_protocol(old_protocol.id, new_protocol.id)
        self.final_uncertainty_reference.replace_protocol(old_protocol.id, new_protocol.id)


class DirectCalculationGraph:
    """A hierarchical structure for storing all protocols that need to be executed.
    """

    class CalculationNode:
        """A node in the calculation graph.

        Each node represents a protocol to be executed.
        """

        def __init__(self, protocol, parent_node, root_directory=None):
            """Constructs a new CalculationNode

            Parameters
            ----------
            protocol: Protocol
                The protocol this node will execute.
            parent_node: CalculationNode
                The parent of this node.
            root_directory: str, optional
                The root directory in which to store all outputs from
                this graph.
            """
            self.protocol = protocol

            self.directory = protocol.id if parent_node is None else \
                path.join(parent_node.directory, protocol.id)

            if root_directory is not None:
                self.directory = path.join(root_directory, self.directory)

            # self.status = self.Status.Queueing

        @property
        def id(self):
            """str: Returns the id of the protocol to execute."""
            return self.protocol.id

        def execute(self, *input_protocols):
            """Execute the protocol's protocol.

            Parameters
            ----------
            input_protocols : list(CalculationNode)
                The input protocols which this node depends on for input.
            """

            input_protocols_by_id = {}

            for input_protocol in input_protocols:
                input_protocols_by_id[input_protocol.id] = input_protocol

            if not path.isdir(self.directory):
                os.makedirs(self.directory)

            for input_reference in self.protocol.input_references:

                if input_reference.output_protocol_id == 'global':
                    continue

                output_protocol = input_protocols_by_id[input_reference.output_protocol_id]

                input_value = output_protocol.get_output_value(input_reference)
                self.protocol.set_input_value(input_reference, input_value)

            self.protocol.execute(self.directory)
            return self.protocol

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

    @property
    def calculations_to_run(self):
        """list(CalculationNode): A list of the calculations to be executed."""
        return self._calculations_to_run

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

        if protocol_name in self._dependants_graph:

            raise RuntimeError('A protocol with id ' + protocol_name + ' has already been inserted'
                                                                       ' into the graph.')

        nodes = self._root_nodes if parent_node_name is None else self._dependants_graph[parent_node_name]

        protocol_to_insert = calculation.protocols[protocol_name]
        existing_node = None

        # Start by checking to see if the starting node of the calculation graph is
        # already present in the full graph.
        for node_id in nodes:

            node = self._nodes_by_id[node_id]

            if not node.protocol.can_merge(protocol_to_insert):
                continue

            existing_node = node
            break

        if existing_node is not None:
            # Make a note that the existing node should be used in place
            # of this calculations version.

            existing_node.protocol.merge(protocol_to_insert)
            calculation.replace_protocol(protocol_to_insert, existing_node.protocol)

        else:

            parent_node = None if parent_node_name is None else self._nodes_by_id[parent_node_name]
            root_directory = self._root_directory if parent_node_name is None else None

            # Add the protocol as a new node in the graph.
            self._nodes_by_id[protocol_name] = self.CalculationNode(protocol_to_insert,
                                                                    parent_node,
                                                                    root_directory)

            existing_node = self._nodes_by_id[protocol_name]
            self._dependants_graph[protocol_name] = []

            nodes.append(protocol_name)

        # Add all of the dependants to the existing node
        for dependant_name in calculation.dependants_graph[existing_node.id]:
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
            raise ValueError('A calculation with the same uuid (' +
                             calculation.uuid + ') is trying to run twice.')

        self._calculations_to_run[calculation.uuid] = calculation

        for starting_protocol_name in calculation.starting_protocols:
            self._insert_node(starting_protocol_name, calculation)

    def submit(self, backend):
        """Executes the protocol graph, employing a dask backend.

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

        submission_order = graph.topological_sort(self._dependants_graph)
        dependencies = {}

        for node_id in self._nodes_by_id:

            node = self._nodes_by_id[node_id]

            dependencies[node_id] = []

            for input_reference in node.protocol.input_references:

                if input_reference.output_protocol_id == 'global':
                    continue

                if input_reference.output_protocol_id not in dependencies[node_id]:
                    dependencies[node_id].append(input_reference.output_protocol_id)

        for node_id in submission_order:

            node = self._nodes_by_id[node_id]

            dependency_futures = []

            for dependency in dependencies[node_id]:
                dependency_futures.append(submitted_futures[dependency])

            submitted_futures[node_id] = backend.submit_task(node.execute, *dependency_futures)

        for calculation_id in self._calculations_to_run:

            calculation = self._calculations_to_run[calculation_id]

            value_node_id = calculation.final_value_reference.output_protocol_id
            uncertainty_node_id = calculation.final_uncertainty_reference.output_protocol_id

            value_futures.append(backend.submit_task(return_args,
                                                     submitted_futures[value_node_id],
                                                     submitted_futures[uncertainty_node_id]))

        return value_futures


# =============================================================================================
# Simulation Layer
# =============================================================================================

@register_calculation_layer()
class SimulationLayer(PropertyCalculationLayer):

    @staticmethod
    def _build_calculation_graph(properties, force_field, schemas):
        """ Construct a graph of the protocols needed to calculate a set of properties.

        Parameters
        ----------
        properties : list of PhysicalProperty
            The properties to attempt to compute.
        force_field : ForceField
            The force field parameters to use in the calculation.
        schemas : dict of str and CalculationSchema
            A list of the schemas to use when performing the calculations.
        """
        calculation_graph = DirectCalculationGraph('property-data')

        for property_to_calculate in properties:

            property_type = type(property_to_calculate).__name__

            if property_type not in schemas:

                logging.warning('The property calculator does not support {} '
                                'calculations.'.format(property_type))

                continue

            schema = schemas[property_type]

            calculation = DirectCalculation(property_to_calculate,
                                            force_field,
                                            schema)

            calculation_graph.add_calculation(calculation)

        return calculation_graph

    @staticmethod
    def perform_calculation(backend, data_model, existing_data, callback, synchronous=False):

        parameter_set = ForceField([])
        parameter_set.__setstate__(data_model.parameter_set)

        calculation_graph = SimulationLayer._build_calculation_graph(data_model.queued_properties,
                                                                     parameter_set,
                                                                     data_model.options.calculation_schemas)

        simulation_futures = calculation_graph.submit(backend)

        PropertyCalculationLayer._await_results(backend, data_model, callback, simulation_futures, synchronous)

    @staticmethod
    def finalise_results(physical_property):
        """A placeholder method that would be used to attempt
        to reweight previous calculations to yield the desired
        property.
        """

        # For now the return tuple indicates that the reweighting
        # was not sufficiently accurate to estimate the property (False)
        # and simply returns the property back to be passed to the next layer.

        # """Extract any properties calculated by a calculation graph.
        #
        # Parameters
        # ----------
        # calculation_graph : DirectCalculationGraph
        #     An already executed calculation graph.
        # """
        #
        # for calculation_uuid in calculation_graph.calculations_to_run:
        #
        #     calculation = calculation_graph.calculations_to_run[calculation_uuid]
        #
        #     value_protocol = calculation.protocols[calculation.final_value_reference.output_protocol_id]
        #     uncertainty_protocol = calculation.protocols[calculation.final_uncertainty_reference.output_protocol_id]
        #
        #     value = value_protocol.get_output_value(calculation.final_value_reference)
        #     uncertainty = uncertainty_protocol.get_output_value(calculation.final_uncertainty_reference)
        #
        #     calculated_property = PhysicalProperty()
        #
        #     calculated_property.substance = calculation.physical_property.substance
        #     calculated_property.thermodynamic_state = calculation.physical_property.thermodynamic_state
        #
        #     calculated_property.phase = PropertyPhase.Liquid
        #
        #     calculated_property.source = CalculationSource(CalculationFidelity.DirectSimulation)
        #
        #     calculated_property.value = value
        #     calculated_property.uncertainty = uncertainty
        #
        #     self.calculated_properties[str(calculated_property.substance)].append(calculated_property)

        return False, physical_property
