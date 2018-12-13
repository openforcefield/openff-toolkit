#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Property calculator 'server' side API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

TODO: Need to version as discussed with Jeff.

"""


# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import os
import copy
import logging
import multiprocessing
import math

import uuid

from os import path
from xml.etree import ElementTree

from openforcefield.utils import graph
from openforcefield.utils.exceptions import XmlNodeMissingException

from openforcefield.properties.estimator.components import protocols, groups
from openforcefield.properties.estimator.components.protocols import ProtocolInputReference

from openforcefield.properties import PropertyType, PropertyPhase, CalculationFidelity, \
                                      PhysicalProperty, CalculatedPhysicalProperty

from openforcefield.properties.datasets.property_dataset import CalculatedPropertySet

from dask.threaded import get as dask_threaded_get
from dask import distributed


# =============================================================================================
# Calculation Classes
# =============================================================================================

class DirectCalculationTemplate:
    """Defines the set of protocols required to calculate a certain property.
    """
    def __init__(self):

        self.property_type = PropertyType.Undefined
        self.id = None

        self.protocols = {}
        self.groups = {}

        self.final_value_reference = None
        self.final_uncertainty_reference = None

        # A list of protocols which have zero or only global inputs.
        # These will be the first protocols to be executed.
        self.starting_protocols = []

        # A key-value pair dictionary of dependants, where
        # each key represents a protocol and each value is a list
        # of the protocols which take said protocol as an input.
        self.dependants_graph = {}

    @classmethod
    def from_xml(cls, xml_node, existing_templates):
        """ Imports a set of protocols from an xml definition.

        Parameters
        ----------
        xml_node : xml.etree.Element
            The element containing the xml to read from.
        existing_templates : list(DirectPropertyCalculationTemplate)
            A list of already loaded calculation templates

        Returns
        ----------
        DirectCalculationTemplate
            The calculated template created from the xml node.
        """
        return_value = cls()

        # Determine which property this calculation will yeild.
        property_node = xml_node.find('property')

        if property_node is None:
            raise XmlNodeMissingException('property')

        return_value.property_type = PropertyType[property_node.text]

        # Find the unique id of this calculation.
        id_node = xml_node.find('id')

        if id_node is None:
            raise XmlNodeMissingException('id')

        return_value.id = id_node.text

        # Make sure the id is really unique.
        if return_value.id in existing_templates:

            raise Exception('A calculation with id {} has been '
                            'defined more than once.'.format(return_value.id))

        return_value.id = id_node.text

        # Load in the protocol definitions.
        cls._parse_protocols_node(existing_templates, return_value, xml_node)
        # Load in the group definitions.
        cls._parse_groups_node(return_value, xml_node)

        cls.validate_inputs(return_value)

        # Find the output value of this calculation.
        output_value_node = xml_node.find('output-value')

        if output_value_node is not None:

            output_value_split = output_value_node.text.split(':')

            if len(output_value_split) != 2:
                raise Exception('The format of the output-value node '
                                'should be protocol_id:property_name')

            protocol_id = output_value_split[0]
            property_name = output_value_split[1]

            if protocol_id not in return_value.protocols:

                raise Exception('The {} calculation does not contain '
                                'a protocol with an id {}'.format(return_value.id, protocol_id))

            if property_name not in return_value.protocols[protocol_id].provided_outputs:

                raise Exception('The {} protocol does not output a '
                                'property named {}'.format(protocol_id, property_name))

            return_value.final_value_reference = ProtocolInputReference('',
                                                                        protocol_id,
                                                                        property_name)

        # Find the output uncertainty of this calculation.
        output_uncertainty_node = xml_node.find('output-uncertainty')

        if output_uncertainty_node is not None:

            output_uncertainty_split = output_uncertainty_node.text.split(':')

            if len(output_uncertainty_split) != 2:
                raise Exception('The format of the output-uncertainty node '
                                'should be protocol_id:property_name')

            protocol_id = output_uncertainty_split[0]
            property_name = output_uncertainty_split[1]

            if protocol_id not in return_value.protocols:

                raise Exception('The {} calculation does not contain '
                                'a protocol with an id {}'.format(return_value.id, protocol_id))

            if property_name not in return_value.protocols[protocol_id].provided_outputs:

                raise Exception('The {} protocol does not output a '
                                'property named {}'.format(protocol_id, property_name))

            return_value.final_uncertainty_reference = ProtocolInputReference('',
                                                                              protocol_id,
                                                                              property_name)

        return_value.apply_groups()
        return_value.build_dependants_graph()

        if not graph.is_acyclic(return_value.dependants_graph):

            raise Exception('The {} contains a cycle. Only acyclic workflows '
                            'are supported'.format(return_value.id))

        return return_value

    @classmethod
    def _parse_protocols_node(cls, existing_templates, calculation_template, xml_node):
        """

        Parameters
        ----------
        existing_templates : list(DirectPropertyCalculationTemplate)
            A list of already loaded calculation templates.
        calculation_template: DirectCalculationTemplate
            The calculation to add the parsed protocols to.
        xml_node : xml.etree.Element
            The element containing the xml to read from.

        Returns
        -------

        """

        protocol_nodes = xml_node.find('protocols')

        if protocol_nodes is None:
            raise XmlNodeMissingException('protocols')

        # Load in all of the protocols.
        for protocol_node in protocol_nodes:

            # Handle the case where we import the protocols
            # from an already defined calculation.
            if protocol_node.tag == 'import':

                if protocol_node.text not in existing_templates:

                    raise RuntimeError('The {} template tries to inherit from a non-existent '
                                       'template: {}'.format(calculation_template.id, protocol_node.text))

                calculation_template.protocols.update(copy.deepcopy(existing_templates[protocol_node.text].protocols))

                continue

            # Try to automatically find which class corresponds to the xml tag.
            protocol_class = getattr(protocols, protocol_node.tag)

            if protocol_class is None:

                raise RuntimeError('The {} calculation contains an undefined protocol of '
                                   'type: {}'.format(calculation_template.id, protocol_node.tag))

            # Make sure the found class is actually a protocol.
            if isinstance(protocol_class, protocols.BaseProtocol):

                raise RuntimeError('{} is not a child of Protocol'.format(protocol_node.tag))

            protocol = protocol_class.from_xml(protocol_node)

            if protocol.id in calculation_template.protocols:

                raise Exception('The {} calculation defines two protocols with the same '
                                'id: {}'.format(calculation_template.id, protocol.id))

            calculation_template.protocols[protocol.id] = protocol

    @classmethod
    def _parse_groups_node(cls, calculation_template, xml_node):
        """Used to parse the xml <groups> node.

        Parameters
        ----------
        calculation_template: DirectCalculationTemplate
            The template to add the parsed groups to.
        xml_node
            The xml nodes that contains the groups node to parse.
        """
        protocol_group_nodes = xml_node.find('groups')

        # Load in all of the protocol groups.
        for protocol_group_node in [] if protocol_group_nodes is None else protocol_group_nodes:

            # Try to automatically find which class corresponds to the xml tag.
            protocol_group_class = getattr(groups, protocol_group_node.tag)

            if protocol_group_class is None:

                raise RuntimeError('The {} calculation contains an undefined protocol '
                                   'group of type: {}'.format(calculation_template.id, protocol_group_node.tag))

            # Make sure the found class is actually a protocol_group.
            if isinstance(protocol_group_class, groups.ProtocolGroup):

                raise RuntimeError('{} is not a child of ProtocolGroup'.format(protocol_group_node.tag))

            protocol_group = protocol_group_class.from_xml(protocol_group_node, calculation_template.protocols)

            if protocol_group.id in calculation_template.groups:

                raise Exception('The {} calculation defines two protocol groups with the '
                                'same id: {}'.format(calculation_template.id, protocol_group.id))

            if protocol_group.id in calculation_template.protocols:

                raise Exception('The {} calculation defines a protocol group {} with the same id as an existing '
                                'protocol.'.format(calculation_template.id, protocol_group.id))

            # Make sure protocols haven't been added to more than one group.
            for existing_group_id in calculation_template.groups:

                existing_group = calculation_template.groups[existing_group_id]

                for existing_group_protocol_id in existing_group:

                    if existing_group_protocol_id not in protocol_group:
                        continue

                    raise Exception('The {} protocol has been added to multiple groups in the '
                                    '{} calculation'.format(calculation_template.id, existing_group_protocol_id))

            calculation_template.groups[protocol_group.id] = protocol_group

    @staticmethod
    def validate_inputs(template):
        """Validates the flow of the data between protocols, ensuring
        that inputs and outputs correctly match up.

        Parameters
        ----------
        template: DirectCalculationTemplate
            The template to validate.
        """

        for protocol_name in template.protocols:

            protocol = template.protocols[protocol_name]

            for input_reference in protocol.input_references:

                if input_reference.output_protocol_id == 'global':
                    # We handle global inputs separately
                    continue

                # Make sure the other protocol whose output we are interested
                # in actually exists.
                if input_reference.output_protocol_id not in template.protocols:

                    raise Exception('The {} protocol of the {} template tries to take input from a non-existent '
                                    'protocol: {}'.format(protocol.id, template.id, input_reference.protocol_id))

                other_protocol = template.protocols[input_reference.output_protocol_id]

                # Make sure the other protocol definitely has the requested output.
                if input_reference.output_property_name not in other_protocol.provided_outputs:

                    raise Exception('The {} protocol does not provide an {} output.'.format(
                        other_protocol.id, input_reference.output_property_name))

    def build_dependants_graph(self):
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
                    # template dependency graph.
                    continue

                if dependant_protocol.id in self.dependants_graph[input_reference.output_protocol_id]:
                    continue

                self.dependants_graph[input_reference.output_protocol_id].append(dependant_protocol.id)

        self.starting_protocols = graph.find_root_nodes(self.dependants_graph)

        # Remove all implicit dependencies from the dependants graph.
        graph.apply_transitive_reduction(self.dependants_graph)

    def apply_groups(self):
        """Groups protocols together into a set of user defined groups."""

        if len(self.groups) == 0:
            # Nothing to do here.
            return

        for group_id in self.groups:

            group = self.groups[group_id]

            # Remove all grouped protocols from the protocol list.
            for grouped_protocol_id in group.protocols:
                self.protocols.pop(grouped_protocol_id)

            # Point the other protocols to the groups rather than
            # the removed protocols
            for grouped_protocol_id in group.protocols:

                for protocol_id in self.protocols:

                    protocol = self.protocols[protocol_id]
                    protocol.replace_protocol(grouped_protocol_id, group.id)

                    for input_reference in protocol.input_references:

                        if input_reference.output_protocol_id != group.id:
                            continue

                        input_reference.grouped_protocol_id = grouped_protocol_id

                if self.final_value_reference.output_protocol_id == grouped_protocol_id:

                    self.final_value_reference.output_protocol_id = group.id
                    self.final_value_reference.grouped_protocol_id = grouped_protocol_id

                if self.final_uncertainty_reference.output_protocol_id == grouped_protocol_id:

                    self.final_uncertainty_reference.output_protocol_id = group.id
                    self.final_uncertainty_reference.grouped_protocol_id = grouped_protocol_id

            # Add the group in their place.
            self.protocols[group.id] = group


class DirectCalculation:
    """Defines the property to calculate and the calculation
    workflow needed to calculate it.
    """
    def __init__(self, physical_property, force_field, template):
        """
        Constructs a new DirectCalculation object.

        Parameters
        ----------
        physical_property: PhysicalProperty
            The protocol this node will execute.
        force_field: openforcefield.typing.engines.smirnoff.ForceField
            The force field to use for this calculation.
        template: DirectCalculationTemplate, optional
            The template on which this calculation is based.
        """
        self.physical_property = physical_property
        self.uuid = str(uuid.uuid4())

        self.property_type = template.property_type

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

        for protocol_name in template.protocols:

            protocol = copy.deepcopy(template.protocols[protocol_name])

            # Try to set global properties on each of the protocols
            protocol.set_global_properties(global_properties)

            protocol.set_uuid(self.uuid)
            id_map[protocol_name] = protocol.id

            self.protocols[protocol.id] = protocol

        # Copy the starting protocols, making sure to use the uuid appended names.
        for starting_protocol in template.starting_protocols:
            self.starting_protocols.append(id_map[starting_protocol])

        # Copy the dependants graph, making sure to use the uuid appended names.
        for dependant_name in template.dependants_graph:

            self.dependants_graph[id_map[dependant_name]] = []

            for input_name in template.dependants_graph[dependant_name]:
                self.dependants_graph[id_map[dependant_name]].append(id_map[input_name])

        self.final_value_reference = copy.deepcopy(template.final_value_reference)
        self.final_value_reference.set_uuid(self.uuid)

        self.final_uncertainty_reference = copy.deepcopy(template.final_uncertainty_reference)
        self.final_uncertainty_reference.set_uuid(self.uuid)

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
        # TODO: This logic is nowhere near robust enough yet.
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

    def execute(self, number_of_workers, threads_per_worker):
        """Executes the protocol graph, employing a dask backend."""

        dask_graph = {}
        leaf_nodes = []

        dependencies = {}

        for node_id in self._nodes_by_id:

            node = self._nodes_by_id[node_id]
            dependencies[node_id] = []

            if len(self._dependants_graph[node_id]) == 0:
                leaf_nodes.append(node_id)

            input_tuple = (node.execute, )

            for input_reference in node.protocol.input_references:

                if input_reference.output_protocol_id == 'global':
                    continue

                if input_reference.output_protocol_id not in dependencies[node_id]:
                    dependencies[node_id].append(input_reference.output_protocol_id)

                input_tuple = input_tuple + (input_reference.output_protocol_id, )

            dask_graph[node_id] = input_tuple

        print('Started')

        value = dask_threaded_get(dask_graph, leaf_nodes, num_workers=number_of_workers)

        # cluster = distributed.LocalCluster(number_of_workers, threads_per_worker, processes=False)
        # client = distributed.Client(cluster, processes=False)
        #
        # value = client.get(dask_graph, leaf_nodes)

        print('Finished', value)


# =============================================================================================
# Property Estimator
# =============================================================================================

class PropertyCalculationRunner:
    """The class responsible for deciding at which fidelity a property
    is calculated, as well as launching the calculations themselves.

    Parameters
    ----------
    worker_threads : int
        The number of threads to run the calculations on.
    threads_per_simulation : int
        The number of threads to use per multi-threaded simulation.
    """

    def __init__(self, worker_threads=1, threads_per_simulation=None):

        self.calculated_properties = {}

        self._worker_threads = worker_threads
        self._threads_per_simulation = threads_per_simulation

        maximum_threads = multiprocessing.cpu_count()

        if self._threads_per_simulation is None:
            self._threads_per_simulation = math.floor(maximum_threads / worker_threads)

        if worker_threads * self._threads_per_simulation > maximum_threads:

            raise ValueError('The total number of requested threads (' +
                             str(worker_threads * self._threads_per_simulation) +
                             ') must be less than the available ' + str(maximum_threads))

        logging.info(str(self._threads_per_simulation) + ' threads will be used per'
                                                         ' openMM simulation')

        os.environ["OPENMM_NUM_THREADS"] = str(threads_per_simulation)

        self.calculator_templates = {}

        self.load_calculator_templates()

    def load_calculator_templates(self):
        """Loads the property calculation templates from the default xml file.
        """
        templates = {}

        with open('../data/properties/templates.xml') as templates_file:

            xml = templates_file.read()
            root_node = ElementTree.fromstring(xml)

            for template_node in root_node:

                calculator_template = DirectCalculationTemplate.from_xml(template_node,
                                                                         templates)

                templates[calculator_template.id] = calculator_template

        for template_id in templates:

            template = templates[template_id]

            if template.property_type is PropertyType.Undefined:
                continue

            self.calculator_templates[template.property_type] = template

    def run(self, measured_properties, parameter_set):
        """Schedules the calculation of the given properties using the passed parameters.

        Parameters
        ----------
        measured_properties : PhysicalPropertyDataSet
            The properties to attempt to compute.
        parameter_set : PhysicalPropertyDataSet
            The force field parameters used to perform the calculations.
        """

        # Make a copy of the properties array so we
        # can easily safely changes to it.
        properties_to_measure = {}

        for substance_hash in measured_properties:

            if substance_hash not in properties_to_measure:
                properties_to_measure[substance_hash] = []

            properties_to_measure[substance_hash].extend(measured_properties[substance_hash])

        # Reset the calculation arrays
        self.calculated_properties = {}

        # Set up the array to track any pending or
        # finished calculations
        for substance_hash in properties_to_measure:
                self.calculated_properties[substance_hash] = []

        # First see if any of these properties could be calculated
        # from a surrogate model. If they can, remove them the list
        # of remaining properties to be calculated.
        self.attempt_surrogate_extrapolation(parameter_set, properties_to_measure)

        # Exit early if everything is covered by the surrogate
        if len(properties_to_measure) == 0:
            return self.calculated_properties

        # Do a similar thing to check whether reweighting could
        # be employed here for any of the properties.
        self.attempt_simulation_reweighting(parameter_set, properties_to_measure)

        # Exit early if everything is covered by reweighting
        if len(properties_to_measure) == 0:
            return self.calculated_properties

        # Next, for each property we need to figure out which protocols
        # are required to calculate the properties in properties_by_mixture
        #
        # Take for example the calculation of the density of water at 298K. One would need
        # a set of protocols to:
        #
        #   1) generate the coordinates of a water box (GenerateMixtureProtocol).
        #   2) do any necessary EM / NVT equilibration (EnergyMinimisationProtocol) + (NVTProtocol).
        #   3) perform the actual NPT production run (NPTProtocol).
        #
        # By splitting up the calculation in this way, one could easily
        # set up the calculation of even complex properties by piecing together
        # simple component protocols.
        #
        # With these protocols defined, the protocols to be run for
        # each mixture type need to be compared to one another.
        #
        # Let's assume in the above example we also wish to calculate the
        # heat capacity of water at 298K. This calculation would again
        # need a set of protocols to:
        #
        #   1) generate the coordinates of a water box (GenerateMixtureProtocol).
        #   2) do any necessary EM / NVT equilibration (EnergyMinimisationProtocol) + (NVTProtocol).
        #   3) perform the actual NPT production run (NPTProtocol).
        #
        # clearly we don't want to calculate these two properties separately as
        # this would lead to unnecessary repetition.
        #
        # Instead, for each property to be calculated for water, all the protocols
        # to be run would be compared to see which are duplicates - these duplicates
        # would then be collapsed down to the bare minimum.
        #
        # This is done by storing all protocols to be run in a graph like structure,
        # where chain branches indicate a divergence in protocol
        #
        # More than likely most graphs will diverge at the analysis stage.
        calculation_graph = self.build_calculation_graph(properties_to_measure, parameter_set)

        # Once a list of unique list of protocols to be executed has been constructed
        # the protocols can be fired.
        self.execute_calculation_graph(calculation_graph)

        # Finally, extract the calculated properties.
        self.extract_calculated_properties(calculation_graph)

        # The results from these calculations would then be stored in some way so that
        # the data required for reweighting is retained.
        #
        # At this point, two things could happen:
        #
        #   1) The simulation results could be fed into the surrogate model
        #   2) The results would be passed back to the 'client'

        property_set = CalculatedPropertySet(self.calculated_properties, parameter_set)
        return property_set

    def attempt_surrogate_extrapolation(self, parameter_set, measured_properties):
        """Attempts to calculate properties from a surrogate model"""

        # This loop would be distributed over a number of threads /
        # nodes depending on the available architecture.
        for substance_hash in measured_properties:

            for measured_property in measured_properties[substance_hash]:

                surrogate_result = self.perform_surrogate_extrapolation(measured_property, parameter_set)

                if surrogate_result is None:
                    continue

                substance_tag = str(measured_property.substance)

                self.calculated_properties[substance_tag].append(surrogate_result)
                measured_properties[substance_hash].remove(measured_property)

    def perform_surrogate_extrapolation(self, measured_property, parameter_set):
        """A placeholder method that would be used to spawn the surrogate
        model backend.
        """
        return None

    def attempt_simulation_reweighting(self, parameter_set, measured_properties):
        """Attempts to calculate properties by reweighting existing data"""

        for substance_hash in measured_properties:

            for measured_property in measured_properties[substance_hash]:

                reweighting_result = self.perform_reweighting(measured_property, parameter_set)

                if reweighting_result is None:
                    continue

                substance_tag = str(measured_property.substance)

                self.calculated_properties[substance_tag].append(reweighting_result)
                measured_properties[substance_hash].remove(measured_property)

    def perform_reweighting(self, measured_property, parameter_set):
        """A placeholder method that would be used to spawn the property
        reweighting backend.
        """
        return None

    def build_calculation_graph(self, properties_by_mixture, force_field):
        """ Construct a graph of the protocols needed to calculate a set of properties.

        Parameters
        ----------
        properties_by_mixture : PhysicalPropertyDataSet
            The properties to attempt to compute.
        """
        calculation_graph = DirectCalculationGraph('property-data')

        for mixture in properties_by_mixture:

            for property_to_calculate in properties_by_mixture[mixture]:

                if property_to_calculate.type not in self.calculator_templates:

                    logging.warning('The property calculator does not support ' +
                                    property_to_calculate.type + ' calculations.')

                    continue

                calculation = DirectCalculation(property_to_calculate,
                                                force_field,
                                                self.calculator_templates[property_to_calculate.type])

                calculation_graph.add_calculation(calculation)

        return calculation_graph

    def execute_calculation_graph(self, calculation_graph):
        """Execute a graph of calculation protocols using a given parameter set.

        Parameters
        ----------
        calculation_graph : DirectCalculationGraph
            The calculation graph to execute.
        """
        calculation_graph.execute(self._worker_threads, self._threads_per_simulation)

        # task_queue = Queue()
        #
        # logging.info('Spawning ' + str(self._worker_threads) + ' worker threads.')
        #
        # # Set up the pool of threads.
        # for i in range(self._worker_threads):
        #
        #     worker = Thread(target=DirectCalculationGraph.worker_thread, args=(task_queue,))
        #     worker.setDaemon(True)
        #     worker.start()
        #
        # calculation_graph.execute(task_queue)
        #
        # task_queue.join()

    def extract_calculated_properties(self, calculation_graph):
        """Extract any properties calculated by a calculation graph.

        Parameters
        ----------
        calculation_graph : DirectCalculationGraph
            An already executed calculation graph.
        """

        for calculation_uuid in calculation_graph.calculations_to_run:

            calculation = calculation_graph.calculations_to_run[calculation_uuid]

            value_protocol = calculation.protocols[calculation.final_value_reference.output_protocol_id]
            uncertainty_protocol = calculation.protocols[calculation.final_uncertainty_reference.output_protocol_id]

            value = value_protocol.get_output_value(calculation.final_value_reference)
            uncertainty = uncertainty_protocol.get_output_value(calculation.final_uncertainty_reference)

            property_type = calculation.property_type

            calculated_property = CalculatedPhysicalProperty()

            calculated_property.substance = calculation.physical_property.substance
            calculated_property.thermodynamic_state = calculation.physical_property.thermodynamic_state

            calculated_property.type = property_type
            calculated_property.phase = PropertyPhase.Liquid

            calculated_property.fidelity = CalculationFidelity.DirectSimulation

            calculated_property.value = value
            calculated_property.uncertainty = uncertainty

            self.calculated_properties[str(calculated_property.substance)].append(calculated_property)
