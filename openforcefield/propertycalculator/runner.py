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

import os
import copy
import logging

from queue import Queue
from threading import Thread
from os import path
from xml.etree import ElementTree

from openforcefield.propertycalculator import protocols
from openforcefield.properties import PropertyType


# =============================================================================================
# Calculation Classes
# =============================================================================================

# TODO: Move into better home.
class XmlNodeMissingException(Exception):

    def __init__(self, node_name):

        message = 'The calculation template does not contain a <' + str(node_name) + '> node.'
        super().__init__(message)


class DirectPropertyCalculationTemplate:
    """
    Defines the set of protocols required to calculate a certain property.
    """
    def __init__(self):

        self.property_type = PropertyType.Undefined
        self.label = ''

        self.protocols = []

    @classmethod
    def from_xml(cls, xml_node, existing_templates):

        return_value = cls()

        # Gather up all possible identifiers
        property_node = xml_node.find('property')

        if property_node is None:
            raise XmlNodeMissingException('property')

        return_value.property_type = PropertyType[property_node.text]

        label_node = xml_node.find('label')

        if label_node is None:
            raise XmlNodeMissingException('label')

        return_value.label = label_node.text

        inherits_node = xml_node.find('inherits')

        if inherits_node is not None:

            if inherits_node.text not in existing_templates:

                raise RuntimeError('The ' + return_value.label + ' template tries to inherit from '
                                                                 'a non-existent template.')

            return_value.protocols.extend(existing_templates[inherits_node.text].protocols)

        protocol_nodes = xml_node.find('protocols')

        if label_node is None:
            raise XmlNodeMissingException('protocols')

        for protocol_node in protocol_nodes:

            protocol = None

            # TODO: There must be a cleaner way to do this...
            if protocol_node.tag == protocols.BuildLiquidCoordinates.__name__:
                protocol = protocols.BuildLiquidCoordinates.from_xml(protocol_node)
            elif protocol_node.tag == protocols.BuildSmirnoffTopology.__name__:
                protocol = protocols.BuildSmirnoffTopology.from_xml(protocol_node)
            elif protocol_node.tag == protocols.RunEnergyMinimisation.__name__:
                protocol = protocols.RunEnergyMinimisation.from_xml(protocol_node)
            elif protocol_node.tag == protocols.RunNVTSimulation.__name__:
                protocol = protocols.RunNVTSimulation.from_xml(protocol_node)
            elif protocol_node.tag == protocols.RunNPTSimulation.__name__:
                protocol = protocols.RunNPTSimulation.from_xml(protocol_node)
            elif protocol_node.tag == protocols.ExtractStatistics.__name__:
                protocol = protocols.ExtractStatistics.from_xml(protocol_node)
            else:

                raise RuntimeError('The ' + return_value.label + ' contains an undefined '
                                                                 'protocol type: ' + protocol_node.tag)

            return_value.protocols.append(protocol)

        return return_value


class DirectCalculationGraph:
    """
    A hierarchical structure for storing all protocols that
    need to be executed.
    """
    class TreeNode:

        def __init__(self, index, depth):

            self.index = index
            self.depth = depth

            self.children = []

    def __init__(self):

        self._master_protocols_list = []
        self._nodes = []

        self._leaf_nodes = []

    @property
    def nodes(self):
        return self._nodes

    @property
    def leaf_protocols(self):
        return [self._master_protocols_list[node.index] for node in self._leaf_nodes]

    def _insert_root_node(self, protocol):

        for node in self._nodes:

            existing_protocol = self._master_protocols_list[node.index]

            if existing_protocol.compare_to(protocol) is False:
                continue

            return node

        self._master_protocols_list.append(protocol)

        new_node = self.TreeNode(len(self._master_protocols_list) - 1, 0)
        self._nodes.append(new_node)

        return new_node

    def _insert_node(self, protocol, tree_node):

        for node in tree_node.children:

            existing_protocol = self._master_protocols_list[node.index]

            if existing_protocol.compare_to(protocol) is False:
                continue

            return node

        self._master_protocols_list.append(protocol)

        new_node = self.TreeNode(len(self._master_protocols_list) - 1, tree_node.depth + 1)
        tree_node.children.append(new_node)

        return new_node

    def add_calculation(self, calculation):

        tree_node = self._insert_root_node(calculation.protocols[0])

        for protocol in calculation.protocols[1:]:
            tree_node = self._insert_node(protocol, tree_node)

        if tree_node is not None and tree_node not in self._leaf_nodes:
            self._leaf_nodes.append(tree_node)

    @staticmethod
    def worker_thread(task_queue):

        while True:

            parent_tree, tree_node, input_data = task_queue.get()
            parent_tree.execute(task_queue, input_data, tree_node)

            task_queue.task_done()

    def execute(self, task_queue, protocol_data, tree_node=None):

        nodes = self._nodes if tree_node is None else tree_node.children

        if not path.isdir(protocol_data.root_directory):
            os.makedirs(protocol_data.root_directory)

        if tree_node is not None:

            protocol = self._master_protocols_list[tree_node.index]
            protocol_data = protocol.execute(protocol_data)

        for index, node in enumerate(nodes):

            input_data = protocol_data

            if len(nodes) > 1:
                input_data = protocols.ProtocolData.clone(protocol_data)

            input_data.root_directory = path.join(input_data.root_directory, str(index))

            task_queue.put((self, node, input_data))


# =============================================================================================
# Property Estimator
# =============================================================================================

class PropertyCalculationRunner:
    """
    The class responsible for deciding at which fidelity a property
    is calculated, as well as launching the calculations themselves.
    """

    def __init__(self, number_of_threads=1):

        self._pending_calculations = {}
        self._finished_calculations = {}

        self._number_of_threads = number_of_threads

        self.calculator_templates = {}

        self.load_calculator_templates()

    def load_calculator_templates(self):
        """
        Loads the property calculation templates from the default xml file.
        """
        templates = {}

        with open('../data/properties/templates.xml') as templates_file:

            xml = templates_file.read()
            root_node = ElementTree.fromstring(xml)

            for template_node in root_node:

                calculator_template = DirectPropertyCalculationTemplate.from_xml(template_node,
                                                                                 templates)

                templates[calculator_template.label] = calculator_template

        for template_label in templates:

            template = templates[template_label]

            if template.property_type is PropertyType.Undefined:
                continue

            self.calculator_templates[template.property_type] = template

    def run(self, measured_properties, parameter_set):

        # Make a copy of the properties array so we
        # can easily safely changes to it.
        properties_to_measure = []
        properties_to_measure.extend(measured_properties)

        # Reset the calculation arrays
        self._pending_calculations = {}
        self._finished_calculations = {}

        # Set up the array to track any pending or
        # finished calculations
        for measured_property in properties_to_measure:

            substance_tag = measured_property.substance.to_tag()

            if substance_tag not in self._pending_calculations:
                self._pending_calculations[substance_tag] = []
            if substance_tag not in self._finished_calculations:
                self._finished_calculations[substance_tag] = []

        # First see if any of these properties could be calculated
        # from a surrogate model. If they can, remove them the list
        # of remaining properties to be calculated.
        self.attempt_surrogate_extrapolation(parameter_set, properties_to_measure)

        # Exit early if everything is covered by the surrogate
        if len(properties_to_measure) == 0:
            return self._finished_calculations

        # Do a similar thing to check whether reweighting could
        # be employed here for any of the properties.
        self.attempt_simulation_reweighting(parameter_set, properties_to_measure)

        # Exit early if everything is covered by reweighting
        if len(properties_to_measure) == 0:
            return self._finished_calculations

        # We are now left with a list of properties that are
        # going to have to be calculated by direct simulation.
        #
        # First split the desired properties into a dict where
        # the mixture tag is the key. This structure would be the
        # starting point in grouping together similar calculations
        # that need to be performed.
        properties_by_mixture = {}

        for measured_property in properties_to_measure:

            substance_tag = measured_property.substance.to_tag()

            if substance_tag not in properties_by_mixture:
                properties_by_mixture[substance_tag] = []

            if measured_property.type in properties_by_mixture[substance_tag]:
                continue

            properties_by_mixture[substance_tag].append(measured_property)

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
        calculation_graphs = self.build_calculation_graph(properties_by_mixture)

        # Once a list of unique list of protocols to be executed has been constructed
        # the protocols can be fired.

        self.execute_calculation_graph(calculation_graphs,
                                       parameter_set)

        for mixture in calculation_graphs:

            calculation_graph = calculation_graphs[mixture]

            for leaf_protocol in calculation_graph.leaf_protocols:

                if type(leaf_protocol) is not protocols.ExtractStatistics:
                    continue

                for calculated_property in leaf_protocol.calculated_properties:
                    self._finished_calculations[mixture].append(leaf_protocol.calculated_properties[calculated_property])

        # TODO : Extract calculated properties from ProtocolData's.

        # The results from these calculations would then be stored in some way so that
        # the data required for reweighting is retained.
        #
        # At this point, two things could happen:
        #
        #   1) The simulation results could be fed into the surrogate model
        #   2) The results would be passed back to the 'client'

        return self._finished_calculations

    def build_calculation_graph(self, properties_by_mixture):

        calculation_graphs = {}

        for mixture in properties_by_mixture:

            calculation_graph = DirectCalculationGraph()

            for property_to_calculate in properties_by_mixture[mixture]:

                if property_to_calculate.type not in self.calculator_templates:

                    logging.warning('The property calculator does not support ' +
                                    property_to_calculate.type + ' calculations.')

                    continue

                calculation = copy.deepcopy(self.calculator_templates[property_to_calculate.type])

                for protocol in calculation.protocols:
                    protocol.set_measured_property(property_to_calculate)

                calculation_graph.add_calculation(calculation)

            calculation_graphs[mixture] = calculation_graph

        return calculation_graphs

    def execute_calculation_graph(self, calculation_graphs, force_field):

        task_queue = Queue()

        logging.info('Spawning ' + str(self._number_of_threads) + ' worker threads.')

        # Set up the pool of threads.
        for i in range(self._number_of_threads):

            worker = Thread(target=DirectCalculationGraph.worker_thread, args=(task_queue,))
            worker.setDaemon(True)
            worker.start()

        for mixture in calculation_graphs:

            calculation_graph = calculation_graphs[mixture]

            if len(calculation_graph.nodes) == 0:
                continue

            protocol_data = protocols.ProtocolData()

            protocol_data.root_directory = path.join('property-data', mixture)

            protocol_data.substance_tag = mixture
            protocol_data.force_field = force_field

            calculation_graph.execute(task_queue, protocol_data)

        task_queue.join()

    def attempt_surrogate_extrapolation(self, parameter_set, properties_to_measure):
        """Attempts to calculate properties from a surrogate model"""

        # This loop would be distributed over a number of threads /
        # nodes depending on the available architecture.
        for measured_property in properties_to_measure[:]:

            surrogate_result = self.perform_surrogate_extrapolation(measured_property, parameter_set)

            if surrogate_result is None:
                continue

            substance_tag = measured_property.substance.to_tag()

            self._finished_calculations[substance_tag].append(surrogate_result)
            properties_to_measure.remove(measured_property)

    def perform_surrogate_extrapolation(self, measured_property, parameter_set):
        """
        A placeholder method that would be used to spawn the surrogate
        model backend.
        """
        return None

    def attempt_simulation_reweighting(self, parameter_set, properties_to_measure):
        """Attempts to calculate properties by reweighting existing data"""

        for measured_property in properties_to_measure[:]:

            reweighting_result = self.perform_reweighting(measured_property, parameter_set)

            if reweighting_result is None:
                continue

            substance_tag = measured_property.substance.to_tag()

            self._finished_calculations[substance_tag].append(reweighting_result)
            properties_to_measure.remove(measured_property)

    def perform_reweighting(self, measured_property, parameter_set):
        """
        A placeholder method that would be used to spawn the property
        reweighting backend.
        """
        return None
