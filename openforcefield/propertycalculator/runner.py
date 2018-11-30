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

from enum import Enum
from queue import Queue
from threading import Thread
from os import path
from xml.etree import ElementTree

from openforcefield.utils.exceptions import XmlNodeMissingException
from openforcefield.propertycalculator import protocols
from openforcefield.properties import PropertyType, PropertyPhase, CalculationFidelity, \
                                      CalculatedPhysicalProperty

from openforcefield.datasets.property_dataset import CalculatedPropertySet


# =============================================================================================
# Calculation Classes
# =============================================================================================

class DirectCalculationTemplate:
    """Defines the set of protocols required to calculate a certain property.
    """
    def __init__(self):

        self.property_type = PropertyType.Undefined
        self.id = None

        self.protocols = []
        self.protocol_ids = set()

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

            raise Exception('A calculation with id ' + return_value.id +
                            ' has been defined more than once.')

        protocol_nodes = xml_node.find('protocols')

        if protocol_nodes is None:
            raise XmlNodeMissingException('protocols')

        # Load in all of the protocols.
        for protocol_node in protocol_nodes:

            # Handle the case where we import the protocols
            # from an already defined calculation.
            if protocol_node.tag == 'import':

                if protocol_node.text not in existing_templates:

                    raise RuntimeError('The ' + return_value.id + ' template tries to inherit from '
                                                                  'a non-existent template.')

                return_value.protocols.extend(copy.deepcopy(existing_templates[protocol_node.text].protocols))

                for protocol in return_value.protocols:
                    return_value.protocol_ids.add(protocol.id)

                continue

            # Try to automatically find which class corresponds to the xml tag.
            protocol_class = getattr(protocols, protocol_node.tag)

            if protocol_class is None:

                raise RuntimeError('The ' + return_value.id + ' calculation contains an undefined'
                                                              ' protocol of type: ' + protocol_node.tag)

            # Make sure the found class is actually a protocol.
            if isinstance(protocol_class, protocols.Protocol):
                raise RuntimeError(protocol_node.tag + ' is not a child of Protocol')

            protocol = protocol_class.from_xml(protocol_node)
            # protocol.id = return_value.id + '|' + protocol.id

            if protocol.id in return_value.protocol_ids:

                raise Exception('The ' + return_value.id + ' calculation defines two protocols'
                                                           ' with the same id: ' + protocol.id)

            return_value.protocol_ids.add(protocol.id)
            return_value.protocols.append(protocol)

        cls.validate_inputs(return_value)
        return return_value

    @staticmethod
    def validate_inputs(template):

        for protocol in template.protocols:

            for input_reference in protocol.input_references:

                if input_reference.protocol_id == 'global':
                    # We handle global inputs separately
                    continue

                # Make sure the other protocol whose output we are interested
                # in actually exists.
                if input_reference.protocol_id not in template.protocol_ids:

                    raise Exception('The ' + protocol.id + ' protocol of the ' + template.id + ' template tries to ' +
                                    'take input from a non-existent protocol: ' + input_reference.protocol_id)

                other_protocol = next(x for x in template.protocols if x.id == input_reference.protocol_id)

                # Make sure the other protocol definitely has the requested output.
                if input_reference.property_name not in other_protocol.provided_outputs:

                    raise Exception('The ' + other_protocol.id + ' protocol does not provide an ' +
                                    input_reference.property_name + ' output.')


class DirectCalculation:
    """Defines the property to calculate and template to follow.

    Parameters
    ----------
    physical_property: Protocol
        The protocol this node will execute.
    force_field: ForceField
        The force field to use for this calculation.
    template: CalculationNode, optional
        The parent of this node.
    """
    def __init__(self, physical_property, force_field, template):

        self.physical_property = physical_property
        self.template = copy.deepcopy(template)

        self.leaf_node_ids = []

        self.uuid = str(uuid.uuid4())

        thermodynamic_state = physical_property.thermodynamic_state
        substance = physical_property.substance

        for protocol in self.protocols:

            # Try to set global properties on each of the protocols
            for input_reference in protocol.input_references:

                if input_reference.protocol_id != 'global':
                    continue

                if input_reference.property_name == 'thermodynamic_state':
                    protocol.thermodynamic_state = thermodynamic_state
                elif input_reference.property_name == 'substance':
                    protocol.substance = substance
                elif input_reference.property_name == 'force_field':
                    protocol.force_field = force_field
                else:
                    raise Exception('Invalid global property: ' + input_reference.property_name)

            protocol.set_uuid(self.uuid)

    @property
    def protocols(self):
        """list(Protocol): The list of protocols to follow to calculate a property."""
        return self.template.protocols


class DirectCalculationGraph:
    """A hierarchical structure for storing all protocols that need to be executed.
    """

    class CalculationNode:
        """A node in the calculation graph.

        Each node represents a protocol to be executed.

        Parameters
        ----------
        protocol: Protocol
            The protocol this node will execute.
        parent: str, optional
            The id of the parent of this node.
        """
        class Status(Enum):
            """Describes the state of a calculation node."""
            Queueing = 0
            Running = 1
            Failed = 2
            Finished = 3

        def __init__(self, protocol, parent_node):

            self.protocol = protocol

            self.directory = protocol.id if parent_node is None else \
                path.join(parent_node.directory, protocol.id)

            self.status = self.Status.Queueing

            self.children = []

        @property
        def id(self):
            """str: Returns the id of the protocol to execute."""
            return self.protocol.id

    def __init__(self):

        self._nodes_by_id = {}
        self._node_id_map = {}

        self._nodes = []

    def _insert_node(self, protocol, parent_node_id=None):
        """Insert a protocol node into the tree.

        Parameters
        ----------
        protocol : Protocol
            The protocol to insert.
        parent_node_id : str
            The id of node in which to insert the new protocol.

        Returns
        ----------
        str
            The id of the node that was inserted. If an identical node already exists,
            the id of that node is returned instead.
        """

        nodes = self._nodes if parent_node_id is None else \
            self._nodes_by_id[parent_node_id].children

        for input_reference in protocol.input_references:

            if input_reference.protocol_id == 'global':
                continue

            input_reference.protocol_id = self._node_id_map[input_reference.protocol_id]

        for node_id in nodes:

            node = self._nodes_by_id[node_id]

            if node.protocol.compare_to(protocol) is False:
                continue

            return node_id

        parent_node = None if parent_node_id is None else \
            self._nodes_by_id[parent_node_id]

        new_node = self.CalculationNode(protocol, parent_node)

        self._nodes_by_id[new_node.id] = new_node
        nodes.append(new_node.id)

        return new_node.id

    def add_calculation(self, calculation):
        """Insert a calculation into the calculation graph.

        Parameters
        ----------
        calculation : DirectCalculation
            The calculation to insert.
        """
        parent_node_id = None

        for protocol in calculation.protocols:

            parent_node_id = self._insert_node(protocol, parent_node_id)
            self._node_id_map[protocol.id] = parent_node_id

            if protocol.id != parent_node_id:
                protocol.set_uuid(parent_node_id.split('|')[0])

        if parent_node_id is not None:
            calculation.leaf_node_ids.append(parent_node_id)

    def execute(self, task_queue, node_id=None):
        """Execute the given node on a thread, then queue its children.

        Parameters
        ----------
        task_queue : Queue
            The queue of nodes to execute.
        node_id : str, optional
            The id of the node to execute. If none is passed,
            all root nodes are queued to be executed.
        """
        if node_id is not None:

            node_to_execute = self._nodes_by_id[node_id]

            if not path.isdir(node_to_execute.directory):
                os.makedirs(node_to_execute.directory)

            for input_reference in node_to_execute.protocol.input_references:

                if input_reference.protocol_id == 'global':
                    continue

                output_node = self._nodes_by_id[self._node_id_map[input_reference.protocol_id]]

                if output_node.status != DirectCalculationGraph.CalculationNode.Status.Finished:

                    raise Exception('A node is trying to take input '
                                    'from an unfinished node: ', output_node.id)

                input_value = getattr(output_node.protocol, input_reference.property_name)
                setattr(node_to_execute.protocol, input_reference.input_property, input_value)

            if node_to_execute.protocol.execute(node_to_execute.directory) is False:
                # This node failed to run, can't continue down the tree.
                return

            node_to_execute.status = DirectCalculationGraph.CalculationNode.Status.Finished

        child_ids = self._nodes if node_id is None else self._nodes_by_id[node_id].children

        for child_id in child_ids:
            task_queue.put((self, child_id))

    @staticmethod
    def worker_thread(task_queue):
        """A method to execute a node from the thread queue.

        Parameters
        ----------
        task_queue : Queue
            The queue of nodes to execute.
        """

        while True:

            parent_tree, node_id = task_queue.get()
            parent_tree.execute(task_queue, node_id)

            task_queue.task_done()


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

        if threads_per_simulation is None:
            threads_per_simulation = math.floor(maximum_threads / worker_threads)

        if worker_threads * threads_per_simulation > maximum_threads:

            raise ValueError('The total number of requested threads (' +
                             str(worker_threads * threads_per_simulation) +
                             ') must be less than the available ' + str(maximum_threads))

        logging.info(str(threads_per_simulation) + ' threads will be used per'
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
        # self.extract_calculated_properties(calculation_graphs)

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
        calculation_graph = DirectCalculationGraph()

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
        calculation_graphs : list(DirectCalculationGraph)
            A list of calculation graphs to execute.
        force_field : ForceField
            The force field parameters to use in any simulations.
        """
        task_queue = Queue()

        logging.info('Spawning ' + str(self._worker_threads) + ' worker threads.')

        # Set up the pool of threads.
        for i in range(self._worker_threads):

            worker = Thread(target=DirectCalculationGraph.worker_thread, args=(task_queue,))
            worker.setDaemon(True)
            worker.start()

        calculation_graph.execute(task_queue)

        task_queue.join()

    def extract_calculated_properties(self, calculation_graphs):
        """Extract any properties calculated by a calculation graph.

        Parameters
        ----------
        calculation_graphs : list(DirectCalculationGraph)
            A list of executed calculation graphs.
        """
        for mixture in calculation_graphs:

            calculation_graph = calculation_graphs[mixture]

            for leaf_protocol in calculation_graph.leaf_protocols:

                if not isinstance(leaf_protocol, protocols.AveragePropertyProtocol):
                    continue

                property_type = PropertyType.Undefined

                if isinstance(leaf_protocol, protocols.ExtractAverageDensity):
                    property_type = PropertyType.Density
                elif isinstance(leaf_protocol, protocols.ExtractAverageDielectric):
                    property_type = PropertyType.DielectricConstant

                if property_type is PropertyType.Undefined:
                    continue

                calculated_property = CalculatedPhysicalProperty()

                calculated_property.substance = leaf_protocol.substance
                calculated_property.thermodynamic_state = leaf_protocol.thermodynamic_state

                calculated_property.type = property_type
                calculated_property.phase = PropertyPhase.Liquid

                calculated_property.fidelity = CalculationFidelity.DirectSimulation

                calculated_property.value = leaf_protocol.value
                calculated_property.uncertainty = leaf_protocol.uncertainty

                self.calculated_properties[mixture].append(calculated_property)
