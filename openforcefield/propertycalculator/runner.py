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

from openforcefield.propertycalculator import protocols
from openforcefield.propertycalculator.protocols import ExtractableStatistics

from openforcefield.properties import PropertyType


# =============================================================================================
# Calculation Classes
# =============================================================================================


class DirectPropertyCalculation:
    """
    Defines the set of protocols required to calculate a certain property.
    """
    def __init__(self):

        self._property_type = PropertyType.Undefined
        self.protocols = []

    @property
    def property_type(self):
        return self._property_type


class DirectDensityCalculation(DirectPropertyCalculation):

    def __init__(self):

        super().__init__()

        self._property_type = PropertyType.Density

        build_liquid = protocols.BuildLiquidCoordinates()
        topology_protocol = protocols.BuildSmirnoffTopology()

        self.protocols.append((build_liquid, topology_protocol))

        energy_minimisation = protocols.RunEnergyMinimisation()
        run_nvt = protocols.RunNVTSimulation()

        self.protocols.append((energy_minimisation, run_nvt))

        run_npt = protocols.RunNPTSimulation()

        self.protocols.append(run_npt)

        statistics = protocols.ExtractableStatistics(ExtractableStatistics.Density)

        self.protocols.append(statistics)


class DirectDielectricCalculation(DirectPropertyCalculation):

    def __init__(self):
        super().__init__()

        self._property_type = PropertyType.DielectricConstant

        build_liquid = protocols.BuildLiquidCoordinates()
        topology_protocol = protocols.BuildSmirnoffTopology()

        self.protocols.append((build_liquid, topology_protocol))

        energy_minimisation = protocols.RunEnergyMinimisation()
        run_nvt = protocols.RunNVTSimulation()

        self.protocols.append((energy_minimisation, run_nvt))

        run_npt = protocols.RunNPTSimulation()

        self.protocols.append(run_npt)


class DirectCalculationGraph:
    """
    A hierarchical structure for storing all protocols that
    need to be executed.
    """
    class TreeNode:

        def __init__(self, index):

            self.index = index
            self.children = []

    def __init__(self):

        self._master_protocols_list = []
        self._nodes = []

        self._leaf_nodes = []

    def _insert_root_node(self, protocol):

            for node in self._nodes:

                existing_protocol = self._master_protocols_list[node.index]

                if type(existing_protocol) != type(protocol):
                    continue

                return node

            self._master_protocols_list.append(protocol)

            new_node = self.TreeNode(len(self._master_protocols_list) - 1)
            self._nodes.append(new_node)

    def _insert_node(self, protocol, tree_node):

        for node in tree_node.children:

            existing_protocol = self._master_protocols_list[node.index]

            if type(existing_protocol) != type(protocol):
                continue

            return node

        self._master_protocols_list.append(protocol)

        new_node = self.TreeNode(len(self._master_protocols_list) - 1)
        self._nodes.append(new_node)

    def add_calculation(self, calculation):

        tree_node = self._insert_root_node(calculation.protocols[0])

        for protocol in calculation.protocols[1:]:
            tree_node = self._insert_node(protocol, tree_node)

        if tree_node is not None and tree_node not in self._leaf_nodes:
            self._leaf_nodes.append(tree_node)


# =============================================================================================
# Property Estimator
# =============================================================================================

class PropertyCalculationRunner:
    """
    The class responsible for deciding at which fidelity a property
    is calculated, as well as launching the calculations themselves.
    """

    def __init__(self):

        self._pending_calculations = {}
        self._finished_calculations = {}

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

            properties_by_mixture[substance_tag].append(measured_property.type)

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

        calculation_graphs = {}

        for mixture in properties_by_mixture:

            calculation_graph = DirectCalculationGraph()

            for property_to_calculate in properties_by_mixture[mixture]:

                calculation = None

                if property_to_calculate is PropertyType.Density:
                    calculation = DirectDensityCalculation()
                elif property_to_calculate is PropertyType.DielectricConstant:
                    calculation = DirectDielectricCalculation()
                else:

                    logging.warning('The property calculator does not support ' +
                                    property_to_calculate + ' calculations.')

                    continue

                calculation_graph.add_calculation(calculation)

            calculation_graphs[mixture] = calculation_graph

        # Once a list of unique list of protocols to be executed has been constructed
        # the protocols can be fired.

        # TODO : Traverse graph and execute protocols.
        # TODO : Extract calculated properties from ProtocolData's.

        # The results from these calculations would then be stored in some way so that
        # the data required for reweighting is retained.
        #
        # At this point, two things could happen:
        #
        #   1) The simulation results could be fed into the surrogate model
        #   2) The results would be passed back to the 'client'

        return self._finished_calculations

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
