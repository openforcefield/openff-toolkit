# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Density Definition API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import logging
import mdtraj

from simtk import unit

from openforcefield.utils import doc_inherit

from openforcefield.properties.properties import PhysicalProperty

from openforcefield.properties.datasets import register_thermoml_property

from openforcefield.properties.estimator import CalculationSchema, register_estimable_property
from openforcefield.properties.estimator.components import protocols, groups
from openforcefield.properties.estimator.components.protocols import AverageTrajectoryProperty, \
    ProtocolInputReference, register_calculation_protocol


# =============================================================================================
# Custom Protocol Building Blocks
# =============================================================================================

@register_calculation_protocol()
class ExtractAverageDensity(AverageTrajectoryProperty):
    """Extracts the average density from a simulation trajectory.

    Todo
    ----
    Refactor this into a general 'ExtractAverageStatistic' class which
        can live in the protocols namespace.

    """

    def __init__(self):

        super().__init__()

        self._system = None

    @protocols.BaseProtocol.InputPipe
    def system(self, value):
        pass

    def execute(self, directory):

        logging.info('Extracting densities: ' + directory)

        if super(ExtractAverageDensity, self).execute(directory) is None:
            return False

        mass_list = []

        for atom_index in range(self._system.getNumParticles()):

            mass = self._system.getParticleMass(atom_index)
            mass /= (unit.gram / unit.mole)

            mass_list.append(mass)

        densities = mdtraj.density(self.trajectory, mass_list)

        self._value, self._uncertainty = self.calculate_average_and_error(densities)

        self._value *= unit.kilogram * unit.meter ** -3
        self._uncertainty *= unit.kilogram * unit.meter ** -3

        logging.info('Extracted densities: ' + directory)

        return True


# =============================================================================================
# Density
# =============================================================================================

@register_estimable_property()
@register_thermoml_property(thermoml_string='Mass density, kg/m3')
class Density(PhysicalProperty):
    """A class representation of a density property"""

    @staticmethod
    def get_calculation_schema():
        schema = CalculationSchema(type(Density))

        # Initial coordinate and topology setup.
        build_coordinates = protocols.BuildCoordinatesPackmol()
        build_coordinates.id = 'build_coordinates'

        build_coordinates.input_references = [
            ProtocolInputReference('substance', 'global', 'substance')
        ]

        schema.protocols[build_coordinates.id] = build_coordinates

        assign_topology = protocols.BuildSmirnoffTopology()
        assign_topology.id = 'build_topology'

        assign_topology.input_references = [
            # Globals
            ProtocolInputReference('force_field', 'global', 'force_field'),
            # Locals
            ProtocolInputReference('topology', build_coordinates.id, 'topology'),
            ProtocolInputReference('molecules', build_coordinates.id, 'molecules')
        ]

        schema.protocols[assign_topology.id] = assign_topology

        energy_minimisation = protocols.RunEnergyMinimisation()
        energy_minimisation.id = 'energy_minimisation'

        # Equilibration
        energy_minimisation.input_references = [
            # Locals
            ProtocolInputReference('positions', build_coordinates.id, 'positions'),
            ProtocolInputReference('topology', build_coordinates.id, 'topology'),
            ProtocolInputReference('system', assign_topology.id, 'system')
        ]

        schema.protocols[energy_minimisation.id] = energy_minimisation

        npt_equilibration = protocols.RunOpenMMSimulation()
        npt_equilibration.id = 'npt_equilibration'

        npt_equilibration.ensemble = protocols.RunOpenMMSimulation.Ensemble.NPT

        # Debug settings.
        npt_equilibration.steps = 2
        npt_equilibration.output_frequency = 1

        npt_equilibration.input_references = [
            # Globals
            ProtocolInputReference('thermodynamic_state', 'global', 'thermodynamic_state'),
            # Locals
            ProtocolInputReference('positions', energy_minimisation.id, 'final_positions'),
            ProtocolInputReference('topology', build_coordinates.id, 'topology'),
            ProtocolInputReference('system', assign_topology.id, 'system')
        ]

        schema.protocols[npt_equilibration.id] = npt_equilibration

        # Production

        npt_production = protocols.RunOpenMMSimulation()
        npt_production.id = 'npt_production'

        npt_production.ensemble = protocols.RunOpenMMSimulation.Ensemble.NPT

        # Debug settings.
        npt_production.steps = 200
        npt_production.output_frequency = 20

        npt_production.input_references = [
            # Globals
            ProtocolInputReference('thermodynamic_state', 'global', 'thermodynamic_state'),
            # Locals
            ProtocolInputReference('positions', npt_equilibration.id, 'final_positions'),
            ProtocolInputReference('topology', build_coordinates.id, 'topology'),
            ProtocolInputReference('system', assign_topology.id, 'system')
        ]

        schema.protocols[npt_production.id] = npt_production

        # Analysis

        extract_density = ExtractAverageDensity()
        extract_density.id = 'extract_density'

        extract_density.input_references = [
            # Globals
            ProtocolInputReference('thermodynamic_state', 'global', 'thermodynamic_state'),
            # Locals
            ProtocolInputReference('positions', npt_production.id, 'final_positions'),
            ProtocolInputReference('trajectory_path', npt_production.id, 'trajectory'),
            ProtocolInputReference('topology', build_coordinates.id, 'topology'),
            ProtocolInputReference('system', assign_topology.id, 'system')
        ]

        schema.protocols[npt_production.id] = npt_production

        # Set up a conditional group to ensure convergence of uncertainty

        converge_uncertainty = groups.ConditionalGroup({
            npt_production.id: npt_production,
            extract_density.id: extract_density
        })

        converge_uncertainty.id = 'converge_uncertainty'

        # TODO: Replace with a general global:convergence_criteria
        condition = groups.ConditionalGroup.Condition()

        condition.left_hand_reference = ProtocolInputReference('', extract_density.id, 'uncertainty')
        condition.right_hand_reference = ProtocolInputReference('', 'global', 'uncertainty')

        condition.condition_type = groups.ConditionalGroup.ConditionType.LessThan

        converge_uncertainty.conditions.append(condition)

        schema.groups[converge_uncertainty.id] = converge_uncertainty

        # Define where the final values come from.
        schema.final_value_reference = ProtocolInputReference('', extract_density.id, 'value')
        schema.final_uncertainty_reference = ProtocolInputReference('', extract_density.id, 'uncertainty')

        schema.build()

        return schema
