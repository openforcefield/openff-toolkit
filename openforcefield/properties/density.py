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

from openforcefield.properties.datasets import register_thermoml_property
from openforcefield.properties.estimator import CalculationSchema, register_estimable_property
from openforcefield.properties.estimator.components import protocols, groups
from openforcefield.properties.estimator.components.protocols import AverageTrajectoryProperty, \
    register_calculation_protocol, ProtocolPath, PropertyCalculatorException
from openforcefield.properties.properties import PhysicalProperty
from openforcefield.utils import statistics


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

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        self._system = None

    @protocols.BaseProtocol.InputPipe
    def system(self, value):
        pass

    def execute(self, directory):

        logging.info('Extracting densities: ' + directory)

        base_exception = super(ExtractAverageDensity, self).execute(directory)

        if isinstance(base_exception, PropertyCalculatorException):
            return base_exception

        mass_list = []

        for atom_index in range(self._system.getNumParticles()):

            mass = self._system.getParticleMass(atom_index)
            mass /= (unit.gram / unit.mole)

            mass_list.append(mass)

        densities = mdtraj.density(self.trajectory, mass_list)

        densities, self._equilibration_index, self._statistical_inefficiency = \
            statistics.decorrelate_time_series(densities)

        self._value, self._uncertainty = self._perform_bootstrapping(densities)

        self._value *= unit.kilogram * unit.meter ** -3
        self._uncertainty *= unit.kilogram * unit.meter ** -3

        logging.info('Extracted densities: ' + directory)

        return self._get_output_dictionary()


# =============================================================================================
# Density
# =============================================================================================

@register_estimable_property()
@register_thermoml_property(thermoml_string='Mass density, kg/m3')
class Density(PhysicalProperty):
    """A class representation of a density property"""

    @staticmethod
    def get_default_calculation_schema():

        schema = CalculationSchema(property_type=Density.__name__)
        schema.id = '{}{}'.format(Density.__name__, 'Schema')

        # Initial coordinate and topology setup.
        build_coordinates = protocols.BuildCoordinatesPackmol('build_coordinates')

        build_coordinates.substance = ProtocolPath('substance', 'global')

        schema.protocols[build_coordinates.id] = build_coordinates.schema

        assign_topology = protocols.BuildSmirnoffTopology('build_topology')

        assign_topology.force_field_path = ProtocolPath('force_field_path', 'global')

        assign_topology.coordinate_file = ProtocolPath('coordinate_file', build_coordinates.id)
        assign_topology.molecules = ProtocolPath('molecules', build_coordinates.id)

        schema.protocols[assign_topology.id] = assign_topology.schema

        # Equilibration
        energy_minimisation = protocols.RunEnergyMinimisation('energy_minimisation')

        energy_minimisation.input_coordinate_file = ProtocolPath('coordinate_file', build_coordinates.id)
        energy_minimisation.system = ProtocolPath('system', assign_topology.id)

        schema.protocols[energy_minimisation.id] = energy_minimisation.schema

        npt_equilibration = protocols.RunOpenMMSimulation('npt_equilibration')

        npt_equilibration.ensemble = protocols.RunOpenMMSimulation.Ensemble.NPT

        npt_equilibration.steps = 100000  # Debug settings.
        npt_equilibration.output_frequency = 5000  # Debug settings.

        npt_equilibration.thermodynamic_state = ProtocolPath('thermodynamic_state', 'global')

        npt_equilibration.input_coordinate_file = ProtocolPath('output_coordinate_file', energy_minimisation.id)
        npt_equilibration.system = ProtocolPath('system', assign_topology.id)

        schema.protocols[npt_equilibration.id] = npt_equilibration.schema

        # Production
        npt_production = protocols.RunOpenMMSimulation('npt_production')

        npt_production.ensemble = protocols.RunOpenMMSimulation.Ensemble.NPT

        npt_production.steps = 1000000  # Debug settings.
        npt_production.output_frequency = 50000  # Debug settings.

        npt_production.thermodynamic_state = ProtocolPath('thermodynamic_state', 'global')

        npt_production.input_coordinate_file = ProtocolPath('output_coordinate_file', npt_equilibration.id)
        npt_production.system = ProtocolPath('system', assign_topology.id)

        # Analysis
        extract_density = ExtractAverageDensity('extract_density')

        extract_density.thermodynamic_state = ProtocolPath('thermodynamic_state', 'global')

        extract_density.input_coordinate_file = ProtocolPath('output_coordinate_file', npt_production.id)
        extract_density.trajectory_path = ProtocolPath('trajectory', npt_production.id)
        extract_density.system = ProtocolPath('system', assign_topology.id)

        # Set up a conditional group to ensure convergence of uncertainty
        converge_uncertainty = groups.ConditionalGroup('converge_uncertainty')
        converge_uncertainty.add_protocols(npt_production, extract_density)

        converge_uncertainty.left_hand_value = ProtocolPath('uncertainty',
                                                            converge_uncertainty.id,
                                                            extract_density.id)

        converge_uncertainty.right_hand_value = ProtocolPath('target_uncertainty', 'global')

        converge_uncertainty.condition_type = groups.ConditionalGroup.ConditionType.LessThan

        converge_uncertainty.max_iterations = 1

        schema.protocols[converge_uncertainty.id] = converge_uncertainty.schema

        # Finally, extract uncorrelated data
        extract_uncorrelated_trajectory = protocols.ExtractUncorrelatedTrajectoryData('extract_traj')

        extract_uncorrelated_trajectory.statistical_inefficiency = ProtocolPath('statistical_inefficiency',
                                                                                converge_uncertainty.id,
                                                                                extract_density.id)

        extract_uncorrelated_trajectory.equilibration_index = ProtocolPath('equilibration_index',
                                                                           converge_uncertainty.id,
                                                                           extract_density.id)

        extract_uncorrelated_trajectory.input_coordinate_file = ProtocolPath('output_coordinate_file',
                                                                             converge_uncertainty.id,
                                                                             npt_production.id)

        extract_uncorrelated_trajectory.input_trajectory_path = ProtocolPath('trajectory',
                                                                             converge_uncertainty.id,
                                                                             npt_production.id)

        schema.protocols[extract_uncorrelated_trajectory.id] = extract_uncorrelated_trajectory.schema

        # Define where the final values come from.
        schema.final_value_source = ProtocolPath('value', converge_uncertainty.id, extract_density.id)
        schema.final_uncertainty_source = ProtocolPath('uncertainty', converge_uncertainty.id, extract_density.id)

        schema.final_trajectory_source = ProtocolPath('output_trajectory_path', extract_uncorrelated_trajectory.id)

        return schema
