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

from openforcefield.properties.properties import PhysicalProperty

from openforcefield.properties.datasets import register_thermoml_property

from openforcefield.properties.estimator import CalculationSchema, register_estimable_property
from openforcefield.properties.estimator.components import protocols, groups
from openforcefield.properties.estimator.components.protocols import AverageTrajectoryProperty, \
    ProtocolInputReference, register_calculation_protocol, PropertyCalculatorException


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

        base_exception = super(ExtractAverageDensity, self).execute(directory)

        if isinstance(base_exception, PropertyCalculatorException):
            return base_exception

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
        build_coordinates = protocols.BuildCoordinatesPackmol()
        build_coordinates.id = 'build_coordinates'

        build_coordinates.input_references = [
            # Globals
            ProtocolInputReference(input_property_name='substance',
                                   output_protocol_id='global',
                                   output_property_name='substance')
        ]

        schema.protocols[build_coordinates.id] = build_coordinates.schema

        assign_topology = protocols.BuildSmirnoffTopology()
        assign_topology.id = 'build_topology'

        assign_topology.input_references = [
            # Globals
            ProtocolInputReference(input_property_name='force_field_path',
                                   output_protocol_id='global',
                                   output_property_name='force_field_path'),
            # Locals
            ProtocolInputReference(input_property_name='coordinate_file',
                                   output_protocol_id=build_coordinates.id,
                                   output_property_name='coordinate_file'),

            ProtocolInputReference(input_property_name='molecules',
                                   output_protocol_id=build_coordinates.id,
                                   output_property_name='molecules')
        ]

        schema.protocols[assign_topology.id] = assign_topology.schema

        energy_minimisation = protocols.RunEnergyMinimisation()
        energy_minimisation.id = 'energy_minimisation'

        # Equilibration
        energy_minimisation.input_references = [
            # Locals
            ProtocolInputReference(input_property_name='input_coordinate_file',
                                   output_protocol_id=build_coordinates.id,
                                   output_property_name='coordinate_file'),

            ProtocolInputReference(input_property_name='system',
                                   output_protocol_id=assign_topology.id,
                                   output_property_name='system')
        ]

        schema.protocols[energy_minimisation.id] = energy_minimisation.schema

        npt_equilibration = protocols.RunOpenMMSimulation()
        npt_equilibration.id = 'npt_equilibration'

        npt_equilibration.ensemble = protocols.RunOpenMMSimulation.Ensemble.NPT

        # Debug settings.
        npt_equilibration.steps = 2
        npt_equilibration.output_frequency = 1

        npt_equilibration.input_references = [
            # Globals
            ProtocolInputReference(input_property_name='thermodynamic_state',
                                   output_protocol_id='global',
                                   output_property_name='thermodynamic_state'),
            # Locals
            ProtocolInputReference(input_property_name='input_coordinate_file',
                                   output_protocol_id=energy_minimisation.id,
                                   output_property_name='output_coordinate_file'),

            ProtocolInputReference(input_property_name='system',
                                   output_protocol_id=assign_topology.id,
                                   output_property_name='system')
        ]

        schema.protocols[npt_equilibration.id] = npt_equilibration.schema

        # Production

        npt_production = protocols.RunOpenMMSimulation()
        npt_production.id = 'npt_production'

        npt_production.ensemble = protocols.RunOpenMMSimulation.Ensemble.NPT

        # Debug settings.
        npt_production.steps = 200
        npt_production.output_frequency = 20

        npt_production.input_references = [
            # Globals
            ProtocolInputReference(input_property_name='thermodynamic_state',
                                   output_protocol_id='global',
                                   output_property_name='thermodynamic_state'),
            # Locals
            ProtocolInputReference(input_property_name='input_coordinate_file',
                                   output_protocol_id=npt_equilibration.id,
                                   output_property_name='output_coordinate_file'),

            ProtocolInputReference(input_property_name='system',
                                   output_protocol_id=assign_topology.id,
                                   output_property_name='system')
        ]

        schema.protocols[npt_production.id] = npt_production.schema

        # Analysis

        extract_density = ExtractAverageDensity()
        extract_density.id = 'extract_density'

        extract_density.input_references = [
            # Globals
            ProtocolInputReference(input_property_name='thermodynamic_state',
                                   output_protocol_id='global',
                                   output_property_name='thermodynamic_state'),
            # Locals
            ProtocolInputReference(input_property_name='input_coordinate_file',
                                   output_protocol_id=npt_production.id,
                                   output_property_name='output_coordinate_file'),

            ProtocolInputReference(input_property_name='trajectory_path',
                                   output_protocol_id=npt_production.id,
                                   output_property_name='trajectory'),

            ProtocolInputReference(input_property_name='system',
                                   output_protocol_id=assign_topology.id,
                                   output_property_name='system')
        ]

        schema.protocols[extract_density.id] = extract_density.schema

        # Set up a conditional group to ensure convergence of uncertainty
        converge_uncertainty = groups.ConditionalGroup([
            npt_production.id,
            extract_density.id
        ])
        converge_uncertainty.id = 'converge_uncertainty'

        converge_uncertainty.input_references = [
            # Locals
            ProtocolInputReference(input_property_name='left_hand_value',
                                   output_protocol_id='extract_density',
                                   output_property_name='uncertainty'),
            # Globals
            ProtocolInputReference(input_property_name='right_hand_value',
                                   output_protocol_id='global',
                                   output_property_name='uncertainty'),
        ]

        schema.groups[converge_uncertainty.id] = converge_uncertainty.schema

        # Define where the final values come from.
        schema.final_value_reference = ProtocolInputReference(input_property_name=None,
                                                              output_protocol_id=extract_density.id,
                                                              output_property_name='value')

        schema.final_uncertainty_reference = ProtocolInputReference(input_property_name=None,
                                                                    output_protocol_id=extract_density.id,
                                                                    output_property_name='uncertainty')

        return schema
