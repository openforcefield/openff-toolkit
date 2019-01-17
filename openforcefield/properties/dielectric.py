# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Dielectric Definition API.

Authors
-------
* Simon Boothroyd <simon.boothroyd@choderalab.org>

"""
# =============================================================================================
# GLOBAL IMPORTS
# =============================================================================================

import logging
import mdtraj

import arch.bootstrap

import numpy as np

from simtk import openmm, unit

from openforcefield.properties.properties import PhysicalProperty

from openforcefield.properties.datasets import register_thermoml_property

from openforcefield.properties.estimator import CalculationSchema, register_estimable_property
from openforcefield.properties.estimator.components import protocols, groups
from openforcefield.properties.estimator.components.protocols import AverageTrajectoryProperty, \
    register_calculation_protocol, ProtocolPath


# =============================================================================================
# Custom Protocol Building Blocks
# =============================================================================================

@register_calculation_protocol()
class ExtractAverageDielectric(AverageTrajectoryProperty):
    """Extracts the average dielectric constant from a simulation trajectory.
    """
    def __init__(self):
        super().__init__()

        self._system = None
        self._thermodynamic_state = None

    @protocols.BaseProtocol.InputPipe
    def system(self, value):
        pass

    @protocols.BaseProtocol.InputPipe
    def thermodynamic_state(self, value):
        pass

    def _find_block_size(self, charges, temperature, block_sizes_to_try=12, num_bootstrap=15):
        """Taken from https://github.com/MobleyLab/SMIRNOFF_paper_code/tree/master/FreeSolv"""

        block_size_grid = np.logspace(0, np.log10(len(self.trajectory)), block_sizes_to_try).astype('int')
        # The float -> int conversion sometimes leads to duplicate values, so avoid this
        block_size_grid = np.unique(block_size_grid)

        epsilon_grid = np.array([self._bootstrap(charges,
                                                 temperature,
                                                 block_length,
                                                 num_bootstrap) for block_length in block_size_grid])

        return block_size_grid[epsilon_grid.argmax()]

    def _bootstrap(self, charges, temperature, block_length, num_bootstrap):
        """Taken from https://github.com/MobleyLab/SMIRNOFF_paper_code/tree/master/FreeSolv"""

        bootstrap = arch.bootstrap.CircularBlockBootstrap(block_length, trajectory=self.trajectory)

        def bootstrap_func(trajectory):
            return mdtraj.geometry.static_dielectric(trajectory, charges, temperature)

        results = bootstrap.apply(bootstrap_func, num_bootstrap)
        epsilon_err = results.std()

        return epsilon_err

    def execute(self, directory):

        logging.info('Extracting dielectrics: ' + directory)

        base_exception = super(ExtractAverageDielectric, self).execute(directory)

        if isinstance(base_exception, ExtractAverageDielectric):
            return base_exception

        charge_list = []

        for force_index in range(self._system.getNumForces()):

            force = self._system.getForce(force_index)

            if not isinstance(force, openmm.NonbondedForce):
                continue

            for atom_index in range(force.getNumParticles()):

                charge = force.getParticleParameters(atom_index)[0]
                charge /= unit.elementary_charge

                charge_list.append(charge)

        temperature = self._thermodynamic_state.temperature / unit.kelvin

        # TODO: Pull out only equilibrated data.
        block_length = self._find_block_size(charge_list, temperature)
        dielectric_sigma = self._bootstrap(charge_list, temperature, block_length, block_length)

        dielectric = mdtraj.geometry.static_dielectric(self.trajectory, charge_list, temperature)

        self._value = unit.Quantity(dielectric, None)
        self._uncertainty = unit.Quantity(dielectric_sigma, None)

        logging.info('Extracted dielectrics: ' + directory)

        return self._get_output_dictionary()


# =============================================================================================
# Dielectric Constant
# =============================================================================================

@register_estimable_property()
@register_thermoml_property(thermoml_string='Relative permittivity at zero frequency')
class DielectricConstant(PhysicalProperty):
    """A class representation of a dielectric property"""

    @staticmethod
    def get_default_calculation_schema():
        
        schema = CalculationSchema(property_type=DielectricConstant.__name__)
        schema.id = '{}{}'.format(DielectricConstant.__name__, 'Schema')

        # Initial coordinate and topology setup.
        build_coordinates = protocols.BuildCoordinatesPackmol()
        build_coordinates.id = 'build_coordinates'

        build_coordinates.substance = ProtocolPath('substance', 'global')

        schema.protocols[build_coordinates.id] = build_coordinates.schema

        assign_topology = protocols.BuildSmirnoffTopology()
        assign_topology.id = 'build_topology'

        assign_topology.force_field_path = ProtocolPath('force_field_path', 'global')

        assign_topology.coordinate_file = ProtocolPath('coordinate_file', build_coordinates.id)
        assign_topology.molecules = ProtocolPath('molecules', build_coordinates.id)

        schema.protocols[assign_topology.id] = assign_topology.schema

        # Equilibration
        energy_minimisation = protocols.RunEnergyMinimisation()
        energy_minimisation.id = 'energy_minimisation'

        energy_minimisation.input_coordinate_file = ProtocolPath('coordinate_file', build_coordinates.id)
        energy_minimisation.system = ProtocolPath('system', assign_topology.id)

        schema.protocols[energy_minimisation.id] = energy_minimisation.schema

        npt_equilibration = protocols.RunOpenMMSimulation()
        npt_equilibration.id = 'npt_equilibration'

        npt_equilibration.ensemble = protocols.RunOpenMMSimulation.Ensemble.NPT

        npt_equilibration.steps = 2  # Debug settings.
        npt_equilibration.output_frequency = 1  # Debug settings.

        npt_equilibration.thermodynamic_state = ProtocolPath('thermodynamic_state', 'global')

        npt_equilibration.input_coordinate_file = ProtocolPath('output_coordinate_file', energy_minimisation.id)
        npt_equilibration.system = ProtocolPath('system', assign_topology.id)

        schema.protocols[npt_equilibration.id] = npt_equilibration.schema

        # Production

        npt_production = protocols.RunOpenMMSimulation()
        npt_production.id = 'npt_production'

        npt_production.ensemble = protocols.RunOpenMMSimulation.Ensemble.NPT

        npt_production.steps = 200  # Debug settings.
        npt_production.output_frequency = 20  # Debug settings.

        npt_production.thermodynamic_state = ProtocolPath('thermodynamic_state', 'global')

        npt_production.input_coordinate_file = ProtocolPath('output_coordinate_file', npt_equilibration.id)
        npt_production.system = ProtocolPath('system', assign_topology.id)

        schema.protocols[npt_production.id] = npt_production.schema

        # Analysis
        extract_dielectric = ExtractAverageDielectric()
        extract_dielectric.id = 'extract_dielectric'

        extract_dielectric.thermodynamic_state = ProtocolPath('thermodynamic_state', 'global')

        extract_dielectric.input_coordinate_file = ProtocolPath('output_coordinate_file', npt_production.id)
        extract_dielectric.trajectory_path = ProtocolPath('trajectory', npt_production.id)
        extract_dielectric.system = ProtocolPath('system', assign_topology.id)

        schema.protocols[extract_dielectric.id] = extract_dielectric.schema

        # Set up a conditional group to ensure convergence of uncertainty
        converge_uncertainty = groups.ConditionalGroup([
            npt_production.id,
            extract_dielectric.id
        ])
        converge_uncertainty.id = 'converge_uncertainty'

        # converge_uncertainty.input_references = [
        #     # # Locals
        #     # ProtocolInputReference(input_property_name='left_hand_value',
        #     #                        output_protocol_id='extract_dielectric',
        #     #                        output_property_name='uncertainty'),
        #     # # Globals
        #     # ProtocolInputReference(input_property_name='right_hand_value',
        #     #                        output_protocol_id='global',
        #     #                        output_property_name='uncertainty'),
        # ]

        schema.groups[converge_uncertainty.id] = converge_uncertainty.schema

        # Define where the final values come from.
        schema.final_value_source = ProtocolPath('value', extract_dielectric.id)
        schema.final_uncertainty_source = ProtocolPath('uncertainty', extract_dielectric.id)

        return schema
