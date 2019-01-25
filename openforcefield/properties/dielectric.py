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
import numpy as np
from simtk import openmm, unit
from simtk.openmm import System

from openforcefield.properties.datasets import register_thermoml_property
from openforcefield.properties.estimator import CalculationSchema, register_estimable_property
from openforcefield.properties.estimator.workflow import protocols, groups, protocol_input
from openforcefield.properties.estimator.workflow.protocols import AverageTrajectoryProperty, \
    register_calculation_protocol, ProtocolPath
from openforcefield.properties.properties import PhysicalProperty
from openforcefield.properties.thermodynamics import ThermodynamicState
from openforcefield.utils import statistics


# =============================================================================================
# Custom Protocol Building Blocks
# =============================================================================================

@register_calculation_protocol()
class ExtractAverageDielectric(AverageTrajectoryProperty):
    """Extracts the average dielectric constant from a simulation trajectory.
    """
    def __init__(self, protocol_id):
        super().__init__(protocol_id)

        self._system = None
        self._thermodynamic_state = None

    @protocol_input(System)
    def system(self, value):
        """The system object which defines the forces present in the system."""
        pass

    @protocol_input(ThermodynamicState)
    def thermodynamic_state(self, value):
        """The thermodynamic state at which the trajectory was generated."""
        pass

    def _bootstrap_function(self, sample_data):
        """Calculates the static dielectric constant from an
        array of dipoles and volumes.

        Notes
        -----
        The static dielectric constant is taken from for Equation 7 of [1]

        References
        ----------
        [1] A. Glattli, X. Daura and W. F. van Gunsteren. Derivation of an improved simple point charge
            model for liquid water: SPC/A and SPC/L. J. Chem. Phys. 116(22):9811-9828, 2002

        Parameters
        ----------
        sample_data: np.ndarray, shape=(num_frames, 4), dtype=float
            The dataset to bootstap. The data stored by trajectory frame is
            four dimensional (Mx, My, Mz, V) where M is dipole moment and
            V is volume.

        Returns
        -------
        float
            The unitless static dielectric constant
        """

        temperature = self._thermodynamic_state.temperature

        dipoles = np.zeros([sample_data.shape[0], 3])
        volumes = np.zeros([sample_data.shape[0], 1])

        for index in range(sample_data.shape[0]):

            dipoles[index][0] = sample_data[index][0]
            dipoles[index][1] = sample_data[index][1]
            dipoles[index][2] = sample_data[index][2]

            volumes[index] = sample_data[index][3]

        dipole_mu = dipoles.mean(0)
        shifted_dipoles = dipoles - dipole_mu

        dipole_variance = (shifted_dipoles * shifted_dipoles).sum(-1).mean(0) * \
                          (unit.elementary_charge * unit.nanometers) ** 2

        volume = volumes.mean() * unit.nanometer**3

        e0 = 8.854187817E-12 * unit.farad / unit.meter  # Taken from QCElemental

        dielectric_constant = 1.0 + dipole_variance / (3 *
                                                       unit.BOLTZMANN_CONSTANT_kB *
                                                       temperature *
                                                       volume *
                                                       e0)

        return dielectric_constant

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

        dipole_moments = mdtraj.geometry.dipole_moments(self.trajectory, charge_list)
        volumes = self.trajectory.unitcell_volumes

        dipole_moments, self._equilibration_index, self._statistical_inefficiency = \
            statistics.decorrelate_time_series(dipole_moments)

        dipole_moments_and_volume = np.zeros([dipole_moments.shape[0], 4])

        for index in range(dipole_moments.shape[0]):

            dipole = dipole_moments[index]
            volume = volumes[index]

            dipole_moments_and_volume[index] = np.array([dipole[0], dipole[1], dipole[2], volume])

        self._value, self._uncertainty = self._perform_bootstrapping(dipole_moments_and_volume)

        self._value = unit.Quantity(self._value, None)
        self._uncertainty = unit.Quantity(self._uncertainty, None)

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
        build_coordinates = protocols.BuildCoordinatesPackmol('build_coordinates')

        build_coordinates.substance = ProtocolPath('substance', 'global')

        schema.protocols[build_coordinates.id] = build_coordinates.schema

        assign_topology = protocols.BuildSmirnoffTopology('build_topology')

        assign_topology.force_field_path = ProtocolPath('force_field_path', 'global')

        assign_topology.coordinate_file_path = ProtocolPath('coordinate_file_path', build_coordinates.id)
        assign_topology.substance = ProtocolPath('substance', 'global')

        schema.protocols[assign_topology.id] = assign_topology.schema

        # Equilibration
        energy_minimisation = protocols.RunEnergyMinimisation('energy_minimisation')

        energy_minimisation.input_coordinate_file = ProtocolPath('coordinate_file_path', build_coordinates.id)
        energy_minimisation.system = ProtocolPath('system', assign_topology.id)

        schema.protocols[energy_minimisation.id] = energy_minimisation.schema

        npt_equilibration = protocols.RunOpenMMSimulation('npt_equilibration')

        npt_equilibration.ensemble = protocols.RunOpenMMSimulation.Ensemble.NPT

        npt_equilibration.steps = 2  # Debug settings.
        npt_equilibration.output_frequency = 1  # Debug settings.

        npt_equilibration.thermodynamic_state = ProtocolPath('thermodynamic_state', 'global')

        npt_equilibration.input_coordinate_file = ProtocolPath('output_coordinate_file', energy_minimisation.id)
        npt_equilibration.system = ProtocolPath('system', assign_topology.id)

        schema.protocols[npt_equilibration.id] = npt_equilibration.schema

        # Production
        npt_production = protocols.RunOpenMMSimulation('npt_production')

        npt_production.ensemble = protocols.RunOpenMMSimulation.Ensemble.NPT

        npt_production.steps = 200  # Debug settings.
        npt_production.output_frequency = 20  # Debug settings.

        npt_production.thermodynamic_state = ProtocolPath('thermodynamic_state', 'global')

        npt_production.input_coordinate_file = ProtocolPath('output_coordinate_file', npt_equilibration.id)
        npt_production.system = ProtocolPath('system', assign_topology.id)

        # Analysis
        extract_dielectric = ExtractAverageDielectric('extract_dielectric')

        extract_dielectric.thermodynamic_state = ProtocolPath('thermodynamic_state', 'global')

        extract_dielectric.input_coordinate_file = ProtocolPath('output_coordinate_file', npt_production.id)
        extract_dielectric.trajectory_path = ProtocolPath('trajectory_file_path', npt_production.id)
        extract_dielectric.system = ProtocolPath('system', assign_topology.id)

        # Set up a conditional group to ensure convergence of uncertainty
        converge_uncertainty = groups.ConditionalGroup('converge_uncertainty')
        converge_uncertainty.add_protocols(npt_production, extract_dielectric)

        converge_uncertainty.left_hand_value = ProtocolPath('uncertainty',
                                                            converge_uncertainty.id,
                                                            extract_dielectric.id)

        converge_uncertainty.right_hand_value = ProtocolPath('target_uncertainty', 'global')

        converge_uncertainty.condition_type = groups.ConditionalGroup.ConditionType.LessThan

        converge_uncertainty.max_iterations = 1

        schema.protocols[converge_uncertainty.id] = converge_uncertainty.schema

        # Finally, extract uncorrelated data
        extract_uncorrelated_trajectory = protocols.ExtractUncorrelatedTrajectoryData('extract_traj')

        extract_uncorrelated_trajectory.statistical_inefficiency = ProtocolPath('statistical_inefficiency',
                                                                                converge_uncertainty.id,
                                                                                extract_dielectric.id)

        extract_uncorrelated_trajectory.equilibration_index = ProtocolPath('equilibration_index',
                                                                           converge_uncertainty.id,
                                                                           extract_dielectric.id)

        extract_uncorrelated_trajectory.input_coordinate_file = ProtocolPath('output_coordinate_file',
                                                                             converge_uncertainty.id,
                                                                             npt_production.id)

        extract_uncorrelated_trajectory.input_trajectory_path = ProtocolPath('trajectory_file_path',
                                                                             converge_uncertainty.id,
                                                                             npt_production.id)

        schema.protocols[extract_uncorrelated_trajectory.id] = extract_uncorrelated_trajectory.schema

        # Define where the final values come from.
        schema.final_value_source = ProtocolPath('value', converge_uncertainty.id, extract_dielectric.id)
        schema.final_uncertainty_source = ProtocolPath('uncertainty', converge_uncertainty.id, extract_dielectric.id)

        schema.final_coordinate_source = ProtocolPath('output_coordinate_file', converge_uncertainty.id,
                                                                                npt_production.id)

        schema.final_trajectory_source = ProtocolPath('output_trajectory_path', extract_uncorrelated_trajectory.id)

        return schema
