#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================

"""
Protocol API.

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

import mdtraj

import numpy as np

from os import path

from enum import Enum, IntFlag, unique

from pymbar import timeseries

from openeye import oechem, oeomega

from openforcefield.utils import packmol
from openforcefield.typing.engines import smirnoff

from simtk import openmm, unit
from simtk.openmm import app


# =============================================================================================
# Protocols
# =============================================================================================

class ProtocolData:
    """
    Stores all of the information which can be passed to, and between
    calculation protocols.
    """

    def __init__(self):

        self.substance_tag = ''

        self.root_directory = ''

        self.molecules = None
        self.force_field = None

        self.positions = None
        self.topology = None

        self.trajectory_path = None
        self.statistics_path = None

        self.system = None

    @classmethod
    def clone(cls, existing_instance):

        return_value = cls()

        return_value.substance_tag = existing_instance.substance_tag
        return_value.root_directory = existing_instance.root_directory

        return_value.molecules = existing_instance.molecules
        return_value.force_field = existing_instance.force_field

        return_value.positions = copy.deepcopy(existing_instance.positions)
        return_value.topology = existing_instance.topology

        return_value.trajectory_path = existing_instance.trajectory_path
        return_value.statistics_path = existing_instance.statistics_path

        return_value.system = copy.deepcopy(existing_instance.system)

        return return_value


class Protocol:
    """
    The base class for a protocol which would form one
    step of a property calculation.

    A protocol may for example:

        create the coordiantes of a mixed simulation box
        set up a bound ligand-protein system
        build the simulation topology
        perform an energy minimisation

    Protocols may be chained together, this modularly defining
    a larger property calculation.

    """

    def set_measured_property(self, measured_property):
        pass

    def execute(self, protocol_data):
        """
        Allow protocols to be daisy chained together by passing the output
        of the previous protocol (coordinates + topol + stats?) to the next
        in line.
        """

        # Return the results of this protocol, ready to pass down the line.
        return None

    @classmethod
    def from_xml(cls, xml_node):
        raise NotImplementedError()

    def compare_to(self, protocol):
        return type(self) == type(protocol)


class BuildLiquidCoordinates(Protocol):

    _cached_molecules = {}

    # TODO: Determine the maximum number of molecules automatically
    def __init__(self):
        """
            Parameters
            ----------
            max_molecules : int, optional, default=True
                The maxmimum number of molecules in the system to be created.
            mass_density : float, simtk.unit.Quantity, or None; optional, default=None
                If provided, will aid in the selecting an initial box size.
        """

        self._substance = None

        self.max_molecules = 100
        self.mass_density = 1.0 * unit.grams / unit.milliliters

    # TODO: Replace with the toolkit function when finished.
    def _create_molecule(self, smiles):
        """
        Create molecule from a smiles pattern.

         Parameters
        ----------
        smiles : str
            Smiles pattern
         Returns
        -------
        molecule : OEMol
            OEMol with 3D coordinates, but no charges
         """

        # Check cache
        if smiles in self._cached_molecules:
            return copy.deepcopy(self._cached_molecules[smiles])

        # Create molecule from smiles.
        molecule = oechem.OEMol()
        parse_smiles_options = oechem.OEParseSmilesOptions(quiet=True)

        if not oechem.OEParseSmiles(molecule, smiles, parse_smiles_options):

            logging.warning('Could not parse SMILES: ' + smiles)
            return None

        # Normalize molecule
        oechem.OEAssignAromaticFlags(molecule, oechem.OEAroModelOpenEye)
        oechem.OEAddExplicitHydrogens(molecule)
        oechem.OETriposAtomNames(molecule)

        # Create configuration
        omega = oeomega.OEOmega()

        omega.SetMaxConfs(1)
        omega.SetIncludeInput(False)
        omega.SetCanonOrder(False)
        omega.SetSampleHydrogens(True)
        omega.SetStrictStereo(True)
        omega.SetStrictAtomTypes(False)

        status = omega(molecule)

        if not status:

            logging.warning('Could not generate a conformer for ' + smiles)
            return None

        self._cached_molecules[smiles] = molecule

        return molecule

    def set_measured_property(self, measured_property):
        self._substance = measured_property.substance

    def execute(self, protocol_data):

        logging.info('Generating coordinates: ' + protocol_data.root_directory)

        if self._substance is None:

            logging.warning('The BuildLiquidCoordinatesProtocol requires a Mixture as'
                            'input.')

            return None

        if protocol_data.positions is not None and \
           protocol_data.topology is not None:

            # The positions have already been built
            return protocol_data

        molecules = []

        for component in self._substance.components:

            molecule = self._create_molecule(component.smiles)

            if molecule is None:
                return None

            molecules.append(molecule)

        # Determine how many molecules of each type will be present in the system.
        mole_fractions = np.array([component.mole_fraction for component in self._substance.components])

        n_copies = np.random.multinomial(self.max_molecules - self._substance.number_of_impurities,
                                         pvals=mole_fractions)

        # Each impurity must have exactly one molecule
        for (index, component) in enumerate(self._substance.components):

            if component.impurity:
                n_copies[index] = 1

        # Create packed box
        topology, positions = packmol.pack_box(molecules, n_copies, mass_density=self.mass_density)

        if topology is None or positions is None:
            return None

        protocol_data.molecules = molecules

        protocol_data.positions = positions
        protocol_data.topology = topology

        with open(path.join(protocol_data.root_directory, 'output.pdb'), 'w+') as minimised_file:
            app.PDBFile.writeFile(topology, positions, minimised_file)

        logging.info('Coordinates generated: ' + self._substance.to_tag())

        return protocol_data

    @classmethod
    def from_xml(cls, xml_node):

        return_value = cls()

        max_molecules_node = xml_node.find('max_molecules')

        if max_molecules_node is not None:
            return_value.max_molecules = int(max_molecules_node.text)

        mass_density_node = xml_node.find('mass_density')

        if mass_density_node is not None:
            return_value.mass_density = float(mass_density_node.text) * unit.grams / unit.milliliters

        return return_value

    def compare_to(self, protocol):

        return super(BuildLiquidCoordinates, self).compare_to(protocol) and \
               self.max_molecules == protocol.max_molecules and \
               self.mass_density == protocol.mass_density


class BuildSmirnoffTopology(Protocol):

    def execute(self, protocol_data):

        logging.info('Generating topology: ' + protocol_data.root_directory)

        system = protocol_data.force_field.createSystem(protocol_data.topology,
                                                        protocol_data.molecules,
                                                        nonbondedMethod=smirnoff.PME,
                                                        chargeMethod='OECharges_AM1BCCSym')

        if system is None:

            logging.warning('Failed to create a system from the'
                            'provided topology and molecules')

            return None

        protocol_data.system = system

        logging.info('Topology generated: ' + protocol_data.root_directory)

        return protocol_data

    @classmethod
    def from_xml(cls, xml_node):

        return_value = cls()
        return return_value

    def compare_to(self, protocol):
        return super(BuildSmirnoffTopology, self).compare_to(protocol)


class RunEnergyMinimisation(Protocol):

    def __init__(self):

        # TODO: Add arguments for max iter + tolerance
        pass

    def execute(self, protocol_data):

        logging.info('Minimising energy: ' + protocol_data.root_directory)

        integrator = openmm.VerletIntegrator(0.002 * unit.picoseconds)

        simulation = app.Simulation(protocol_data.topology,
                                    protocol_data.system, integrator)

        simulation.context.setPositions(protocol_data.positions)

        simulation.minimizeEnergy()

        positions = simulation.context.getState(getPositions=True).getPositions()

        with open(path.join(protocol_data.root_directory, 'minimised.pdb'), 'w+') as minimised_file:
            app.PDBFile.writeFile(simulation.topology, positions, minimised_file)

        protocol_data.positions = positions

        logging.info('Energy minimised: ' + protocol_data.root_directory)

        return protocol_data

    @classmethod
    def from_xml(cls, xml_node):

        return_value = cls()
        return return_value

    def compare_to(self, protocol):

        # TODO: Properly implement comparison
        return super(RunEnergyMinimisation, self).compare_to(protocol)


class RunOpenMMSimulation(Protocol):

    class Ensemble(Enum):

        NVT = 0
        NPT = 1

    def __init__(self):

        self.thermodynamic_state = None

        self.steps = 1000

        self.thermostat_friction = 1.0 / unit.picoseconds
        self.timestep = 0.002 * unit.picoseconds

        self.output_frequency = 1000

        self.ensemble = self.Ensemble.NPT

    def set_measured_property(self, measured_property):
        self.thermodynamic_state = measured_property.thermodynamic_state

    def execute(self, protocol_data):

        temperature = self.thermodynamic_state.temperature
        pressure = self.thermodynamic_state.pressure

        substance_tag = protocol_data.root_directory

        if temperature is None:
            logging.error('A temperature must be set to perform a simulation in any ensemble: ' + substance_tag)
            return None
        if self.ensemble is self.Ensemble.NPT and pressure is None:
            logging.error('A pressure must be set to perform an NPT simulation: ' + substance_tag)
            return None

        logging.info('Performing a simulation in the ' + str(self.ensemble) + ' ensemble: ' + substance_tag)

        # For now set some 'best guess' thermostat parameters.
        integrator = openmm.LangevinIntegrator(temperature,
                                               self.thermostat_friction,
                                               self.timestep)

        system = protocol_data.system

        if self.ensemble is self.Ensemble.NPT:
            barostat = openmm.MonteCarloBarostat(pressure, temperature)

            system = copy.deepcopy(system)
            system.addForce(barostat)

        simulation = app.Simulation(protocol_data.topology, system, integrator)
        simulation.context.setPositions(protocol_data.positions)

        simulation.context.setVelocitiesToTemperature(temperature)

        trajectory_path = path.join(protocol_data.root_directory, 'trajectory.dcd')
        statistics_path = path.join(protocol_data.root_directory, 'statistics.dat')

        configuration_path = path.join(protocol_data.root_directory, 'input.pdb')

        with open(configuration_path, 'w+') as configuration_file:
            app.PDBFile.writeFile(protocol_data.topology, protocol_data.positions, configuration_file)

        simulation.reporters.append(app.DCDReporter(trajectory_path, self.output_frequency))

        simulation.reporters.append(app.StateDataReporter(statistics_path, self.output_frequency, step=True,
                                                          potentialEnergy=True, temperature=True, volume=True))

        try:
            simulation.step(self.steps)
        except Exception:
            logging.warning('Failed to run in ' + protocol_data.root_directory)
            return None

        positions = simulation.context.getState(getPositions=True).getPositions()

        protocol_data.positions = positions

        protocol_data.trajectory_path = trajectory_path
        protocol_data.statistics_path = statistics_path

        logging.info('Simulation performed in the ' + str(self.ensemble) + ' ensemble: ' + substance_tag)

        configuration_path = path.join(protocol_data.root_directory, 'output.pdb')

        with open(configuration_path, 'w+') as configuration_file:
            app.PDBFile.writeFile(protocol_data.topology, protocol_data.positions, configuration_file)

        return protocol_data

    @classmethod
    def from_xml(cls, xml_node):

        return_value = cls()

        steps_node = xml_node.find('steps')

        if steps_node is not None:
            return_value.steps = int(steps_node.text)

        thermostat_friction_node = xml_node.find('thermostat_friction')

        if thermostat_friction_node is not None:
            return_value.thermostat_friction = float(thermostat_friction_node.text) / unit.picoseconds

        timestep_node = xml_node.find('timestep')

        if timestep_node is not None:
            return_value.timestep = float(timestep_node.text) * unit.picoseconds

        output_frequency_node = xml_node.find('output_frequency')

        if output_frequency_node is not None:
            return_value.output_frequency = int(output_frequency_node.text)

        ensemble_node = xml_node.find('ensemble')

        if ensemble_node is not None and ensemble_node.text in cls.Ensemble:
                return_value.ensemble = cls.Ensemble[ensemble_node.text]

        return return_value

    def compare_to(self, protocol):

        return super(RunOpenMMSimulation, self).compare_to(protocol) and \
               self.thermodynamic_state.temperature == protocol.thermodynamic_state.temperature and \
               self.thermodynamic_state.pressure == protocol.thermodynamic_state.pressure and \
               self.ensemble == protocol.ensemble


class AveragePropertyProtocol(Protocol):

    def __init__(self):

        self.thermodynamic_state = None
        self.substance = None

        self.value = None
        self.uncertainty = None

    def set_measured_property(self, measured_property):

        self.substance = measured_property.substance
        self.thermodynamic_state = measured_property.thermodynamic_state

    def execute(self, protocol_data):

        return protocol_data

    @staticmethod
    def calculate_average_and_error(correlated_data):

        # Compute the indices of the uncorrelated timeseries
        [equilibration_index, inefficiency, Neff_max] = timeseries.detectEquilibration(correlated_data)
        equilibrated_data = correlated_data[equilibration_index:]

        # Extract a set of uncorrelated data points.
        indices = timeseries.subsampleCorrelatedData(equilibrated_data, g=inefficiency)
        uncorrelated_data = equilibrated_data[indices]

        average = uncorrelated_data.mean()
        uncertainty = uncorrelated_data.std() * len(uncorrelated_data) ** -0.5

        return average, uncertainty

    @classmethod
    def from_xml(cls, xml_node):
        raise NotImplementedError('AveragePropertyProtocol is an abstract class.')


class AverageTrajectoryProperty(AveragePropertyProtocol):

    def __init__(self):
        super().__init__()

        self.trajectory = None

    def execute(self, protocol_data):

        if protocol_data.trajectory_path is None:
            logging.warning('The AverageTrajectoryProperty protocol '
                            'requires a previously calculated trajectory')

            return None

        configuration_path = path.join(protocol_data.root_directory, 'configuration.pdb')

        with open(configuration_path, 'w+') as output_file:
            app.PDBFile.writeFile(protocol_data.topology, protocol_data.positions, output_file)

        self.trajectory = mdtraj.load_dcd(filename=protocol_data.trajectory_path, top=configuration_path)

        return protocol_data

    @classmethod
    def from_xml(cls, xml_node):
        raise NotImplementedError('AverageTrajectoryProperty is an abstract class.')


class ExtractAverageDensity(AverageTrajectoryProperty):

    def __init__(self):
        super().__init__()

    def execute(self, protocol_data):

        logging.info('Extracting densities: ' + protocol_data.root_directory)

        if super(ExtractAverageDensity, self).execute(protocol_data) is None:
            return None

        mass_list = []

        for atom_index in range(protocol_data.system.getNumParticles()):

            mass = protocol_data.system.getParticleMass(atom_index)
            mass /= (unit.gram / unit.mole)

            mass_list.append(mass)

        densities = mdtraj.density(self.trajectory, mass_list)

        self.value, self.uncertainty = self.calculate_average_and_error(densities)

        self.value *= unit.kilogram * unit.meter ** -3
        self.uncertainty *= unit.kilogram * unit.meter ** -3

        logging.info('Extracted densities: ' + protocol_data.root_directory)

        return protocol_data

    @classmethod
    def from_xml(cls, xml_node):

        return_value = cls()
        return return_value


class ExtractAverageDielectric(AverageTrajectoryProperty):

    def __init__(self):
        super().__init__()

    def execute(self, protocol_data):

        logging.info('Extracting dielectrics: ' + protocol_data.root_directory)

        if super(ExtractAverageDielectric, self).execute(protocol_data) is None:
            return None

        charge_list = []

        for force_index in range(protocol_data.system.getNumForces()):

            force = protocol_data.system.getForce(force_index)

            if not isinstance(force, openmm.NonbondedForce):
                continue

            for atom_index in range(force.getNumParticles()):

                charge = force.getParticleParameters(atom_index)[0]
                charge /= unit.elementary_charge

                charge_list.append(charge)

        temperature = self.thermodynamic_state.temperature / unit.kelvin

        # Determine the frame index at which the dipole moment has
        # reached equilibrium.
        # dipoles = mdtraj.geometry.dipole_moments(self.trajectory, charge_list)
        #
        # [equilibration_index, inefficiency, Neff_max] = timeseries.detectEquilibration(dipoles)
        # equilibrated_trajectory = self.trajectory[equilibration_index:]

        dielectric = mdtraj.geometry.static_dielectric(self.trajectory, charge_list, temperature)

        # TODO: Calculate uncertainty in dielectric constant
        self.value = dielectric
        self.uncertainty = 0.0

        logging.info('Extracted dielectrics: ' + protocol_data.root_directory)

        return protocol_data

    @classmethod
    def from_xml(cls, xml_node):

        return_value = cls()
        return return_value


