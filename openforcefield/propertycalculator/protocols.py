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

import copy
import logging

from enum import IntFlag

import numpy as np

from openeye import oechem, oeomega

from openforcefield import packmol
from openforcefield.propertycalculator.client import CalculatedPhysicalProperty

from openforcefield.typing.engines import smirnoff

from simtk import openmm, unit
from simtk.openmm import app


# =============================================================================================
# Protocols
# =============================================================================================

class ExtractableStatistics(IntFlag):
    """
    An array describing which properties can be extracted from
    a standard molecular simulations.
    """

    Undefined = 0x00
    Density = 0x02
    DielectricConstant = 0x03


class ProtocolData:
    """
    Stores all of the information which can be passed to, and between
    calculation protocols.
    """

    def __init__(self):

        self.root_directory = ''

        self.substance = None
        self.molecules = None

        self.force_field = None

        self.thermodynamic_state = None

        self.positions = None
        self.topology = None

        self.system = None

        self.results_directory = None

        self.extracted_statistics = None


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
    def execute(self, protocol_data):
        """
        Allow protocols to be daisy chained together by passing the output
        of the previous protocol (coordinates + topol + stats?) to the next
        in line.
        """

        # Return the results of this protocol, ready to pass down the line.
        return None


class BuildLiquidCoordinates(Protocol):

    _cached_molecules = {}

    # TODO: Determine the maximum number of molecules automatically
    def __init__(self, max_molecules=1000, mass_density=1.0 * unit.grams / unit.milliliters):
        """
            Parameters
            ----------
            max_molecules : int, optional, default=True
                The maxmimum number of molecules in the system to be created.
            mass_density : float, simtk.unit.Quantity, or None; optional, default=None
                If provided, will aid in the selecting an initial box size.
        """

        self._max_molecules = max_molecules
        self._mass_density = mass_density

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

    def execute(self, protocol_data):

        logging.info('Generating coordinates: ' + protocol_data.substance.to_tag())

        if protocol_data.substance is None:

            logging.warning('The BuildLiquidCoordinatesProtocol requires a Mixture as'
                            'input.')

            return None

        if protocol_data.positions is not None and \
           protocol_data.topology is not None:

            # The positions have already been built
            return protocol_data

        molecules = []

        for component in protocol_data.substance.components:

            molecule = self._create_molecule(component.smiles)

            if molecule is None:
                return None

            molecules.append(molecule)

        # Determine how many molecules of each type will be present in the system.
        mole_fractions = np.array([component.mole_fraction for component in protocol_data.substance.components])

        n_copies = np.random.multinomial(self._max_molecules - protocol_data.substance.number_of_impurities,
                                         pvals=mole_fractions)

        # Each impurity must have exactly one molecule
        for (index, component) in enumerate(protocol_data.substance.components):

            if component.impurity:
                n_copies[index] = 1

        # Create packed box
        topology, positions = packmol.pack_box(molecules, n_copies, mass_density=self._mass_density)

        if topology is None or positions is None:
            return None

        protocol_data.molecules = molecules

        protocol_data.positions = positions
        protocol_data.topology = topology

        logging.info('Coordinates generated: ' + protocol_data.substance.to_tag())

        return protocol_data


class BuildSmirnoffTopology(Protocol):

    def execute(self, protocol_data: ProtocolData):

        logging.info('Generating topology: ' + protocol_data.substance.to_tag())

        system = protocol_data.force_field.createSystem(protocol_data.topology,
                                                        protocol_data.molecules,
                                                        nonbondedMethod=smirnoff.PME,
                                                        chargeMethod='BCC')

        if system is None:

            logging.warning('Failed to create a system from the'
                            'provided topology and molecules')

            return None

        protocol_data.system = system

        logging.info('Topology generated: ' + protocol_data.substance.to_tag())

        return protocol_data


class RunEnergyMinimisation(Protocol):

    def __init__(self):

        # TODO: Add arguments for max iter + tolerance
        pass

    def execute(self, protocol_data):

        substance_tag = protocol_data.substance.to_tag()
        logging.info('Minimising energy: ' + substance_tag)

        integrator = openmm.VerletIntegrator(0.002 * unit.picoseconds)

        simulation = app.Simulation(protocol_data.topology,
                                    protocol_data.system, integrator)

        simulation.context.setPositions(protocol_data.positions)

        with open(protocol_data.root_directory + substance_tag + '_PRE_EM.pdb', 'w+') as minimised_file:
            app.PDBFile.writeFile(simulation.topology, protocol_data.positions, minimised_file)

        simulation.minimizeEnergy()

        positions = simulation.context.getState(getPositions=True).getPositions()

        with open(protocol_data.root_directory + substance_tag + '_EM.pdb', 'w+') as minimised_file:
            app.PDBFile.writeFile(simulation.topology, positions, minimised_file)

        protocol_data.positions = positions

        substance_tag = protocol_data.substance.to_tag()
        logging.info('Energy minimised: ' + substance_tag)

        return protocol_data


class RunNVTSimulation(Protocol):

    def __init__(self):

        # TODO: Add parameters for steps, timestep etc..
        pass

    def execute(self, protocol_data):

        temperature = protocol_data.thermodynamic_state.temperature
        substance_tag = protocol_data.substance.to_tag()

        if temperature is None:
            logging.warning('A temperature must be set to perform an NVT simulation: ' + substance_tag)
            return protocol_data

        logging.info('Performing NVT simulation: ' + substance_tag)

        # For now set some 'best guess' thermostat parameters.
        integrator = openmm.LangevinIntegrator(temperature,
                                               1.000 / unit.picosecond,
                                               0.002 * unit.picoseconds)

        simulation = app.Simulation(protocol_data.topology, protocol_data.system, integrator)
        simulation.context.setPositions(protocol_data.positions)

        simulation.context.setVelocitiesToTemperature(temperature)

        trajectory_path = protocol_data.root_directory + substance_tag + '_NVT.dcd'
        statistics_path = protocol_data.root_directory + substance_tag + '_NVT.dat'

        simulation.reporters.append(app.DCDReporter(trajectory_path, 1000))

        simulation.reporters.append(app.StateDataReporter(statistics_path, 1000, step=True,
                                                          potentialEnergy=True, temperature=True))

        simulation.step(10000)

        positions = simulation.context.getState(getPositions=True).getPositions()
        protocol_data.positions = positions

        logging.info('NVT simulation performed: ' + substance_tag)

        return protocol_data


class RunNPTSimulation(Protocol):

    def __init__(self):

        # TODO: Add parameters for steps, timestep etc..
        pass

    def execute(self, protocol_data):

        temperature = protocol_data.thermodynamic_state.temperature
        pressure = protocol_data.thermodynamic_state.pressure

        substance_tag = protocol_data.substance.to_tag()

        if temperature is None:
            logging.warning('A temperature must be set to perform an NPT simulation: ' + substance_tag)
            return protocol_data
        if pressure is None:
            logging.warning('A pressure must be set to perform an NPT simulation: ' + substance_tag)
            return protocol_data

        logging.info('Performing NPT simulation: ' + substance_tag)

        # For now set some 'best guess' thermostat parameters.
        integrator = openmm.LangevinIntegrator(temperature,
                                               1.000 / unit.picosecond,
                                               0.002 * unit.picoseconds)

        barostat = openmm.MonteCarloBarostat(pressure, temperature)

        cloned_system = copy.deepcopy(protocol_data.system)
        cloned_system.addForce(barostat)

        simulation = app.Simulation(protocol_data.topology, cloned_system, integrator)
        simulation.context.setPositions(protocol_data.positions)

        simulation.context.setVelocitiesToTemperature(temperature)

        trajectory_path = protocol_data.root_directory + substance_tag + '_NPT.dcd'
        statistics_path = protocol_data.root_directory + substance_tag + '_NPT.dat'

        simulation.reporters.append(app.DCDReporter(trajectory_path, 1000))

        simulation.reporters.append(app.StateDataReporter(statistics_path, 1000, step=True,
                                                          potentialEnergy=True, temperature=True))

        simulation.step(100000)

        positions = simulation.context.getState(getPositions=True).getPositions()
        protocol_data.positions = positions

        logging.info('NPT simulation performed: ' + substance_tag)

        return protocol_data


class ExtractStatistics(Protocol):

    def __init__(self, types):

        self.types = types

    def execute(self, protocol_data: ProtocolData):

        logging.info('Extracting statistics.')

        if protocol_data.extracted_statistics is None:
            protocol_data.extracted_statistics = []

        protocol_data.extracted_statistics.append(0.0)

        logging.info('Extracted statistics.')

        return protocol_data
