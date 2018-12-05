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

import mdtraj

import arch.bootstrap

import numpy as np

from os import path

from enum import Enum

from pymbar import timeseries

from openeye import oechem, oeomega

from openforcefield.utils import packmol
from openforcefield.utils.exceptions import XmlNodeMissingException
from openforcefield.typing.engines import smirnoff

from simtk import openmm, unit
from simtk.openmm import app


# =============================================================================================
# Protocols
# =============================================================================================

class ProtocolInputReference:
    """Stores a reference to a required input from another protocol.

    Each node represents a protocol to be executed.

    Parameters
    ----------
    protocol_id: str
        The identity of the other protocol.
    property_name: str
        The name of the property to inherit from the other protocol.
    """

    def __init__(self, input_property, protocol_id, property_name):

        self.input_property = input_property

        self.protocol_id = protocol_id
        self.property_name = property_name

    def __hash__(self):
        """returns the hash key of this ProtocolInputReference."""
        return hash((self.input_property, self.protocol_id, self.property_name))

    def __eq__(self, other):
        """Returns true if the two inputs are equal."""
        return (self.input_property == other.input_property and
                self.protocol_id == other.protocol_id and
                self.property_name == other.property_name)

    def __ne__(self, other):
        """Returns true if the two inputs are not equal."""
        return not (self == other)


class Protocol:
    """
    The base class for a protocol which would form one
    step of a property calculation.

    A protocol may for example:

        - create the coordinates of a mixed simulation box
        - set up a bound ligand-protein system
        - build the simulation topology
         - perform an energy minimisation

    Protocols may be chained together, thus defining
    a larger property calculation from simple building blocks.

    """

    class ProtocolPipe(object):
        """A custom decorator used to mark class attributes as either
         a required input, or output, of a protocol.

        Parameters
        ----------
        class_attribute: function
            The attribute to mark as a pipe.
        documentation: str, optional
            Documentation for this pipe.
        private_attribute: str, optional
            The name of the underlying private attribute. The default is
            '_attribute'
        """

        def __init__(self, class_attribute, documentation=None, private_attribute=None):
            if class_attribute.__doc__:
                documentation = class_attribute.__doc__

            self.__doc__ = documentation  # Set the documentation of the instance.

            if private_attribute is None:
                self.attribute = '_' + class_attribute.__name__  # Add leading underscore to the attribute name
            else:
                self.attribute = private_attribute

        def __get__(self, instance, owner=None):

            if not hasattr(instance, self.attribute):
                raise RuntimeError('Missing ' + self.attribute + 'attribute.')

            return getattr(instance, self.attribute)

        def __set__(self, instance, value):

            if not hasattr(instance, self.attribute):
                raise RuntimeError('Missing ' + self.attribute + 'attribute.')

            return setattr(instance, self.attribute, value)

    class InputPipe(ProtocolPipe):
        """A custom decorator used to mark properties as a required input to
        the protocol.

        Examples
        ----------
        To mark a property as an input pipe:

        >>> @Protocol.InputPipe
        >>> def substance(self):
        >>>     pass
        """
        def __init__(self, attribute, documentation=None):
            super().__init__(attribute, documentation)

    class OutputPipe(ProtocolPipe):
        """A custom decorator used to mark properties as an output of the
        the protocol.

        Examples
        ----------
        To mark a property as an output pipe:

        >>> @Protocol.OutputPipe
        >>> def positions(self):
        >>>     pass
        """
        def __init__(self, attribute, documentation=None):
            super().__init__(attribute, documentation)

    def __init__(self):

        # A unique identifier for this node.
        self.id = None

        self.required_inputs = self._find_types_with_decorator(Protocol.InputPipe)
        self.provided_outputs = self._find_types_with_decorator(Protocol.OutputPipe)

        # Defines where to pull the values from.
        self.input_references = []

    def _find_types_with_decorator(self, decorator_type):
        """ A method to collect all attributes marked by a specified
        decorator type (e.g. InputProperty).

        Parameters
        ----------
        decorator_type
            The type of decorator to search for.

        Returns
        ----------
        The names of the attributes decorated with the specified decorator.
        """
        inputs = [attribute_name for attribute_name in type(self).__dict__ if
                  type(type(self).__dict__[attribute_name]) is decorator_type]

        for base in type(self).__bases__:

            inputs.extend([attribute_name for attribute_name in base.__dict__ if
                           type(base.__dict__[attribute_name]) is decorator_type])

        return inputs

    def set_uuid(self, value):

        old_id_split = self.id.split('|')

        if len(old_id_split) == 1:

            self.id = value + '|' + old_id_split[0]

            for input_reference in self.input_references:

                if input_reference.protocol_id == 'global':
                    continue

                input_reference.protocol_id = value + '|' + input_reference.protocol_id

            return

        old_uuid = old_id_split[0]
        self.id = self.id.replace(old_uuid, value)

        for input_reference in self.input_references:

            if input_reference.protocol_id == 'global':
                continue

            input_reference.protocol_id = input_reference.protocol_id.replace(old_uuid, value)

    def execute(self, directory):
        """ Execute the protocol.

        Protocols may be chained together by passing the output
        of previous protocols as input to the current one.

        Parameters
        ----------
        directory : str
            The directory to store output data in.

        Returns
        ----------
        bool
            True if the command successfully executes.
        """

        # Return the results of this protocol, ready to pass down the line.
        return True

    @classmethod
    def from_xml(cls, xml_node):
        """ Creates a protocol from an xml definition.

        Parameters
        ----------
        xml_node : xml.etree.Element
            The element containing the xml to create the protocol from.

        Returns
        ----------
        Protocol
            The protocol created from the xml node.
        """
        return_value = cls()

        # Find the unique id of this protocol.
        id_node = xml_node.find('id')

        if id_node is None:
            raise XmlNodeMissingException('id')

        return_value.id = id_node.text

        inputs_node = xml_node.find('inputs')

        if inputs_node is not None:

            for input_node in inputs_node.findall('input'):

                if 'property' not in input_node.attrib:
                    raise Exception('Protocol inputs must define a property attribute.')

                input_property = input_node.attrib['property']

                text_split = input_node.text.split('.')

                if len(text_split) != 2:
                    raise Exception('Protocol inputs must be of the form: node_id.property_name')

                # Only process inputs which are actually required.
                if input_property not in return_value.required_inputs:
                    continue

                protocol_input = ProtocolInputReference(input_property,
                                                        text_split[0],
                                                        text_split[1])

                # Don't add multiple of the same input.
                if protocol_input in return_value.input_references:
                    continue

                return_value.input_references.append(protocol_input)

        # Make sure this protocol is being passed all the required inputs.
        for required_input in return_value.required_inputs:

            if not hasattr(return_value, required_input):

                raise Exception('A ' + type(return_value).__name__ + 'protocol must receive ' +
                                required_input + ' as an input.')

        return return_value

    def compare_to(self, other, id_maps):
        """ Compares this protocol with another.

        Parameters
        ----------
        other : Protocol
            The protocol to compare against.
        id_maps : dict(str, str)
            A dictionary that maps original protocol ids to the id of the protocol
            they have been merged into (the key and value may be the same).

        Returns
        ----------
        bool
            True if the protocols would essentialy perform the same task.
        """
        if not isinstance(self, type(other)):
            return False

        for protocol_input in self.input_references:

            self_protocol_id = protocol_input.protocol_id

            if self_protocol_id in id_maps:
                self_protocol_id = id_maps[self_protocol_id]

            shares_input = False

            for other_reference in other.input_references:

                other_protocol_id = other_reference.protocol_id

                if other_protocol_id in id_maps:
                    other_protocol_id = id_maps[other_protocol_id]

                if protocol_input.input_property != other_reference.input_property or \
                   self_protocol_id != other_protocol_id or \
                   protocol_input.property_name != other_reference.property_name:

                    continue

                shares_input = True
                break

            if shares_input is False:
                return False

            self_value = getattr(self, protocol_input.input_property)
            other_value = getattr(other, protocol_input.input_property)

            if self_value != other_value:
                return False

        return True


class BuildLiquidCoordinates(Protocol):
    """Create 3D coordinates and bond information for a given Substance

    The coordinates are created using packmol.

    Attributes
    ----------
    max_molecules : int, optional, default=True
        The maxmimum number of molecules in the system to be created.
    mass_density : float, simtk.unit.Quantity, or None; optional, default=None
        If provided, will aid in the selecting an initial box size.
    """

    _cached_molecules = {}

    def __init__(self):

        super().__init__()

        # inputs
        self._substance = None

        # outputs
        self._topology = None
        self._positions = None
        self._molecules = None

        # TODO: Determine the maximum number of molecules automatically
        self.max_molecules = 100
        self.mass_density = 1.0 * unit.grams / unit.milliliters

    @Protocol.InputPipe
    def substance(self):
        pass

    @Protocol.OutputPipe
    def topology(self):
        pass

    @Protocol.OutputPipe
    def positions(self):
        pass

    @Protocol.OutputPipe
    def molecules(self):
        pass

    def _create_molecule(self, smiles):
        """
        Create molecule from a smiles pattern.

        Todo
        ----------

        * Replace with the toolkit function when finished.

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
            return False

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
            return False

        self._cached_molecules[smiles] = molecule

        return molecule

    def execute(self, directory):

        logging.info('Generating coordinates: ' + directory)

        if self._substance is None:

            logging.warning('The BuildLiquidCoordinatesProtocol requires a '
                            'non-optional Substance as input.')

            return False

        molecules = []

        for component in self._substance.components:

            molecule = self._create_molecule(component.smiles)

            if molecule is None:
                return False

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
            return False

        self._molecules = molecules

        self._positions = positions
        self._topology = topology

        with open(path.join(directory, 'output.pdb'), 'w+') as minimised_file:
            app.PDBFile.writeFile(topology, positions, minimised_file)

        logging.info('Coordinates generated: ' + str(self._substance))

        return True

    @classmethod
    def from_xml(cls, xml_node):

        return_value = super(BuildLiquidCoordinates, cls).from_xml(xml_node)

        max_molecules_node = xml_node.find('max_molecules')

        if max_molecules_node is not None:
            return_value.max_molecules = int(max_molecules_node.text)

        mass_density_node = xml_node.find('mass_density')

        if mass_density_node is not None:
            return_value.mass_density = float(mass_density_node.text) * unit.grams / unit.milliliters

        return return_value

    def compare_to(self, protocol, id_maps):

        return super(BuildLiquidCoordinates, self).compare_to(protocol, id_maps) and \
               self.max_molecules == protocol.max_molecules and \
               self.mass_density == protocol.mass_density


class BuildSmirnoffTopology(Protocol):
    """Parametrise a set of molecules with a given smirnoff force field.
    """
    def __init__(self):

        super().__init__()

        # inputs
        self._force_field = None
        self._topology = None
        self._molecules = None
        # outputs
        self._system = None

    @Protocol.InputPipe
    def force_field(self, value):
        pass

    @Protocol.InputPipe
    def molecules(self, value):
        pass

    @Protocol.InputPipe
    def topology(self, value):
        pass

    @Protocol.OutputPipe
    def system(self):
        pass

    def execute(self, directory):

        logging.info('Generating topology: ' + directory)

        system = self._force_field.createSystem(self._topology,
                                                self._molecules,
                                                nonbondedMethod=smirnoff.PME,
                                                chargeMethod='OECharges_AM1BCCSym')

        if system is None:

            logging.warning('Failed to create a system from the'
                            'provided topology and molecules')

            return False

        self._system = system

        logging.info('Topology generated: ' + directory)

        return True


class RunEnergyMinimisation(Protocol):
    """Minimises the energy of a passed in system.
    """

    def __init__(self):

        super().__init__()

        # inputs
        self._topology = None
        self._positions = None
        self._system = None
        # outputs
        self._final_positions = None

        # TODO: Add arguments for max iter + tolerance
        pass

    @Protocol.InputPipe
    def topology(self, value):
        pass

    @Protocol.InputPipe
    def positions(self, value):
        pass

    @Protocol.InputPipe
    def system(self, value):
        pass

    @Protocol.OutputPipe
    def final_positions(self):
        return self._final_positions

    def execute(self, directory):

        logging.info('Minimising energy: ' + directory)

        integrator = openmm.VerletIntegrator(0.002 * unit.picoseconds)

        simulation = app.Simulation(self._topology,
                                    self._system, integrator)

        simulation.context.setPositions(self._positions)

        simulation.minimizeEnergy()

        positions = simulation.context.getState(getPositions=True).getPositions()

        with open(path.join(directory, 'minimised.pdb'), 'w+') as minimised_file:
            app.PDBFile.writeFile(simulation.topology, positions, minimised_file)

        self._final_positions = positions

        logging.info('Energy minimised: ' + directory)

        return True

    def compare_to(self, protocol, id_maps):

        # TODO: Properly implement comparison
        return super(RunEnergyMinimisation, self).compare_to(protocol, id_maps)


class RunOpenMMSimulation(Protocol):
    """Performs a molecular dynamics simulation in a given ensemble using OpenMM

    Attributes
    ----------
    steps : int
        The number of steps to run the simulation for
    timestep : float
        The timestep of the integrator.
    output_frequency : int
        The frequency with which to store simulation data.
    ensemble : RunOpenMMSimulation.Ensemble
        The ensemble to run the simulation in.
    """

    class Ensemble(Enum):
        """An enum describing the available ensembles.
        """
        NVT = 0
        NPT = 1

    def __init__(self):

        super().__init__()

        self.steps = 1000

        self.thermostat_friction = 1.0 / unit.picoseconds
        self.timestep = 0.002 * unit.picoseconds

        self.output_frequency = 1000

        self.ensemble = self.Ensemble.NPT

        # inputs
        self._thermodynamic_state = None
        self._topology = None
        self._positions = None
        self._system = None

        # outputs
        self._final_positions = None
        self._trajectory = None
        self._statistics = None

        # TODO: Add arguments for max iter + tolerance
        pass

    @Protocol.InputPipe
    def thermodynamic_state(self, value):
        pass

    @Protocol.InputPipe
    def topology(self, value):
        pass

    @Protocol.InputPipe
    def positions(self, value):
        pass

    @Protocol.InputPipe
    def system(self, value):
        pass

    @Protocol.OutputPipe
    def final_positions(self):
        pass

    @Protocol.OutputPipe
    def trajectory(self):
        pass

    @Protocol.OutputPipe
    def statistics(self):
        pass

    def execute(self, directory):

        temperature = self._thermodynamic_state.temperature
        pressure = self._thermodynamic_state.pressure

        if temperature is None:
            logging.error('A temperature must be set to perform a simulation in any ensemble: ' + directory)
            return False
        if self.ensemble is self.Ensemble.NPT and pressure is None:
            logging.error('A pressure must be set to perform an NPT simulation: ' + directory)
            return False

        logging.info('Performing a simulation in the ' + str(self.ensemble) + ' ensemble: ' + directory)

        # For now set some 'best guess' thermostat parameters.
        integrator = openmm.LangevinIntegrator(temperature,
                                               self.thermostat_friction,
                                               self.timestep)

        system = self._system

        if self.ensemble is self.Ensemble.NPT:
            barostat = openmm.MonteCarloBarostat(pressure, temperature)

            # inputs are READONLY! Never directly alter an input
            system = copy.deepcopy(system)
            system.addForce(barostat)

        simulation = app.Simulation(self._topology, system, integrator)
        simulation.context.setPositions(self._positions)

        simulation.context.setVelocitiesToTemperature(temperature)

        trajectory_path = path.join(directory, 'trajectory.dcd')
        statistics_path = path.join(directory, 'statistics.dat')

        configuration_path = path.join(directory, 'input.pdb')

        with open(configuration_path, 'w+') as configuration_file:

            app.PDBFile.writeFile(self._topology,
                                  self._positions, configuration_file)

        simulation.reporters.append(app.DCDReporter(trajectory_path, self.output_frequency))

        simulation.reporters.append(app.StateDataReporter(statistics_path, self.output_frequency, step=True,
                                                          potentialEnergy=True, temperature=True, volume=True))

        try:
            simulation.step(self.steps)
        except Exception:
            logging.warning('Failed to run in ' + directory)
            return False

        positions = simulation.context.getState(getPositions=True).getPositions()

        self._final_positions = positions

        self._trajectory = trajectory_path
        self._statistics = statistics_path

        logging.info('Simulation performed in the ' + str(self.ensemble) + ' ensemble: ' + directory)

        configuration_path = path.join(directory, 'output.pdb')

        with open(configuration_path, 'w+') as configuration_file:

            app.PDBFile.writeFile(self._topology,
                                  positions, configuration_file)

        return True

    @classmethod
    def from_xml(cls, xml_node):

        return_value = super(RunOpenMMSimulation, cls).from_xml(xml_node)

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

    def compare_to(self, protocol, id_maps):

        return super(RunOpenMMSimulation, self).compare_to(protocol, id_maps) and \
               self.ensemble == protocol.ensemble


class AveragePropertyProtocol(Protocol):
    """Calculates the average of a property and its uncertainty.
    """

    def __init__(self):

        super().__init__()

        self._value = None
        self._uncertainty = None

    @Protocol.OutputPipe
    def value(self):
        return self._value

    @Protocol.OutputPipe
    def uncertainty(self):
        return self._uncertainty

    def execute(self, directory):
        return True

    @staticmethod
    def calculate_average_and_error(correlated_data):
        """Calculates the average of a property and its uncertainty from
        a list of possibly correlated data.

        Parameters
        ----------
        correlated_data : list(float)
            The data to average over.

        Returns
        ----------
        float
            The average value
        float
            The uncertainty in the average.
        """

        # Compute the indices of the uncorrelated timeseries
        [equilibration_index, inefficiency, effictive_samples] = timeseries.detectEquilibration(correlated_data)
        equilibrated_data = correlated_data[equilibration_index:]

        # Extract a set of uncorrelated data points.
        indices = timeseries.subsampleCorrelatedData(equilibrated_data, g=inefficiency)
        uncorrelated_data = equilibrated_data[indices]

        average = uncorrelated_data.mean()
        uncertainty = uncorrelated_data.std() * len(uncorrelated_data) ** -0.5

        return average, uncertainty


class AverageTrajectoryProperty(AveragePropertyProtocol):
    """Calculates the average of a property from a simulation trajectory.
    """

    def __init__(self):

        super().__init__()

        self._topology = None
        self._positions = None
        self._trajectory_path = None

        self.trajectory = None

    @Protocol.InputPipe
    def topology(self, value):
        pass

    @Protocol.InputPipe
    def positions(self, value):
        pass

    @Protocol.InputPipe
    def trajectory_path(self, value):
        pass

    def execute(self, directory):

        if self._trajectory_path is None:

            logging.warning('The AverageTrajectoryProperty protocol '
                            'requires a previously calculated trajectory')

            return False

        configuration_path = path.join(directory, 'configuration.pdb')

        with open(configuration_path, 'w+') as output_file:
            app.PDBFile.writeFile(self._topology, self._positions, output_file)

        self.trajectory = mdtraj.load_dcd(filename=self._trajectory_path, top=configuration_path)

        return True


class ExtractAverageDensity(AverageTrajectoryProperty):
    """Extracts the average density from a simulation trajectory.
    """

    def __init__(self):

        super().__init__()

        self._system = None

    @Protocol.InputPipe
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


class ExtractAverageDielectric(AverageTrajectoryProperty):
    """Extracts the average dielectric constant from a simulation trajectory.
    """
    def __init__(self):
        super().__init__()

        self._system = None
        self._thermodynamic_state = None

    @Protocol.InputPipe
    def system(self, value):
        pass

    @Protocol.InputPipe
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

        if super(ExtractAverageDielectric, self).execute(directory) is None:
            return False

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

        self._value = dielectric
        self._uncertainty = dielectric_sigma

        logging.info('Extracted dielectrics: ' + directory)

        return True
