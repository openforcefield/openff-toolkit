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

import numpy as np

from os import path
from enum import Enum

from pymbar import timeseries

from pydantic import BaseModel
from typing import Dict, List, Any, Optional

from openeye import oechem, oeomega

from openforcefield.utils import packmol, graph
from openforcefield.utils.exceptions import XmlNodeMissingException
from openforcefield.utils.serialization import quantity_to_json

from openforcefield.typing.engines import smirnoff

from simtk import openmm, unit
from simtk.openmm import app


# =============================================================================================
# Registration Decorators
# =============================================================================================

available_protocols = []


def register_calculation_protocol():
    """A decorator which registers a class as being a
    protocol which may be used in calculation schemas.
    """

    def decorator(cls):
        available_protocols.append(cls)
        return cls

    return decorator


# =============================================================================================
# Protocols
# =============================================================================================

class ProtocolInputReference(BaseModel):
    """Stores a reference to a required input from another protocol.

    Each node represents a protocol to be executed.

    .. warning::

        This class is still heavily under development and is subject to rapid changes.

    .. todo::

        Rename this class to something more obvious - ProtocolDependency?

    Parameters
    ----------
    input_property_name: str
        The name of the property which will take the an output value and use it as input.
    output_protocol_id: str
        The identity of the protocol whose output will be passed to input_property_name.
    output_property_name: str
        The name of the property which will output the required input.

    Attributes
    ----------
    input_property_name: str
        The name of the property which will take the an output value and use it as input.
    output_protocol_id: str
        The identity of the protocol whose output will be passed to input_property_name.
    output_property_name: str
        The name of the property which will output the required input.
    grouped_protocol_id: str, optional
        The name of the protocol which has been grouped to take output from. When set,
        `output_protocol_id` should refer to the name of the ProtocolGroup which contains
        the protocol identified by `grouped_protocol_id`.
    """

    input_property_name: str = None

    output_protocol_id: str = None
    output_property_name: str = None

    grouped_protocol_id: Optional[str] = None

    def __hash__(self):
        """Returns the hash key of this ProtocolInputReference."""
        return hash((self.input_property_name,
                     self.output_protocol_id,
                     self.output_property_name,
                     self.grouped_protocol_id))

    def __eq__(self, other):
        """Returns true if the two inputs are equal."""
        return (self.input_property_name == other.input_property_name and
                self.output_protocol_id == other.output_protocol_id and
                self.output_property_name == other.output_property_name and
                self.grouped_protocol_id == other.grouped_protocol_id)

    def __ne__(self, other):
        """Returns true if the two inputs are not equal."""
        return not (self == other)

    def set_uuid(self, uuid):
        """Appends a uuid to each of the protocol ids

        Notes
        ----------
        Existing uuid's will be overwritten.

        Parameters
        ----------
        uuid : str
            The uuid to append.
        """

        if self.output_protocol_id is not None and self.output_protocol_id != 'global':
            self.output_protocol_id = graph.append_uuid(self.output_protocol_id, uuid)
        if self.grouped_protocol_id is not None and self.grouped_protocol_id != 'global':
            self.grouped_protocol_id = graph.append_uuid(self.grouped_protocol_id, uuid)

    def replace_protocol(self, old_id, new_id):
        """Redirect the input to point at a new protocol.

        The main use of this method is when merging multiple protocols
        into one.

        Parameters
        ----------
        old_id : str
            The id of the protocol to replace.
        new_id : str
            The id of the new protocol to use.
        """
        if self.output_protocol_id is not None:
            self.output_protocol_id = self.output_protocol_id.replace(old_id, new_id)
        if self.grouped_protocol_id is not None:
            self.grouped_protocol_id = self.grouped_protocol_id.replace(old_id, new_id)


class ProtocolSchema(BaseModel):

    id: str = None
    type: str = None

    input_references: List[ProtocolInputReference] = []
    parameters: Dict[str, Any] = {}


class BaseProtocol:
    """The base class for a protocol which would form one
    step of a property calculation.

    A protocol may for example:

        - create the coordinates of a mixed simulation box
        - set up a bound ligand-protein system
        - build the simulation topology
         - perform an energy minimisation

    Protocols may be chained together, thus defining
    a larger property calculation from simple building blocks.

    .. warning::

        This class is still heavily under development and is subject to rapid changes.

    Attributes
    ----------
    id : str, optional
        The unique identity of the protocol
    input_references : list of ProtocolInputReference
        A list of the inputs which this protocol will receive.
    self.required_inputs : list of str
        A list of the inputs that must be passed to this protocol.
    self.provided_outputs : list of str
        A list of the outputs that this protocol produces.
    """

    class ProtocolArgumentDecorator(object):
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

    class InputPipe(ProtocolArgumentDecorator):
        """A custom decorator used to mark properties as a required input to
        the protocol.

        Examples
        ----------
        To mark a property as an input pipe:

        >>> @BaseProtocol.InputPipe
        >>> def substance(self):
        >>>     pass
        """
        def __init__(self, attribute, documentation=None):
            super().__init__(attribute, documentation)

    class OutputPipe(ProtocolArgumentDecorator):
        """A custom decorator used to mark properties as an output of the
        the protocol.

        Examples
        ----------
        To mark a property as an output pipe:

        >>> @BaseProtocol.OutputPipe
        >>> def positions(self):
        >>>     pass
        """
        def __init__(self, attribute, documentation=None):
            super().__init__(attribute, documentation)

    class Parameter(ProtocolArgumentDecorator):
        """A custom decorator used to mark arguments as a settable parameter of
        the protocol.

        Examples
        ----------
        To mark an attribute as a parameter:

        >>> @BaseProtocol.Parameter
        >>> def number_of_steps(self):
        >>>     pass
        """
        def __init__(self, attribute, documentation=None):
            super().__init__(attribute, documentation)

    def __init__(self):

        # A unique identifier for this node.
        self.id = None

        # Defines where to pull the values from.
        self.input_references = []

        # Find the required inputs and outputs.
        self.parameters = self._find_types_with_decorator(BaseProtocol.Parameter)

        self.required_inputs = self._find_types_with_decorator(BaseProtocol.InputPipe)
        self.provided_outputs = self._find_types_with_decorator(BaseProtocol.OutputPipe)

    @property
    def schema(self):
        """ProtocolSchema: Returns a serializable schema for this object."""
        schema = ProtocolSchema()

        schema.id = self.id
        schema.type = str(type(self))

        schema.input_references = self.input_references

        for parameter in self.parameters:

            value = getattr(self, parameter)

            if isinstance(value, unit.Quantity):
                value = quantity_to_json(value)

            schema.parameters[parameter] = value

        return schema

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

    def set_uuid(self, value):
        """Store the uuid of the calculation this protocol belongs to

        Parameters
        ----------
        value : str
            The uuid of the parent calculation.
        """
        if self.id.find(value) >= 0:
            return

        self.id = graph.append_uuid(self.id, value)

        for input_reference in self.input_references:
            input_reference.set_uuid(value)

    def replace_protocol(self, old_id, new_id):
        """Finds each input which came from a given protocol
         and redirects it to instead take input from a new one.

        The main use of this method is when merging multiple protocols
        into one.

        Parameters
        ----------
        old_id : str
            The id of the old input protocol.
        new_id : str
            The id of the new input protocol.
        """
        for input_reference in self.input_references:
            input_reference.replace_protocol(old_id, new_id)

    def set_global_properties(self, global_properties):
        """Set the value of any inputs which takes values
        from the 'global' (i.e property to calculate) scope

        Parameters
        ----------
        global_properties: dict of str to object
            The list of global properties to draw from.
        """
        for input_reference in self.input_references:

            if input_reference.output_protocol_id != 'global':
                continue

            if input_reference.output_property_name not in global_properties:
                raise Exception('Invalid global property: {}'.format(input_reference.output_property_name))

            self.set_input_value(input_reference, global_properties[input_reference.output_property_name])

    def can_merge(self, other):
        """Determines whether this protocol can be merged with another.

        Parameters
        ----------
        other : BaseProtocol
            The protocol to compare against.

        Returns
        ----------
        bool
            True if the two protocols are safe to merge.
        """
        if not isinstance(self, type(other)):
            return False

        for input_reference in self.input_references:

            if input_reference not in other.input_references:
                return False

            self_value = self.get_input_value(input_reference)
            other_value = other.get_input_value(input_reference)

            if self_value != other_value:
                return False

        return True

    def merge(self, other):
        """Merges another BaseProtocol with this one. The id
        of this protocol will remain unchanged.

        It is assumed that can_merge has already returned that
        these protocols are compatible to be merged together.

        Parameters
        ----------
        other: BaseProtocol
            The protocol to merge into this one.
        """

        pass

    def set_input_value(self, input_reference, value):
        """Set the value of one of the protocols inputs.

        Parameters
        ----------
        input_reference: ProtocolInputReference
            The input to set.
        value
            The value to set the input to.
        """
        setattr(self, input_reference.input_property_name, value)

    def get_input_value(self, input_reference):
        """Gets the value that was set on one of this protocols inputs.

        Parameters
        ----------
        input_reference: ProtocolInputReference
            The input to get.

        Returns
        ----------
        object:
            The value of the input
        """
        return getattr(self, input_reference.input_property_name)

    def get_output_value(self, input_reference):
        """Returns the value of one of this protocols outputs.

        Parameters
        ----------
        input_reference: ProtocolInputReference
            An input reference which points to the output to return.

        Returns
        ----------
        object:
            The value of the input
        """
        return getattr(self, input_reference.output_property_name)

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
        inputs = []

        def get_bases(current_base_type):

            bases = [current_base_type]

            for base_type in current_base_type.__bases__:
                bases.extend(get_bases(base_type))

            return bases

        all_bases = get_bases(type(self))

        for base in all_bases:

            inputs.extend([attribute_name for attribute_name in base.__dict__ if
                           type(base.__dict__[attribute_name]) is decorator_type])

        return inputs

    @classmethod
    def from_xml(cls, xml_node, existing_protocols=None):
        """ Creates a protocol from an xml definition.

        Parameters
        ----------
        xml_node : xml.etree.Element
            The element containing the xml to create the protocol from.
        existing_protocols : dict(str, BaseProtocol)
            A list of already created protocols.
        Returns
        ----------
        BaseProtocol
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

                text_split = input_node.text.split(':')

                if len(text_split) != 2:
                    raise Exception('Protocol inputs must be of the form node_id:property_name')

                # Only process inputs which are actually required.
                if input_property not in return_value.required_inputs:
                    continue

                protocol_id = text_split[0]
                property_name = text_split[1]

                protocol_input = ProtocolInputReference(input_property,
                                                        protocol_id,
                                                        property_name)

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


@register_calculation_protocol()
class BuildLiquidCoordinates(BaseProtocol):
    """Create 3D coordinates and bond information for a given Substance

    The coordinates are created using packmol.

    Attributes
    ----------
    _max_molecules : int, optional, default=True
        The maxmimum number of molecules in the system to be created.
    _mass_density : float, simtk.unit.Quantity, or None; optional, default=None
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
        self._max_molecules = 128
        self._mass_density = 1.0 * unit.grams / unit.milliliters

    @BaseProtocol.Parameter
    def max_molecules(self):
        pass

    @BaseProtocol.Parameter
    def mass_density(self):
        pass

    @BaseProtocol.InputPipe
    def substance(self):
        pass

    @BaseProtocol.OutputPipe
    def topology(self):
        pass

    @BaseProtocol.OutputPipe
    def positions(self):
        pass

    @BaseProtocol.OutputPipe
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

        n_copies = np.random.multinomial(self._max_molecules - self._substance.number_of_impurities,
                                         pvals=mole_fractions)

        # Each impurity must have exactly one molecule
        for (index, component) in enumerate(self._substance.components):

            if component.impurity:
                n_copies[index] = 1

        # Create packed box
        topology, positions = packmol.pack_box(molecules, n_copies, mass_density=self._mass_density)

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

    def can_merge(self, protocol):

        return super(BuildLiquidCoordinates, self).can_merge(protocol) and \
               self._max_molecules == protocol.max_molecules and \
               self._mass_density == protocol.mass_density


@register_calculation_protocol()
class BuildSmirnoffTopology(BaseProtocol):
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

    @BaseProtocol.InputPipe
    def force_field(self, value):
        pass

    @BaseProtocol.InputPipe
    def molecules(self, value):
        pass

    @BaseProtocol.InputPipe
    def topology(self, value):
        pass

    @BaseProtocol.OutputPipe
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


@register_calculation_protocol()
class RunEnergyMinimisation(BaseProtocol):
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

    @BaseProtocol.InputPipe
    def topology(self, value):
        pass

    @BaseProtocol.InputPipe
    def positions(self, value):
        pass

    @BaseProtocol.InputPipe
    def system(self, value):
        pass

    @BaseProtocol.OutputPipe
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

    def can_merge(self, protocol):

        # TODO: Properly implement comparison
        return super(RunEnergyMinimisation, self).can_merge(protocol)


@register_calculation_protocol()
class RunOpenMMSimulation(BaseProtocol):
    """Performs a molecular dynamics simulation in a given ensemble using OpenMM

    Attributes
    ----------
    _steps : int
        The number of steps to run the simulation for
    _timestep : float
        The timestep of the integrator.
    _output_frequency : int
        The frequency with which to store simulation data.
    _ensemble : RunOpenMMSimulation.Ensemble
        The ensemble to run the simulation in.
    """

    class Ensemble(Enum):
        """An enum describing the available ensembles.
        """
        NVT = 0
        NPT = 1

    def __init__(self):

        super().__init__()

        self._steps = 1000

        self._thermostat_friction = 1.0 / unit.picoseconds
        self._timestep = 0.001 * unit.picoseconds

        self._output_frequency = 1000

        self._ensemble = self.Ensemble.NPT

        # keep a track of the simulation object in case we need to restart.
        self._simulation_object = None

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

    @BaseProtocol.Parameter
    def steps(self):
        pass

    @BaseProtocol.Parameter
    def thermostat_friction(self):
        pass

    @BaseProtocol.Parameter
    def timestep(self):
        pass

    @BaseProtocol.Parameter
    def output_frequency(self):
        pass

    @BaseProtocol.Parameter
    def ensemble(self):
        pass

    @BaseProtocol.InputPipe
    def thermodynamic_state(self):
        pass

    @BaseProtocol.InputPipe
    def topology(self):
        pass

    @BaseProtocol.InputPipe
    def positions(self):
        pass

    @BaseProtocol.InputPipe
    def system(self):
        pass

    @BaseProtocol.OutputPipe
    def final_positions(self):
        pass

    @BaseProtocol.OutputPipe
    def trajectory(self):
        pass

    @BaseProtocol.OutputPipe
    def statistics(self):
        pass

    def execute(self, directory):

        temperature = self._thermodynamic_state.temperature
        pressure = self._thermodynamic_state.pressure

        if temperature is None:
            logging.error('A temperature must be set to perform a simulation in any ensemble: ' + directory)
            return False
        if self._ensemble is self.Ensemble.NPT and pressure is None:
            logging.error('A pressure must be set to perform an NPT simulation: ' + directory)
            return False

        logging.info('Performing a simulation in the ' + str(self._ensemble) + ' ensemble: ' + directory)

        if self._simulation_object is None:
            self._simulation_object = self._setup_new_simulation(directory, pressure, temperature)

        try:
            self._simulation_object.step(self._steps)
        except Exception:
            logging.warning('Failed to run in ' + directory)
            return False

        positions = self._simulation_object.context.getState(getPositions=True).getPositions()

        self._final_positions = positions

        logging.info('Simulation performed in the ' + str(self._ensemble) + ' ensemble: ' + directory)

        configuration_path = path.join(directory, 'output.pdb')

        with open(configuration_path, 'w+') as configuration_file:

            app.PDBFile.writeFile(self._topology,
                                  positions, configuration_file)

        return True

    def _setup_new_simulation(self, directory, pressure, temperature):
        # For now set some 'best guess' thermostat parameters.
        integrator = openmm.LangevinIntegrator(temperature,
                                               self._thermostat_friction,
                                               self._timestep)

        system = self._system

        if self._ensemble is self.Ensemble.NPT:
            barostat = openmm.MonteCarloBarostat(pressure, temperature)

            # inputs are READONLY! Never directly alter an input
            system = copy.deepcopy(system)
            system.addForce(barostat)

        simulation = app.Simulation(self._topology, system, integrator)
        simulation.context.setPositions(self._positions)

        simulation.context.setVelocitiesToTemperature(temperature)

        trajectory_path = path.join(directory, 'trajectory.dcd')
        statistics_path = path.join(directory, 'statistics.dat')

        self._trajectory = trajectory_path
        self._statistics = statistics_path

        configuration_path = path.join(directory, 'input.pdb')

        with open(configuration_path, 'w+') as configuration_file:
            app.PDBFile.writeFile(self._topology,
                                  self._positions, configuration_file)

        simulation.reporters.append(app.DCDReporter(trajectory_path, self._output_frequency))

        simulation.reporters.append(app.StateDataReporter(statistics_path, self._output_frequency, step=True,
                                                          potentialEnergy=True, temperature=True, volume=True))

        return simulation

    def __getstate__(self):
        state = self.__dict__.copy()

        if self._simulation_object is not None:
            del state['_simulation_object']

        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

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

    def can_merge(self, protocol):

        return super(RunOpenMMSimulation, self).can_merge(protocol) and \
               self._ensemble == protocol.ensemble


@register_calculation_protocol()
class AveragePropertyProtocol(BaseProtocol):
    """Calculates the average of a property and its uncertainty.
    """

    def __init__(self):

        super().__init__()

        self._value = None
        self._uncertainty = None

    @BaseProtocol.OutputPipe
    def value(self):
        pass

    @BaseProtocol.OutputPipe
    def uncertainty(self):
        pass

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


@register_calculation_protocol()
class AverageTrajectoryProperty(AveragePropertyProtocol):
    """Calculates the average of a property from a simulation trajectory.
    """

    def __init__(self):

        super().__init__()

        self._topology = None
        self._positions = None
        self._trajectory_path = None

        self.trajectory = None

    @BaseProtocol.InputPipe
    def topology(self, value):
        pass

    @BaseProtocol.InputPipe
    def positions(self, value):
        pass

    @BaseProtocol.InputPipe
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
