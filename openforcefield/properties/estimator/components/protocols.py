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
import pickle

import numpy as np

import mdtraj

from os import path
from enum import Enum

from pymbar import timeseries

from pydantic import BaseModel
from typing import Dict, Any

from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils import packmol, graph, utils
from openforcefield.utils.serialization import serialize_quantity, deserialize_quantity, TypedBaseModel

from openforcefield.typing.engines import smirnoff

from simtk import openmm, unit
from simtk.openmm import app

# =============================================================================================
# Registration Decorators
# =============================================================================================

available_protocols = {}


def register_calculation_protocol():
    """A decorator which registers a class as being a
    protocol which may be used in calculation schemas.
    """

    def decorator(cls):

        if cls.__name__ in available_protocols:
            raise ValueError('The {} protocol is already registered.'.format(cls.__name__))

        available_protocols[cls.__name__] = cls
        return cls

    return decorator


# =============================================================================================
# Protocols
# =============================================================================================

class PropertyCalculatorException(BaseModel):
    """A json serializable object wrapper containing information about
    a failed property calculation.

    TODO: Flesh out more fully.
    """
    directory: str
    message: str


class ProtocolPath:
    """
    TODO: Document.
    """

    path_separator = '/'  # The character which separates protocol ids.
    property_separator = '.'  # The character which separates the property name from the path.

    @property
    def property_name(self):
        """str: The property name pointed to by the path."""
        property_name, protocol_ids = ProtocolPath.to_components(self._full_path)
        return property_name

    @property
    def start_protocol(self):
        """str: The leading protocol id of the path."""
        property_name, protocol_ids = ProtocolPath.to_components(self._full_path)
        return None if len(protocol_ids) == 0 else protocol_ids[0]

    @property
    def last_protocol(self):
        """str: The leading protocol id of the path."""
        property_name, protocol_ids = ProtocolPath.to_components(self._full_path)
        return None if len(protocol_ids) == 0 else protocol_ids[len(protocol_ids) - 1]

    @property
    def full_path(self):
        """str: The full path referenced by this object."""
        return self._full_path

    @property
    def is_global(self):
        return self.start_protocol == 'global'

    def __init__(self, property_name, *protocol_ids):
        """Constructs a new ProtocolPath object.

        Parameters
        ----------
        property_name: str
            The property name referenced by the path.
        protocol_ids: tuple of str
            An args list of protocol ids in the order in which they will appear in the path.
        """

        self._full_path = ''
        self._from_components(property_name, *protocol_ids)

    def _from_components(self, property_name, *protocol_ids):
        """Sets this components path from individual components.

        Parameters
        ----------
        property_name: str
            The property name referenced by the path.
        protocol_ids: list of str
            A list of protocol ids in the order in which they will appear in the path.
        """

        assert property_name is not None and isinstance(property_name, str)

        assert property_name.find(ProtocolPath.property_separator) < 0 and \
               property_name.find(ProtocolPath.path_separator) < 0

        for protocol_id in protocol_ids:

            assert protocol_id is not None and isinstance(protocol_id, str)

            assert protocol_id.find(ProtocolPath.property_separator) < 0 and \
                   protocol_id.find(ProtocolPath.path_separator) < 0

        protocol_path = ProtocolPath.path_separator.join(protocol_ids)

        if len(protocol_ids) == 0:
            protocol_path = ''

        self._full_path = '{}{}{}'.format(protocol_path,
                                          ProtocolPath.property_separator,
                                          property_name)

    @classmethod
    def from_string(cls, existing_path_string: str):

        existing_path_string = existing_path_string.lstrip().rstrip()
        property_name_index = existing_path_string.find(ProtocolPath.property_separator)

        if property_name_index < 0:

            raise ValueError('A protocol path must contain a {} followed by the '
                             'property name this path represents'.format(ProtocolPath.property_separator))

        if existing_path_string.find(ProtocolPath.property_separator, property_name_index + 1) >= 0:

            raise ValueError('A protocol path must contain at most one '
                             'property separator ({})'.format(ProtocolPath.property_separator))

        property_name, protocol_ids = ProtocolPath.to_components(existing_path_string)

        for protocol_id in protocol_ids:

            if protocol_id is not None and len(protocol_id) > 0:
                continue

            raise ValueError('An invalid protocol id (either None or empty) was found.')

        return ProtocolPath(property_name, *protocol_ids)

    @staticmethod
    def to_components(path_string):
        """Splits a protocol path string into the property
        name, and the individual protocol ids.

        Parameters
        ----------
        path_string: str
            The protocol path to split.

        Returns
        -------
        str, list of str
            A tuple of the property name, and a list of the protocol ids in the path.
        """
        property_name_index = path_string.find(ProtocolPath.property_separator)
        property_name = path_string[property_name_index + 1:]

        protocol_id_path = path_string[:property_name_index]
        protocol_ids = protocol_id_path.split(ProtocolPath.path_separator)

        if len(protocol_id_path) == 0:
            protocol_ids = []

        return property_name, protocol_ids

    def prepend_protocol_id(self, id_to_prepend):
        """Prepend a new protocol id onto the front of the path.

        Parameters
        ----------
        id_to_prepend: str
            The protocol id to prepend to the path
        """
        property_name, protocol_ids = ProtocolPath.to_components(self._full_path)

        if len(protocol_ids) == 0 or (len(protocol_ids) > 0 and protocol_ids[0] != id_to_prepend):
            protocol_ids.insert(0, id_to_prepend)

        self._from_components(property_name, *protocol_ids)

    def pop_next_in_path(self):
        """Pops and then returns the leading protocol id from the path.

        Returns
        -------
        str:
            The previously leading protocol id.
        """
        property_name, protocol_ids = ProtocolPath.to_components(self._full_path)

        if len(protocol_ids) == 0:
            return None

        next_in_path = protocol_ids.pop(0)
        self._from_components(property_name, *protocol_ids)

        return next_in_path

    def append_uuid(self, uuid):
        """Appends a uuid to each of the protocol id's in the path

        Parameters
        ----------
        uuid: str
            The uuid to append.
        """

        if self.is_global:
            # Don't append uuids to global paths.
            return

        property_name, protocol_ids = ProtocolPath.to_components(self._full_path)
        appended_ids = []

        for protocol_id in protocol_ids:

            if protocol_id is None:
                continue

            appended_id = graph.append_uuid(protocol_id, uuid)
            appended_ids.append(appended_id)

        self._from_components(property_name, *appended_ids)

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
        self._full_path = self._full_path.replace(old_id, new_id)

    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):

        if isinstance(v, str):
            return ProtocolPath.from_string(v)

        return v

    def __str__(self):
        return self._full_path

    def __repr__(self):
        return '<ProtocolPath full_path={}>'.format(self._full_path)

    def __hash__(self):
        """Returns the hash key of this ProtocolPath."""
        return hash(self._full_path)

    def __eq__(self, other):
        """Returns true if the two inputs are equal."""
        return self._full_path == other.full_path

    def __ne__(self, other):
        """Returns true if the two inputs are not equal."""
        return not (self == other)

    def __getstate__(self):
        return {'full_path': self._full_path}

    def __setstate__(self, state):
        self._full_path = state['full_path']


class ProtocolSchema(TypedBaseModel):
    """A json serializable representation which stores the
    user definable parameters of a protocol.

    TODO: Merge inputs and parameters.
    """
    id: str = None
    type: str = None

    inputs: Dict[str, ProtocolPath] = {}
    parameters: Dict[str, Any] = {}

    class Config:
        arbitrary_types_allowed = True

        json_encoders = {
            ProtocolPath: lambda v: v.full_path
        }


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
    self.required_inputs : list of ProtocolPath
        A list of the inputs that must be passed to this protocol.
    self.provided_outputs : list of ProtocolPath
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

    @property
    def id(self):
        return self._id

    @property
    def schema(self):
        """ProtocolSchema: Returns a serializable schema for this object."""
        return self._get_schema()

    @schema.setter
    def schema(self, schema_value):
        self._set_schema(schema_value)

    @property
    def dependencies(self):

        return_dependencies = []

        for input_path in self.required_inputs:

            input_value = self.get_value(input_path)

            if not isinstance(input_value, ProtocolPath):
                continue

            if input_value not in return_dependencies:
                return_dependencies.append(input_value)

        return return_dependencies

    @Parameter
    def allow_merging(self):
        pass

    def __init__(self, protocol_id):

        # A unique identifier for this node.
        self._id = protocol_id

        # Defines whether a protocol is allowed to try and merge with other identical ones.
        self._allow_merging = True

        # Find the required inputs and outputs.
        self.parameters = utils.find_types_with_decorator(type(self), BaseProtocol.Parameter)

        self.provided_outputs = []
        self.required_inputs = []

        output_attributes = utils.find_types_with_decorator(type(self), BaseProtocol.OutputPipe)
        input_attributes = utils.find_types_with_decorator(type(self), BaseProtocol.InputPipe)

        for output_attribute in output_attributes:
            self.provided_outputs.append(ProtocolPath(output_attribute))

        for input_attribute in input_attributes:
            self.required_inputs.append(ProtocolPath(input_attribute))

        # The directory in which to execute the protocol.
        self.directory = None

    def execute(self, directory):
        """ Execute the protocol.

        Protocols may be chained together by passing the output
        of previous protocols as input to the current one.

        Parameters
        ----------
        directory: str
            The directory to store output data in.

        Returns
        ----------
        Dict[str, Any]
            The output of the execution.
        """

        return self._get_output_dictionary()

    def _get_schema(self):

        schema = ProtocolSchema()

        schema.id = self.id
        schema.type = type(self).__name__

        for input_path in self.required_inputs:

            if input_path.start_protocol is None or (input_path.start_protocol == self.id and
                                                     input_path.start_protocol == input_path.last_protocol):
                schema.inputs[input_path.full_path] = self.get_value(input_path)

        for parameter in self.parameters:

            value = getattr(self, parameter)

            if isinstance(value, unit.Quantity):
                value = serialize_quantity(value)

            schema.parameters[parameter] = value

        return schema

    def _set_schema(self, schema_value):
        """Sets this protocols properties (i.e id and parameters)
        from a ProtocolSchema
        """
        self._id = schema_value.id

        if type(self).__name__ != schema_value.type:
            # Make sure this object is the correct type.
            raise ValueError('Cannot convert a {} protocol to a {}.'
                             .format(str(type(self)), schema_value.type))

        for input_full_path in schema_value.inputs:
            input_path = ProtocolPath.from_string(input_full_path)
            self.set_value(input_path, schema_value.inputs[input_full_path])

        for parameter in schema_value.parameters:

            value = schema_value.parameters[parameter]

            if isinstance(value, dict) and 'unit' in value and 'unitless_value' in value:
                value = deserialize_quantity(value)

            setattr(self, parameter, value)

    def _get_output_dictionary(self):
        """Builds a dictionary of the output property names and their values.

        Returns
        -------
        Dict[str, Any]
            A dictionary whose keys are the output property names, and the
            values their associated values.
        """

        return_dictionary = {}

        for output_path in self.provided_outputs:
            return_dictionary[output_path.full_path] = self.get_value(output_path)

        return return_dictionary

    def set_uuid(self, value):
        """Store the uuid of the calculation this protocol belongs to

        Parameters
        ----------
        value : str
            The uuid of the parent calculation.
        """
        if self.id.find(value) >= 0:
            return

        self._id = graph.append_uuid(self.id, value)

        for input_path in self.required_inputs:

            input_path.append_uuid(value)
            input_value = self.get_value(input_path)

            if isinstance(input_value, ProtocolPath):
                input_value.append_uuid(value)

        for output_path in self.provided_outputs:
            output_path.append_uuid(value)

    def replace_protocol(self, old_id, new_id):
        """Finds each input which came from a given protocol
         and redirects it to instead take input from a new one.

        Notes
        -----
        This method is mainly intended to be used only when merging
        multiple protocols into one.

        Parameters
        ----------
        old_id : str
            The id of the old input protocol.
        new_id : str
            The id of the new input protocol.
        """

        for input_path in self.required_inputs:

            input_path.replace_protocol(old_id, new_id)
            input_value = self.get_value(input_path)

            if isinstance(input_value, ProtocolPath):
                input_value.replace_protocol(old_id, new_id)

        for output_path in self.provided_outputs:
            output_path.replace_protocol(old_id, new_id)

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
        if not self.allow_merging:
            return False

        if not isinstance(self, type(other)):
            return False

        for input_path in self.required_inputs:

            # if input_references not in other.input_references:
            #     return False
            if input_path not in other.required_inputs:
                return False

            self_value = self.get_value(input_path)
            other_value = other.get_value(input_path)

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

    def get_value(self, reference_path):
        """Returns the value of one of this protocols parameters / inputs.

        Parameters
        ----------
        reference_path: ProtocolPath
            The path pointing to the value to return.

        Returns
        ----------
        object:
            The value of the input
        """

        if (reference_path.start_protocol is not None and
            reference_path.start_protocol != self.id):

            raise ValueError('The reference path does not target this protocol.')

        if not hasattr(self, reference_path.property_name):

            raise ValueError('This protocol does not have contain a {} '
                             'property.'.format(reference_path.property_name))

        return getattr(self, reference_path.property_name)

    def set_value(self, reference_path, value):
        """Sets the value of one of this protocols parameters / inputs.

        Parameters
        ----------
        reference_path: ProtocolPath
            The path pointing to the value to return.
        value: Any
            The value to set.
        """

        if (reference_path.start_protocol is not None and
            reference_path.start_protocol != self.id):

            raise ValueError('The reference path does not target this protocol.')

        if not hasattr(self, reference_path.property_name):

            raise ValueError('This protocol does not have contain a {} '
                             'property.'.format(reference_path.property_name))

        if reference_path in self.provided_outputs:
            raise ValueError('Output values cannot be set by this method.')

        setattr(self, reference_path.property_name, value)


@register_calculation_protocol()
class BuildCoordinatesPackmol(BaseProtocol):
    """Create 3D coordinates and bond information for a given Substance

    The coordinates are created using packmol.

    Attributes
    ----------
    _max_molecules : int, optional, default=True
        The maximum number of molecules in the system to be created.
    _mass_density : float, simtk.unit.Quantity, or None; optional, default=None
        If provided, will aid in the selecting an initial box size.
    """

    _cached_molecules = {}

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        # inputs
        self._substance = None

        # outputs
        self._coordinate_file = None
        self._positions = None
        self._molecules = None

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
    def coordinate_file(self):
        pass

    @BaseProtocol.OutputPipe
    def molecules(self):
        pass

    def _create_molecule(self, smiles):
        """
        Create molecule from a smiles pattern.

        Todo
        ----
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

        from openeye import oechem, oeomega

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

            return PropertyCalculatorException(directory=directory,
                                               message='The substance input is non-optional')

        molecules = []

        for component in self._substance.components:

            molecule = self._create_molecule(component.smiles)

            if molecule is None:

                return PropertyCalculatorException(directory=directory,
                                                   message='{} could not be converted to a Molecule'.format(component))

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

            return PropertyCalculatorException(directory=directory,
                                               message='Packmol failed to complete.')

        self._molecules = molecules

        self._coordinate_file = path.join(directory, 'output.pdb')

        with open(self._coordinate_file, 'w+') as minimised_file:
            app.PDBFile.writeFile(topology, positions, minimised_file)

        logging.info('Coordinates generated: ' + str(self._substance))

        return self._get_output_dictionary()

    def can_merge(self, protocol):

        return super(BuildCoordinatesPackmol, self).can_merge(protocol) and \
               self._max_molecules == protocol.max_molecules and \
               self._mass_density == protocol.mass_density


@register_calculation_protocol()
class BuildSmirnoffTopology(BaseProtocol):
    """Parametrise a set of molecules with a given smirnoff force field.
    """
    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        # inputs
        self._force_field_path = None
        self._coordinate_file = None
        self._molecules = None
        # outputs
        self._system = None

    @BaseProtocol.InputPipe
    def force_field_path(self, value):
        pass

    @BaseProtocol.InputPipe
    def molecules(self, value):
        pass

    @BaseProtocol.InputPipe
    def coordinate_file(self, value):
        pass

    @BaseProtocol.OutputPipe
    def system(self):
        pass

    def execute(self, directory):

        logging.info('Generating topology: ' + directory)

        pdb_file = app.PDBFile(self._coordinate_file)

        parameter_set = ForceField([])

        with open(self._force_field_path, 'rb') as file:
            parameter_set.__setstate__(pickle.load(file))

        system = parameter_set.createSystem(pdb_file.topology,
                                            self._molecules,
                                            nonbondedMethod=smirnoff.PME,
                                            chargeMethod='OECharges_AM1BCCSym')

        if system is None:

            return PropertyCalculatorException(directory=directory,
                                               message='Failed to create a system from the'
                                                       'provided topology and molecules')

        self._system = system

        logging.info('Topology generated: ' + directory)

        return self._get_output_dictionary()


@register_calculation_protocol()
class RunEnergyMinimisation(BaseProtocol):
    """Minimises the energy of a passed in system.
    """

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        # inputs
        self._input_coordinate_file = None
        self._system = None
        # outputs
        self._output_coordinate_file = None

        # TODO: Add arguments for max iter + tolerance
        pass

    @BaseProtocol.InputPipe
    def input_coordinate_file(self, value):
        pass

    @BaseProtocol.InputPipe
    def system(self, value):
        pass

    @BaseProtocol.OutputPipe
    def output_coordinate_file(self):
        return self._final_positions

    def execute(self, directory):

        logging.info('Minimising energy: ' + directory)

        integrator = openmm.VerletIntegrator(0.002 * unit.picoseconds)

        input_pdb_file = app.PDBFile(self._input_coordinate_file)

        simulation = app.Simulation(input_pdb_file.topology,
                                    self._system, integrator)

        simulation.context.setPositions(input_pdb_file.positions)

        simulation.minimizeEnergy()

        positions = simulation.context.getState(getPositions=True).getPositions()

        self._output_coordinate_file = path.join(directory, 'minimised.pdb')

        with open(self._output_coordinate_file, 'w+') as minimised_file:
            app.PDBFile.writeFile(simulation.topology, positions, minimised_file)

        logging.info('Energy minimised: ' + directory)

        return self._get_output_dictionary()

    def can_merge(self, protocol):
        return self._get_output_dictionary()


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

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        self._steps = 1000

        self._thermostat_friction = 1.0 / unit.picoseconds
        self._timestep = 0.001 * unit.picoseconds

        self._output_frequency = 1000

        self._ensemble = self.Ensemble.NPT

        # keep a track of the simulation object in case we need to restart.
        self._simulation_object = None

        # inputs
        self._input_coordinate_file = None
        self._thermodynamic_state = None
        self._system = None

        # outputs
        self._output_coordinate_file = None
        self._trajectory = None
        self._statistics = None

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
    def input_coordinate_file(self):
        pass

    @BaseProtocol.InputPipe
    def system(self):
        pass

    @BaseProtocol.OutputPipe
    def output_coordinate_file(self):
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

            return PropertyCalculatorException(directory=directory,
                                               message='A temperature must be set to perform '
                                                       'a simulation in any ensemble')

        if self._ensemble is self.Ensemble.NPT and pressure is None:

            return PropertyCalculatorException(directory=directory,
                                               message='A pressure must be set to perform an NPT simulation')

        logging.info('Performing a simulation in the ' + str(self._ensemble) + ' ensemble: ' + directory)

        if self._simulation_object is None:
            self._simulation_object = self._setup_new_simulation(directory, temperature, pressure)

        try:
            self._simulation_object.step(self._steps)
        except Exception as e:

            return PropertyCalculatorException(directory=directory,
                                               message='Simulation failed: {}'.format(e))

        positions = self._simulation_object.context.getState(getPositions=True).getPositions()

        input_pdb_file = app.PDBFile(self._input_coordinate_file)
        self._output_coordinate_file = path.join(directory, 'output.pdb')

        logging.info('Simulation performed in the ' + str(self._ensemble) + ' ensemble: ' + directory)

        with open(self._output_coordinate_file, 'w+') as configuration_file:

            app.PDBFile.writeFile(input_pdb_file.topology,
                                  positions, configuration_file)

        return self._get_output_dictionary()

    def _setup_new_simulation(self, directory, temperature, pressure):
        """Creates a new OpenMM simulation object.

        Parameters
        ----------
        directory: str
            The directory in which the object will produce output files.
        temperature: unit.Quantiy
            The temperature at which to run the simulation
        pressure: unit.Quantiy
            The pressure at which to run the simulation
        """

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

        input_pdb_file = app.PDBFile(self._input_coordinate_file)

        simulation = app.Simulation(input_pdb_file.topology, system, integrator)
        simulation.context.setPositions(input_pdb_file.positions)

        simulation.context.setVelocitiesToTemperature(temperature)

        trajectory_path = path.join(directory, 'trajectory.dcd')
        statistics_path = path.join(directory, 'statistics.dat')

        self._trajectory = trajectory_path
        self._statistics = statistics_path

        configuration_path = path.join(directory, 'input.pdb')

        with open(configuration_path, 'w+') as configuration_file:

            app.PDBFile.writeFile(input_pdb_file.topology,
                                  input_pdb_file.positions, configuration_file)

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

    def merge(self, other):

        super(RunOpenMMSimulation, self).merge(other)

        self._steps = max(self._steps, other.steps)
        self._timestep = min(self._timestep, other.timestep)

        self._output_frequency = max(self._output_frequency, other.output_frequency)

    def can_merge(self, protocol):

        return super(RunOpenMMSimulation, self).can_merge(protocol) and \
               self._ensemble == protocol.ensemble


@register_calculation_protocol()
class AveragePropertyProtocol(BaseProtocol):
    """Calculates the average of a property and its uncertainty.
    """

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        self._value = None
        self._uncertainty = None

    @BaseProtocol.OutputPipe
    def value(self):
        pass

    @BaseProtocol.OutputPipe
    def uncertainty(self):
        pass

    def execute(self, directory):
        return self._get_output_dictionary()

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

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        self._input_coordinate_file = None
        self._trajectory_path = None

        self.trajectory = None

    @BaseProtocol.InputPipe
    def input_coordinate_file(self):
        pass

    @BaseProtocol.InputPipe
    def trajectory_path(self):
        pass

    def execute(self, directory):

        if self._trajectory_path is None:

            return PropertyCalculatorException(directory=directory,
                                               message='The AverageTrajectoryProperty protocol '
                                                       'requires a previously calculated trajectory')

        self.trajectory = mdtraj.load_dcd(filename=self._trajectory_path, top=self._input_coordinate_file)

        return self._get_output_dictionary()
