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
import math
import pickle
from enum import Enum
from os import path
from typing import Dict, Any

import mdtraj
import numpy as np
from pydantic import BaseModel
from simtk import openmm, unit
from simtk.openmm import app

from openforcefield.typing.engines import smirnoff
from openforcefield.utils import packmol, graph, utils, statistics
from openforcefield.utils.serialization import serialize_quantity, deserialize_quantity, TypedBaseModel

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
    step of a larger property calculation workflow.

    A protocol may for example:

        * create the coordinates of a mixed simulation box
        * set up a bound ligand-protein system
        * build the simulation topology
        * perform an energy minimisation

    An individual protocol may require a set of inputs, which may either be
    set as constants

    >>> npt_equilibration = RunOpenMMSimulation('npt_equilibration')
    >>> npt_equilibration.ensemble = RunOpenMMSimulation.Ensemble.NPT

    or from the output of another protocol, pointed to by a ProtocolPath

    >>> npt_production = RunOpenMMSimulation('npt_production')
    >>> # Use the coordinate file output by the npt_equilibration protocol
    >>> # as the input to the npt_production protocol
    >>> npt_production.input_coordinate_file = ProtocolPath('output_coordinate_file',
    >>>                                                     npt_equilibration.id)

    In this way protocols may be chained together, thus defining a larger property
    calculation workflow from simple, reusable building blocks.

    .. warning:: This class is still heavily under development and is subject to
                 rapid changes.
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
    """Creates a set of 3D coordinates a system defined by a
    ``openforcefield.properties.Substance``.

    Notes
    -----
    The coordinates are created using packmol.
    """

    _cached_molecules = {}

    @BaseProtocol.Parameter
    def max_molecules(self):
        """int: The maximum number of molecules to be added to the system."""
        pass

    @BaseProtocol.Parameter
    def mass_density(self):
        """unit.Quantity: The target density of the created system."""
        pass

    @BaseProtocol.InputPipe
    def substance(self):
        """Substance: The composition of the system to build."""
        pass

    @BaseProtocol.OutputPipe
    def coordinate_file(self):
        """str: A path to the created PDB coordinate file."""
        pass

    @BaseProtocol.OutputPipe
    def molecules(self):
        """list of OEMol: A list of the molecule types in the built system."""
        pass

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

    def _create_molecule(self, smiles):
        """
        Create an ``OEMol`` molecule from a smiles pattern.

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

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        # inputs
        self._force_field_path = None
        self._coordinate_file = None
        self._molecules = None
        # outputs
        self._system = None

    def execute(self, directory):

        logging.info('Generating topology: ' + directory)

        pdb_file = app.PDBFile(self._coordinate_file)

        parameter_set = None

        with open(self._force_field_path, 'rb') as file:
            parameter_set = pickle.load(file)

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

    @BaseProtocol.InputPipe
    def input_coordinate_file(self, value):
        pass

    @BaseProtocol.InputPipe
    def system(self, value):
        pass

    @BaseProtocol.OutputPipe
    def output_coordinate_file(self):
        return self._final_positions

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        # inputs
        self._input_coordinate_file = None
        self._system = None
        # outputs
        self._output_coordinate_file = None

        # TODO: Add arguments for max iter + tolerance
        pass

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

    def execute(self, directory):

        temperature = self._thermodynamic_state.temperature
        pressure = self._thermodynamic_state.pressure

        if temperature is None:

            return PropertyCalculatorException(directory=directory,
                                               message='A temperature must be set to perform '
                                                       'a simulation in any ensemble')

        if self.Ensemble(self._ensemble) == self.Ensemble.NPT and pressure is None:

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

        if self.Ensemble(self._ensemble) == self.Ensemble.NPT:

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
    """An abstract base class for protocols which will calculate the
    average of a property and its uncertainty via bootstrapping.
    """

    @BaseProtocol.Parameter
    def bootstrap_iterations(self):
        pass

    @BaseProtocol.Parameter
    def bootstrap_sample_size(self):
        pass

    @BaseProtocol.OutputPipe
    def value(self):
        pass

    @BaseProtocol.OutputPipe
    def uncertainty(self):
        pass

    @BaseProtocol.OutputPipe
    def equilibration_index(self):
        pass

    @BaseProtocol.OutputPipe
    def statistical_inefficiency(self):
        pass

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        self._bootstrap_iterations = 100
        self._bootstrap_sample_size = 1.0

        self._value = None
        self._uncertainty = None

        self._equilibration_index = None
        self._statistical_inefficiency = None

    def _bootstrap_function(self, sample_data):
        """The function to perform on the data set being sampled by
        bootstrapping.

        Parameters
        ----------
        sample_data: np.ndarray, shape=(num_frames, num_dimensions), dtype=float
            A sample of the full data set.

        Returns
        -------
        float
            The result of evaluating the data.
        """
        return sample_data.mean()

    def _perform_bootstrapping(self, data_set):
        """Performs bootstrapping on a data set to calculate the
        average value, and the standard error in the average,
        bootstrapping.

        Parameters
        ----------
        data_set: np.ndarray, shape=(num_frames, num_dimensions), dtype=float
            The data set to perform bootstrapping on.

        Returns
        -------
        float
            The average of the data.
        float
            The uncertainty in the average.
        """

        if data_set is None:
            raise ValueError('There is no data to bootstrap in protocol {}'.format(self.id))

        # Make a copy of the data so we don't accidentally destroy anything.
        data_to_bootstrap = np.array(data_set)

        data_size = len(data_to_bootstrap)

        # Choose the sample size as a percentage of the full data set.
        sample_size = min(math.floor(data_size * self._bootstrap_sample_size), data_size)

        average_values = np.zeros(self._bootstrap_iterations)

        for bootstrap_iteration in range(self._bootstrap_iterations):

            sample_indices = np.random.choice(data_size, sample_size)
            sample_data = data_to_bootstrap[sample_indices]

            average_values[bootstrap_iteration] = self._bootstrap_function(sample_data)

        average_value = self._bootstrap_function(data_to_bootstrap)
        uncertainty = average_values.std() * len(average_values) ** -0.5

        if isinstance(average_value, np.float32) or isinstance(average_value, np.float64):
            average_value = average_value.item()

        if isinstance(uncertainty, np.float32) or isinstance(uncertainty, np.float64):
            uncertainty = uncertainty.item()

        return average_value, uncertainty

    def execute(self, directory):
        return self._get_output_dictionary()


@register_calculation_protocol()
class AverageTrajectoryProperty(AveragePropertyProtocol):
    """An abstract base class for protocols which will calculate the
    average of a property from a simulation trajectory.
    """

    @BaseProtocol.InputPipe
    def input_coordinate_file(self):
        pass

    @BaseProtocol.InputPipe
    def trajectory_path(self):
        pass

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        self._input_coordinate_file = None
        self._trajectory_path = None

        self.trajectory = None

    def execute(self, directory):

        if self._trajectory_path is None:

            return PropertyCalculatorException(directory=directory,
                                               message='The AverageTrajectoryProperty protocol '
                                                       'requires a previously calculated trajectory')

        self.trajectory = mdtraj.load_dcd(filename=self._trajectory_path, top=self._input_coordinate_file)

        return self._get_output_dictionary()


@register_calculation_protocol()
class ExtractUncorrelatedData(BaseProtocol):
    """An abstract base class for protocols which will subsample
    a data set, yielding only equilibrated, uncorrelated data.
    """

    @BaseProtocol.InputPipe
    def equilibration_index(self):
        pass

    @BaseProtocol.InputPipe
    def statistical_inefficiency(self):
        pass

    def __init__(self, protocol_id):
        super().__init__(protocol_id)

        self._equilibration_index = None
        self._statistical_inefficiency = None

    def execute(self, directory):
        raise NotImplementedError


@register_calculation_protocol()
class ExtractUncorrelatedTrajectoryData(ExtractUncorrelatedData):
    """A protocol which will subsample frames from a trajectory, yielding only uncorrelated 
    frames as determined from a provided statistical inefficiency and equilibration time.
    """

    @BaseProtocol.InputPipe
    def input_coordinate_file(self):
        pass

    @BaseProtocol.InputPipe
    def input_trajectory_path(self):
        pass

    @BaseProtocol.OutputPipe
    def output_trajectory_path(self):
        pass

    def __init__(self, protocol_id):

        super().__init__(protocol_id)

        self._input_coordinate_file = None
        self._input_trajectory_path = None

        self._output_trajectory_path = None

    def execute(self, directory):

        if self._input_trajectory_path is None:

            return PropertyCalculatorException(directory=directory,
                                               message='The ExtractUncorrelatedTrajectoryData protocol '
                                                       'requires a previously calculated trajectory')

        trajectory = mdtraj.load_dcd(filename=self._input_trajectory_path, top=self._input_coordinate_file)

        uncorrelated_indices = statistics.get_uncorrelated_indices(trajectory.n_frames, self._statistical_inefficiency)
        uncorrelated_trajectory = trajectory[uncorrelated_indices]

        self._output_trajectory_path = path.join(directory, 'uncorrelated_trajectory.dcd')
        uncorrelated_trajectory.save_dcd(self._output_trajectory_path)

        return self._get_output_dictionary()