#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================
"""
Parameter assignment tools for the SMIRNOFF (SMIRKS Native Open Force Field) format.

.. codeauthor:: John D. Chodera <john.chodera@choderalab.org>
.. codeauthor:: David L. Mobley <dmobley@mobleylab.org>
.. codeauthor:: Peter K. Eastman <peastman@stanford.edu>

.. todo ::

   * Speed up overall import time by putting non-global imports only where they are needed

"""

__all__ = [
    'MAX_SUPPORTED_VERSION',
    'ParameterHandlerRegistrationError',
    'SMIRNOFFVersionError',
    'SMIRNOFFAromaticityError',
    'ParseError',
    'ForceField',
]


#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import copy
import logging
import os

from collections import OrderedDict


from simtk import openmm, unit

from openforcefield.utils import all_subclasses, MessageException
from openforcefield.topology.molecule import DEFAULT_AROMATICITY_MODEL
from openforcefield.typing.engines.smirnoff.parameters import ParameterList, ParameterHandler
from openforcefield.typing.engines.smirnoff.io import ParameterIOHandler


#=============================================================================================
# CONFIGURE LOGGER
#=============================================================================================

logger = logging.getLogger(__name__)

#=============================================================================================
# PRIVATE METHODS
#=============================================================================================

# Directory paths used by ForceField to discover offxml files.
_installed_offxml_dir_paths = []

def _get_installed_offxml_dir_paths():
    """Return the list of directory paths where to search for offxml files.

    This function load the information by calling all the entry points
    registered in the "openff.forcefielddirs" group. Each entry point
    (i.e. a function) should return a list of paths to directories
    containing the offxml files.

    Returns
    -------
    installed_offxml_dir_paths : List[str]
        All the installed directory paths where ``ForceField`` will
        look for offxml files.

    """
    global _installed_offxml_dir_paths
    if len(_installed_offxml_dir_paths) == 0:
        from pkg_resources import iter_entry_points
        # Find all registered entry points that should return a list of
        # paths to directories where to search for offxml files.
        for entry_point in iter_entry_points(group='openforcefield.smirnoff_forcefield_directory'):
            _installed_offxml_dir_paths.extend(entry_point.load()())
    return _installed_offxml_dir_paths


# TODO: Instead of having a global version number, alow each Force to have a separate version number
MAX_SUPPORTED_VERSION = '1.0'  # maximum version of the SMIRNOFF spec supported by this SMIRNOFF forcefield


class ParameterHandlerRegistrationError(MessageException):
    """
    Exception for errors in ParameterHandler registration
    """
    pass

class SMIRNOFFVersionError(MessageException):
    """
    Exception thrown when an incompatible SMIRNOFF version data structure in attempted to be read.
    """
    pass

class SMIRNOFFAromaticityError(MessageException):
    """
    Exception thrown when an incompatible SMIRNOFF aromaticity model is checked for compatibility.
    """
    pass

class ParseError(MessageException):
    """
    Error for when a SMIRNOFF data structure is not parseable by a ForceField
    """
    pass




#=============================================================================================
# FORCEFIELD
#=============================================================================================

# QUESTION: How should we document private object fields?

# TODO: How do we serialize/deserialize `ForceField`'s object model? Can we rely on pickle?


class ForceField:
    """A factory that assigns SMIRNOFF parameters to a molecular system

    :class:`ForceField` is a factory that constructs an OpenMM :class:`simtk.openmm.System` object from a
    :class:`openforcefield.topology.Topology` object defining a (bio)molecular system containing one or more molecules.

    When a :class:`ForceField` object is created from one or more specified SMIRNOFF serialized representations,
    all :class:`ParameterHandler` subclasses currently imported are identified and registered to handle different
    sections of the SMIRNOFF force field definition file(s).

    All :class:`ParameterIOHandler` subclasses currently imported are identified and registered to handle different
    serialization formats (such as XML).

    The force field definition is processed by these handlers to populate the ``ForceField`` object model data
    structures that can easily be manipulated via the API:

    Processing a :class:`Topology` object defining a chemical system will then call all :class`ParameterHandler`
    objects in an order guaranteed to satisfy the declared processing order constraints of each
    :class`ParameterHandler`.

    Attributes
    ----------
    parameters : dict of str : list of ParameterType
        ``parameters[tagname]`` is the instantiated :class:`ParameterHandler` class that handles parameters associated
        with the force ``tagname``.
        This is the primary means of retrieving and modifying parameters, such as
        ``parameters['vdW'][0].sigma *= 1.1``
    parameter_object_handlers : dict of str : ParameterHandler class
        Registered list of :class:`ParameterHandler` classes that will handle different forcefield tags to create the parameter object model.
        ``parameter_object_handlers[tagname]`` is the :class:`ParameterHandler` that will be instantiated to process the force field definition section ``tagname``.
        :class:`ParameterHandler` classes are registered when the ForceField object is created, but can be manipulated afterwards.
    parameter_io_handlers : dict of str : ParameterIOHandler class
        Registered list of :class:`ParameterIOHandler` classes that will handle serializing/deserializing the parameter object model to string or file representations, such as XML.
        ``parameter_io_handlers[iotype]`` is the :class:`ParameterHandler` that will be instantiated to process the serialization scheme ``iotype``.
        :class:`ParameterIOHandler` classes are registered when the ForceField object is created, but can be manipulated afterwards.

    Examples
    --------

    Create a new ForceField containing the smirnoff99Frosst parameter set:

    >>> from openforcefield.typing.engines.smirnoff import ForceField
    >>> forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')

    Create an OpenMM system from a :class:`openforcefield.topology.Topology` object:

    >>> from openforcefield.topology import Molecule, Topology
    >>> ethanol = Molecule.from_smiles('CCO')
    >>> topology = Topology.from_molecules(molecules=[ethanol])
    >>> system = forcefield.create_openmm_system(topology)

    Modify the long-range electrostatics method:

    >>> forcefield.get_parameter_handler('Electrostatics').method = 'PME'

    Inspect the first few vdW parameters:

    >>> low_precedence_parameters = forcefield.get_parameter_handler('vdW').parameters[0:3]

    Retrieve the vdW parameters by SMIRKS string and manipulate it:

    >>> parameter = forcefield.get_parameter_handler('vdW').parameters['[#1:1]-[#7]']
    >>> parameter.sigma += 0.1 * unit.angstroms
    >>> parameter.epsilon *= 1.02

    Make a child vdW type more specific (checking modified SMIRKS for validity):

    >>> forcefield.get_parameter_handler('vdW').parameters[-1].smirks += '~[#53]'

    .. warning ::

       While we check whether the modified SMIRKS is still valid and has the appropriate valence type,
       we currently don't check whether the typing remains hierarchical, which could result in some types
       no longer being assignable because more general types now come *below* them and preferentially match.

    Delete a parameter:

    >>> del forcefield.get_parameter_handler('vdW').parameters['[#1:1]-[#6X4]']

    Insert a parameter at a specific point in the parameter tree:

    >>> from openforcefield.typing.engines.smirnoff import vdWHandler
    >>> new_parameter = vdWHandler.vdWType(smirks='[*:1]', epsilon=0.0157*unit.kilocalories_per_mole, rmin_half=0.6000*unit.angstroms)
    >>> forcefield.get_parameter_handler('vdW').parameters.insert(0, new_parameter)

    .. warning ::

       We currently don't check whether removing a parameter could accidentally remove the root type, so it's possible to no longer type all molecules this way.

    """

    def __init__(self,
                 *sources,
                 parameter_handler_classes=None,
                 parameter_io_handler_classes=None,
                 disable_version_check=False,
                 allow_cosmetic_attributes=False):
        """Create a new :class:`ForceField` object from one or more SMIRNOFF parameter definition files.

        Parameters
        ----------
        sources : string or file-like object or open file handle or URL (or iterable of these)
            A list of files defining the SMIRNOFF force field to be loaded.
            Currently, only `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_ is supported.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a ``read()`` method from which the forcefield XML data can be loaded.
            If multiple files are specified, any top-level tags that are repeated will be merged if they are compatible,
            with files appearing later in the sequence resulting in parameters that have higher precedence.
            Support for multiple files is primarily intended to allow solvent parameters to be specified by listing them last in the sequence.
        parameter_handler_classes : iterable of ParameterHandler classes, optional, default=None
            If not None, the specified set of ParameterHandler classes will be instantiated to create the parameter object model.
            By default, all imported subclasses of ParameterHandler are automatically registered.
        parameter_io_handler_classes : iterable of ParameterIOHandler classes
            If not None, the specified set of ParameterIOHandler classes will be used to parse/generate serialized parameter sets.
            By default, all imported subclasses of ParameterIOHandler are automatically registered.
        disable_version_check : bool, optional, default=False
            If True, will disable checks against the current highest supported forcefield version.
            This option is primarily intended for forcefield development.
        allow_cosmetic_attributes : bool, optional. Default = False
            Whether to retain non-spec kwargs from data sources.

        Examples
        --------

        Load one SMIRNOFF parameter set in XML format (searching the package data directory by default, which includes some standard parameter sets):

        >>> forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml')

        Load multiple SMIRNOFF parameter sets:

        forcefield = ForceField('test_forcefields/smirnoff99Frosst.offxml', 'test_forcefields/tip3p.offxml')

        Load a parameter set from a string:

        >>> offxml = '<SMIRNOFF version="0.2" aromaticity_model="OEAroModel_MDL"/>'
        >>> forcefield = ForceField(offxml)

        """
        # Clear all object fields
        self._initialize()

        # Store initialization options
        self.disable_version_check = disable_version_check
        # if True, we won't check which SMIRNOFF version number we're parsing

        # Register all ParameterHandler objects that will process SMIRNOFF force definitions
        # TODO: We need to change this to just find all ParameterHandler objects in this file;
        # otherwise, we can't define two different ParameterHandler subclasses to compare for a new type of energy term
        # since both will try to register themselves for the same XML tag and an Exception will be raised.
        if parameter_handler_classes is None:
            parameter_handlers = all_subclasses(ParameterHandler)
        self._register_parameter_handler_classes(parameter_handlers)

        # Register all ParameterIOHandler objects that will process serialized parameter representations
        if parameter_io_handler_classes is None:
            parameter_io_handler_classes = all_subclasses(ParameterIOHandler)
        self._register_parameter_io_handler_classes(
            parameter_io_handler_classes)

        # Parse all sources containing SMIRNOFF parameter definitions
        self.parse_sources(sources, allow_cosmetic_attributes=allow_cosmetic_attributes)

    def _initialize(self):
        """
        Initialize all object fields.
        """
        self._MAX_SUPPORTED_SMIRNOFF_VERSION = 0.2
        self._disable_version_check = False  # if True, will disable checking compatibility version
        self._aromaticity_model = DEFAULT_AROMATICITY_MODEL  # aromaticity model
        self._parameter_handler_classes = OrderedDict()  # Parameter handler classes that _can_ be initialized if needed
        self._parameter_handlers = OrderedDict()  # ParameterHandler classes to be instantiated for each parameter type
        self._parameter_io_handler_classes = OrderedDict()  # ParameterIOHandler classes that _can_ be initialiazed if needed
        self._parameter_io_handlers = OrderedDict()  # ParameterIO classes to be used for each file type
        self._parameters = ParameterList()  # ParameterHandler objects instantiated for each parameter type
        self._aromaticity_model = None


    def _check_smirnoff_version_compatibility(self, version):
        """
        Raise a parsing exception if the given file version is incompatible with this ForceField class.

        Parameters
        ----------
        version : str
            The SMIRNOFF version being read.

        Raises
        ------
        SMIRNOFFVersionError if an incompatible version is passed in.

        """
        import packaging.version
        # Use PEP-440 compliant version number comparison, if requested
        if (not self.disable_version_check) and (
                packaging.version.parse(str(version)) >
                packaging.version.parse(str(self._MAX_SUPPORTED_SMIRNOFF_VERSION))):
            raise SMIRNOFFVersionError(
                'SMIRNOFF offxml file was written with version {}, but this version of ForceField only supports '
                'up to version {}'.format(version, self._MAX_SUPPORTED_SMIRNOFF_VERSION))


    def _set_aromaticity_model(self, aromaticity_model):
        """
        Register that this forcefield is using an aromaticity model. Will check for
        compatibility with other aromaticity model(s) already in use.

        Parameters
        ----------
        aromaticity_model : str
            The aromaticity model to register.

        Raises
        ------
        SMIRNOFFAromaticityError if an incompatible aromaticity model is passed in.

        """
        # Implement better logic here if we ever support another aromaticity model
        if aromaticity_model != 'OEAroModel_MDL':
            raise SMIRNOFFAromaticityError("Read aromaticity model {}. Currently only "
                                           "OEAroModel_MDL is supported.".format(aromaticity_model))

        self._aromaticity_model = aromaticity_model


    def _register_parameter_handler_classes(self, parameter_handler_classes):
        """
        Register multiple ParameterHandler classes, ensuring they specify unique tags to process

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        parameter_handler_classes : iterable of ParameterHandler subclasses
            List of ParameterHandler classes to register for this ForceField.
        """
        for parameter_handler_class in parameter_handler_classes:
            tagname = parameter_handler_class._TAGNAME
            if tagname is not None:
                if tagname in self._parameter_handler_classes:
                    raise Exception(
                        "Attempting to register ParameterHandler {}, which provides a parser for tag"
                        " '{}', but ParameterHandler {} has already been registered to handle that tag.".format(
                        parameter_handler_class, tagname,
                        self._parameter_handler_classes[tagname])
                    )
                self._parameter_handler_classes[tagname] = parameter_handler_class

    def _register_parameter_io_handler_classes(self,
                                               parameter_io_handler_classes):
        """
        Register multiple ParameterIOHandler classes, ensuring they specify unique suffixes

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        parameter_io_handler_classes : iterable of ParameterIOHandler subclasses
            All specified ParameterIOHandler classes will be registered as ways to translate to/from the object model
            to serialized parameter sets.

        Raises
        ------
        Exception if two ParameterIOHandlers are attempted to be registered for the same file format.

        """
        for parameter_io_handler_class in parameter_io_handler_classes:
            serialization_format = parameter_io_handler_class._FORMAT
            if serialization_format is not None:
                if serialization_format in self._parameter_io_handler_classes.keys(
                ):
                    raise Exception(
                        "Attempting to register ParameterIOHandler {}, which provides a IO parser for format "
                        "'{}', but ParameterIOHandler {} has already been registered to handle that tag.".format(
                        parameter_io_handler_class, serialization_format,
                        self._parameter_io_handler_classes[
                        serialization_format])
                    )
                self._parameter_io_handler_classes[
                    serialization_format] = parameter_io_handler_class

    def register_parameter_handler(self, parameter_handler):
        """
        Register a new ParameterHandler for a specific tag, making it
        available for lookup in the ForceField.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        parameter_handler : A ParameterHandler object
            The ParameterHandler to register. The TAGNAME attribute of this object will be used as the key for
            registration.

        """
        tagname = parameter_handler._TAGNAME
        if tagname in self._parameter_handlers.keys():
            raise ParameterHandlerRegistrationError(
                "Tried to register parameter handler '{}' for tag '{}', but "
                "tag is already registered to {}".format(
                    parameter_handler, tagname,
                    self._parameter_handlers[tagname]))

        self._parameter_handlers[parameter_handler._TAGNAME] = parameter_handler

    def register_parameter_io_handler(self, parameter_io_handler):
        """
        Register a new ParameterIOHandler, making it available for lookup in the ForceField.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        parameter_io_handler :  A ParameterIOHandler object
            The ParameterIOHandler to register. The FORMAT attribute of this object will be used
            to associate it to a file format/suffix.

        """
        io_format = parameter_io_handler._FORMAT
        if io_format in self._parameter_io_handlers.keys():
            raise ParameterHandlerRegistrationError(
                "Tried to register parameter IO handler '{}' for tag '{}', but "
                "tag is already registered to {}".format(
                    parameter_io_handler, io_format,
                    self._parameter_io_handlers[io_format]))
        self._parameter_io_handlers[io_format] = parameter_io_handler


    # TODO: Do we want to make this optional?

    @staticmethod
    def _check_for_missing_valence_terms(name, topology, assigned_terms,
                                         topological_terms):
        """
        Check to ensure there are no missing valence terms in the given topology, identifying potential gaps in
        parameter coverage.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        name : str
            Name of the calling force Handler
        topology : openforcefield.topology.Topology
            The Topology object
        assigned_terms : iterable of ints or int tuples
            Atom index tuples defining added valence terms
        topological_terms : iterable of atoms or atom tuples
            Atom tuples defining topological valence atomsets to which forces should be assigned

        """
        # Convert assigned terms and topological terms to lists
        assigned_terms = [item for item in assigned_terms]
        topological_terms = [item for item in topological_terms]

        def ordered_tuple(atoms):
            atoms = list(atoms)
            if atoms[0] < atoms[-1]:
                return tuple(atoms)
            else:
                return tuple(reversed(atoms))

        try:
            topology_set = set([
                ordered_tuple(atom.index for atom in atomset)
                for atomset in topological_terms
            ])
            assigned_set = set([
                ordered_tuple(index for index in atomset)
                for atomset in assigned_terms
            ])
        except TypeError as te:
            topology_set = set([atom.index for atom in topological_terms])
            assigned_set = set([atomset[0] for atomset in assigned_terms])

        def render_atoms(atomsets):
            msg = ""
            for atomset in atomsets:
                msg += '%30s :' % str(atomset)
                try:
                    for atom_index in atomset:
                        atom = atoms[atom_index]
                        msg += ' %5s %3s %3s' % (atom.residue.index,
                                                 atom.residue.name, atom.name)
                except TypeError as te:
                    atom = atoms[atomset]
                    msg += ' %5s %3s %3s' % (atom.residue.index,
                                             atom.residue.name, atom.name)

                msg += '\n'
            return msg

        if set(assigned_set) != set(topology_set):
            # Form informative error message
            msg = '%s: Mismatch between valence terms added and topological terms expected.\n' % name
            atoms = [atom for atom in topology.topology_atoms]
            if len(assigned_set.difference(topology_set)) > 0:
                msg += 'Valence terms created that are not present in Topology:\n'
                msg += render_atoms(assigned_set.difference(topology_set))
            if len(topology_set.difference(assigned_set)) > 0:
                msg += 'Topological atom sets not assigned parameters:\n'
                msg += render_atoms(topology_set.difference(assigned_set))
            msg += 'topology_set:\n'
            msg += str(topology_set) + '\n'
            msg += 'assigned_set:\n'
            msg += str(assigned_set) + '\n'
            raise Exception(
                msg)  # TODO: Should we raise a more specific exception here?

    def get_parameter_handler(self, tagname, handler_kwargs=None, allow_cosmetic_attributes=False):
        """Retrieve the parameter handlers associated with the provided tagname.

        If the parameter handler has not yet been instantiated, it will be created and returned.
        If a parameter handler object already exists, it will be checked for compatibility
        and an Exception raised if it is incompatible with the provided kwargs. If compatible, the
        existing ParameterHandler will be returned.

        Parameters
        ----------
        tagname : str
            The name of the parameter to be handled.
        handler_kwargs : dict, optional. Default = None
            Dict to be passed to the handler for construction or checking compatibility. If None, will be assumed
            to represent handler defaults.
        allow_cosmetic_attributes : bool, optional. Default = False
            Whether to permit non-spec kwargs in smirnoff_data.

        Returns
        -------
        handler : An openforcefield.engines.typing.smirnoff.ParameterHandler

        Raises
        ------
        KeyError if there is no ParameterHandler for the given tagname
        """
        # If there are no kwargs for the handler, initialize handler_kwargs as an empty dict
        if handler_kwargs is None:
            handler_kwargs = dict()

        # Ensure that the ForceField has a ParameterHandler class registered that can handle this tag
        if tagname in self._parameter_handler_classes:
            ph_class = self._parameter_handler_classes[tagname]
        else:
            msg = "Cannot find a registered parameter handler class for tag '{}'\n".format(
                tagname)
            msg += "Known parameter handler class tags are {}".format(self._parameter_handler_classes.keys())
            raise KeyError(msg)

        # Initialize a new instance of this parameter handler class with the given kwargs
        new_handler = ph_class(**handler_kwargs,
                               allow_cosmetic_attributes=allow_cosmetic_attributes)

        if tagname in self._parameter_handlers:
            # If a handler of this class already exists, ensure that the two handlers encode compatible science
            old_handler = self._parameter_handlers[tagname]
            # If no handler kwargs were provided, skip the compatibility check
            if handler_kwargs != {}:
                old_handler.check_handler_compatibility(new_handler)
            return_handler = old_handler
        elif tagname in self._parameter_handler_classes:
            # Otherwise, register this handler in the forcefield
            self.register_parameter_handler(new_handler)
            return_handler = new_handler

        return return_handler

    def get_parameter_io_handler(self, io_format):
        """Retrieve the parameter handlers associated with the provided tagname.
        If the parameter IO handler has not yet been instantiated, it will be created.

        Parameters
        ----------
        io_format : str
            The name of the io format to be handled.

        Returns
        -------
        io_handler : An openforcefield.engines.typing.smirnoff.ParameterIOHandler

        Raises
        ------
        KeyError if there is no ParameterIOHandler for the given tagname
        """
        # Uppercase the format string to avoid case mismatches
        io_format = io_format.upper()

        # Remove "." (ParameterIOHandler tags do not include it)
        io_format = io_format.strip('.')

        # Find or initialize ParameterIOHandler for this format
        io_handler = None
        if io_format in self._parameter_io_handlers.keys():
            io_handler = self._parameter_io_handlers[io_format]
        elif io_format in self._parameter_io_handler_classes.keys():
            new_handler_class = self._parameter_io_handler_classes[io_format]
            io_handler = new_handler_class()
            self.register_parameter_io_handler(io_handler)
        if io_handler is None:
            msg = "Cannot find a registered parameter IO handler for format '{}'\n".format(
                io_format)
            msg += "Registered parameter IO handlers: {}\n".format(
                self._parameter_io_handlers.keys())
            raise KeyError(msg)

        return io_handler


    def parse_sources(self, sources, allow_cosmetic_attributes=True):
        """Parse a SMIRNOFF force field definition.

        Parameters
        ----------
        sources : string or file-like object or open file handle or URL (or iterable of these)
            A list of files defining the SMIRNOFF force field to be loaded.
            Currently, only `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_ is supported.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a ``read()`` method from which the forcefield XML data can be loaded.
            If multiple files are specified, any top-level tags that are repeated will be merged if they are compatible,
            with files appearing later in the sequence resulting in parameters that have higher precedence.
            Support for multiple files is primarily intended to allow solvent parameters to be specified by listing them last in the sequence.
        allow_cosmetic_attributes : bool, optional. Default = False
            Whether to permit non-spec kwargs present in the source.

        .. notes ::
           * New SMIRNOFF sections are handled independently, as if they were specified in the same file.
           * If a SMIRNOFF section that has already been read appears again, its definitions are appended to the end of the previously-read
             definitions if the sections are configured with compatible attributes; otherwise, an ``IncompatibleTagException`` is raised.

        """
        # Ensure that we are working with an iterable
        try:
            sources = iter(sources)
        except TypeError as te:
            # Make iterable object
            sources = [sources]

        # TODO: If a non-first source fails here, the forcefield might be partially modified
        for source in sources:
            smirnoff_data = self.parse_smirnoff_from_source(source)
            self._load_smirnoff_data(smirnoff_data,
                                     allow_cosmetic_attributes=allow_cosmetic_attributes)


    def _to_smirnoff_data(self, discard_cosmetic_attributes=True):
        """
        Convert this ForceField and all related ParameterHandlers to an OrderedDict representing a SMIRNOFF
        data object.

        Returns
        -------
        smirnoff_dict : OrderedDict
            A nested OrderedDict representing this ForceField as a SMIRNOFF data object.
        discard_cosmetic_attributes : bool, optional. Default=True
            Whether to discard any non-spec attributes stored in the ForceField.

        """
        l1_dict = OrderedDict()

        # Assume we will write out SMIRNOFF data in compliance with the max supported spec version
        l1_dict['version'] = self._MAX_SUPPORTED_SMIRNOFF_VERSION

        # Write out the aromaticity model used
        l1_dict['aromaticity_model'] = self._aromaticity_model

        for handler_format, parameter_handler in self._parameter_handlers.items():
            handler_tag = parameter_handler._TAGNAME
            l1_dict[handler_tag] = parameter_handler.to_dict(discard_cosmetic_attributes=discard_cosmetic_attributes)

        smirnoff_dict = OrderedDict()
        smirnoff_dict['SMIRNOFF'] = l1_dict
        return smirnoff_dict

    # TODO: Should we call this "from_dict"?
    def _load_smirnoff_data(self, smirnoff_data, allow_cosmetic_attributes=False):
        """
        Add parameters from a SMIRNOFF-format data structure to this ForceField.

        Parameters
        ----------
        smirnoff_data : OrderedDict
            A representation of a SMIRNOFF-format data structure. Begins at top-level 'SMIRNOFF' key.
        allow_cosmetic_attributes : bool, optional. Default = False
            Whether to permit non-spec kwargs in smirnoff_data.
        """

        # Ensure that SMIRNOFF is a top-level key of the dict
        if not('SMIRNOFF' in smirnoff_data):
            raise ParseError("'SMIRNOFF' must be a top-level key in the SMIRNOFF object model")

        l1_dict = smirnoff_data['SMIRNOFF']
        # Check that the aromaticity model required by this parameter set is compatible with
        # others loaded by this ForceField
        if 'aromaticity_model' in l1_dict:
            aromaticity_model = l1_dict['aromaticity_model']
            self._set_aromaticity_model(aromaticity_model)

        elif self._aromaticity_model is None:
            raise ParseError("'aromaticity_model' attribute must be specified in SMIRNOFF "
                             "tag, or contained in a previously-loaded SMIRNOFF data source")

        # Check that the SMIRNOFF version of this data structure is supported by this ForceField implementation
        if 'version' in l1_dict:
            version = l1_dict['version']
        else:
            raise ParseError("'version' attribute must be specified in SMIRNOFF tag")
        self._check_smirnoff_version_compatibility(str(version))


        # Go through the subsections, delegating each to the proper ParameterHandler

        # Define keys which are expected from the spec, but are not parameter sections
        l1_spec_keys = ['Author', 'Date', 'version', 'aromaticity_model']

        for parameter_name in l1_dict:
            # Skip (for now) cosmetic l1 items
            if parameter_name in l1_spec_keys:
                continue
            # Handle cases where a parameter name has no info (eg. ToolkitAM1BCC)
            if l1_dict[parameter_name] is None:
                handler = self.get_parameter_handler(parameter_name, {})
                continue

            # Otherwise, we expect this l1_key to correspond to a ParameterHandler
            section_dict = l1_dict[parameter_name]

            # Retrieve or create parameter handler, passing in section_dict to check for
            # compatibility if a handler for this parameter name already exists
            handler = self.get_parameter_handler(parameter_name,
                                                 section_dict,
                                                 allow_cosmetic_attributes=allow_cosmetic_attributes)
            handler._add_parameters(section_dict,
                                    allow_cosmetic_attributes=allow_cosmetic_attributes)

    def parse_smirnoff_from_source(self, source):
        """
        Reads a SMIRNOFF data structure from a source, which can be one of many types.

        Parameters
        ----------
        source : str or bytes
            sources : string or file-like object or open file handle or URL (or iterable of these)
            A list of files defining the SMIRNOFF force field to be loaded
            Currently, only `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_ is supported.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a ``read()`` method from which the forcefield XML data can be loaded.

        Returns
        -------
        smirnoff_data : OrderedDict
            A representation of a SMIRNOFF-format data structure. Begins at top-level 'SMIRNOFF' key.

        """
        from openforcefield.utils import get_data_file_path

        # Check whether this could be a file path. It could also be a
        # file handler or a simple XML string.
        if isinstance(source, str):
            # Try first the simple path.
            searched_dirs_paths = ['']
            # Then try a relative file path w.r.t. an installed directory.
            searched_dirs_paths.extend( _get_installed_offxml_dir_paths())
            # Finally, search in openforcefield/data/.
            # TODO: Remove this when smirnoff99Frosst 1.0.9 will be released.
            searched_dirs_paths.append(get_data_file_path(''))

            # Determine the actual path of the file.
            # TODO: What is desired toolkit behavior if two files with the desired name are available?
            for dir_path in searched_dirs_paths:
                file_path = os.path.join(dir_path, source)
                if os.path.isfile(file_path):
                    source = file_path
                    break

        # Process all SMIRNOFF definition files or objects
        # QUESTION: Allow users to specify forcefield URLs so they can pull forcefield definitions from the web too?
        io_formats_to_try = self._parameter_io_handler_classes.keys()

        # Parse content depending on type
        for parameter_io_format in io_formats_to_try:
            parameter_io_handler = self.get_parameter_io_handler(parameter_io_format)

            # Try parsing as a forcefield file or file-like object
            try:
                smirnoff_data = parameter_io_handler.parse_file(source)
                return smirnoff_data
            except ParseError as e:
                exception_msg = str(e)
            except (FileNotFoundError, OSError):
                # If this is not a file path or a file handle, attempt parsing as a string.
                try:
                    smirnoff_data = parameter_io_handler.parse_string(source)
                    return smirnoff_data
                except ParseError as e:
                    exception_msg = str(e)

        # If we haven't returned by now, the parsing was unsuccessful
        valid_formats = [
            input_format
            for input_format in self._parameter_io_handlers.keys()
        ]
        msg = f"Source {source} could not be read. If this is a file, ensure that the path is correct.\n"
        msg += "If the file is present, ensure it is in a known SMIRNOFF encoding.\n"
        msg += f"Valid formats are: {valid_formats}\n"
        msg += f"Parsing failed with the following error:\n{exception_msg}\n"
        raise IOError(msg)


    def to_string(self, io_format='XML', discard_cosmetic_attributes=True):
        """
        Write this Forcefield and all its associated parameters to a string in a given format which
        complies with the SMIRNOFF spec.


        Parameters
        ----------
        io_format : str or ParameterIOHandler, optional. Default='XML'
            The serialization format to write to
        discard_cosmetic_attributes : bool, default=True
            Whether to discard any non-spec attributes stored in the ForceField.

        Returns
        -------
        forcefield_string : str
            The string representation of the serialized forcefield
        """
        # Resolve which IO handler to use
        if isinstance(io_format, ParameterIOHandler):
            io_handler = io_format
        else:
            io_handler = self.get_parameter_io_handler(io_format)

        smirnoff_data = self._to_smirnoff_data(discard_cosmetic_attributes=discard_cosmetic_attributes)
        string_data = io_handler.to_string(smirnoff_data)
        return string_data

    def to_file(self, filename, io_format=None, discard_cosmetic_attributes=True):
        """
        Write this Forcefield and all its associated parameters to a string in a given format which
        complies with the SMIRNOFF spec.


        Parameters
        ----------
        filename : str
            The filename to write to
        io_format : str or ParameterIOHandler, optional. Default=None
            The serialization format to write out. If None, will attempt to be inferred from the filename.
        discard_cosmetic_attributes : bool, default=True
            Whether to discard any non-spec attributes stored in the ForceField.

        Returns
        -------
        forcefield_string : str
            The string representation of the serialized forcefield
        """

        if io_format is None:
            basename, io_format = os.path.splitext(filename)



        # Resolve which IO handler to use
        if isinstance(io_format, ParameterIOHandler):
            io_handler = io_format
        else:
            # Handle the fact that .offxml is the same as .xml
            if io_format.lower() == 'offxml' or io_format.lower() == '.offxml':
                io_format = 'xml'
            io_handler = self.get_parameter_io_handler(io_format)


        # Write out the SMIRNOFF data to the IOHandler
        smirnoff_data = self._to_smirnoff_data(discard_cosmetic_attributes=discard_cosmetic_attributes)
        io_handler.to_file(filename, smirnoff_data)


    def _resolve_parameter_handler_order(self):
        """Resolve the order in which ParameterHandler objects should execute to satisfy constraints.

        Returns
        -------
        Iterable of ParameterHandlers
            The ParameterHandlers in this ForceField, in the order that they should be called to satisfy constraints.
        """

        # Create a DAG expressing dependencies
        import networkx as nx
        G = nx.DiGraph()
        for tagname, parameter_handler in self._parameter_handlers.items():
            G.add_node(tagname)
            if parameter_handler._DEPENDENCIES is not None:
                for dependency in parameter_handler._DEPENDENCIES:
                    G.add_edge(dependency._TAGNAME, parameter_handler._TAGNAME)
        # TODO: Check to make sure DAG isn't cyclic
        # Resolve order
        ordered_parameter_handlers = list()
        for tagname in nx.topological_sort(G):
            if tagname in self._parameter_handlers:
                ordered_parameter_handlers.append(
                    self._parameter_handlers[tagname])
            else:
                # TODO: Is it safe to pass "{}" as the handler_kwargs? If the handler doesn't exist, do we want to assume default values?
                ordered_parameter_handlers.append(
                    self.get_parameter_handler(tagname, {}))
        #ordered_parameter_handlers = [ self.get_parameter_handler(tagname, {}) for tagname in nx.topological_sort(G) ]
        return ordered_parameter_handlers

    # TODO: Should we add convenience methods to parameterize a Topology and export directly to AMBER, gromacs, CHARMM, etc.?
    #       Or should we create an "enhanced" openforcefield System object that knows how to convert to all of these formats?
    #       We could even create a universal applyParameters(format='AMBER') method that allows us to export to whatever system we want.

    # TODO: Should the Topology contain the default box vectors? Or should we require they be specified externally?

    # TODO: How do we know if the system is periodic or not?
    # TODO: Should we also accept a Molecule as an alternative to a Topology?



    def create_openmm_system(self,
                             topology,
                             **kwargs):
        """Create an OpenMM System representing the interactions for the specified Topology with the current force field

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The ``Topology`` corresponding to the system to be parameterized
        charge_from_molecules : List[openforcefield.molecule.Molecule], optional
            If specified, partial charges will be taken from the given molecules
            instead of being determined by the force field.

        Returns
        -------
        system : simtk.openmm.System
            The newly created OpenMM System corresponding to the specified ``topology``

        """
        # Make a deep copy of the topology so we don't accidentally modify it
        topology = copy.deepcopy(topology)

        # Set the topology aromaticity model to that used by the current forcefield
        # TODO: See openforcefield issue #206 for proposed implementation of aromaticity
        #topology.set_aromaticity_model(self._aromaticity_model)

        # Create an empty OpenMM System
        system = openmm.System()

        # Set periodic boundary conditions if specified
        if topology.box_vectors is not None:
            system.setDefaultPeriodicBoxVectors(*topology.box_vectors)

        # Add particles (both atoms and virtual sites) with appropriate masses
        for atom in topology.topology_particles:
            system.addParticle(atom.atom.mass)

        # Determine the order in which to process ParameterHandler objects in order to satisfy dependencies
        parameter_handlers = self._resolve_parameter_handler_order()

        # Check if any kwargs have been provided that aren't handled by force Handlers
        # TODO: Delete this and kwargs from arguments above?
        known_kwargs = set()
        for parameter_handler in parameter_handlers:
            known_kwargs.update(parameter_handler.known_kwargs)
        unknown_kwargs = set(kwargs.keys()).difference(known_kwargs)
        if len(unknown_kwargs) > 0:
            msg = "The following keyword arguments to create_openmm_system() are not used by any registered force Handler: {}\n".format(
                unknown_kwargs)
            msg += "Known keyword arguments: {}".format(known_kwargs)
            raise ValueError(msg)

        # Add forces and parameters to the System
        for parameter_handler in parameter_handlers:
            parameter_handler.create_force(system, topology, **kwargs)

        # Let force Handlers do postprocessing
        for parameter_handler in parameter_handlers:
            parameter_handler.postprocess_system(system, topology, **kwargs)

        return system

    def create_parmed_structure(self,
                                topology,
                                positions,
                                **kwargs):
        """Create a ParmEd Structure object representing the interactions for the specified Topology with the current force field

        This method creates a `ParmEd <http://github.com/parmed/parmed>`_ ``Structure`` object containing a topology, positions, and parameters.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The ``Topology`` corresponding to the ``System`` object to be created.
        positions : simtk.unit.Quantity of dimension (natoms,3) with units compatible with angstroms
            The positions corresponding to the ``System`` object to be created

        Returns
        -------
        structure : parmed.Structure
            The newly created ``parmed.Structure`` object

        """
        raise NotImplementedError
        #import parmed
        # TODO: Automagically handle expansion of virtual sites? Or is Topology supposed to do that?

        # Create OpenMM System
        #system = self.create_openmm_system(
        #    topology, **kwargs)

        # Create a ParmEd Structure object
        #structure = parmed.openmm.topsystem.load_topology(
        #    topology.to_openmm(), system, positions)
        #
        #return structure

    def label_molecules(self, topology):
        """Return labels for a list of molecules corresponding to parameters from this force field.
        For each molecule, a dictionary of force types is returned, and for each force type,
        each force term is provided with the atoms involved, the parameter id assigned, and the corresponding SMIRKS.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            A Topology object containing one or more unique molecules to be labeled

        Returns
        -------
        molecule_labels : list
            List of labels for unique molecules. Each entry in the list corresponds to
            one unique molecule in the Topology and is a dictionary keyed by force type,
            i.e., ``molecule_labels[0]['HarmonicBondForce']`` gives details for the harmonic
            bond parameters for the first molecule. Each element is a list of the form:
            ``[ ( [ atom1, ..., atomN], parameter_id, SMIRKS), ... ]``.

        .. todo ::

           What is the most useful API for this method?
           Should we instead accept :class:`Molecule` objects as input and individually return labels?
           Should we attach the labels to the :class:`Molecule` object?
           Or should we label all interactions in a :class:`Topology` instead of just labeling its ``unique_molecules``?

        """
        from openforcefield.topology import Topology
        # Loop over molecules and label
        molecule_labels = list()
        for molecule_idx, molecule in enumerate(topology.reference_molecules):
            top_mol = Topology.from_molecules([molecule])
            current_molecule_labels = dict()
            for tag, parameter_handler in self._parameter_handlers.items():
                matches = parameter_handler.find_matches(top_mol)
                current_molecule_labels[tag] = matches

            molecule_labels.append(current_molecule_labels)
        return molecule_labels
