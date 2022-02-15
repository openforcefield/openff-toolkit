#!/usr/bin/env python

# =============================================================================================
# MODULE DOCSTRING
# =============================================================================================
"""
Parameter assignment tools for the SMIRNOFF (SMIRKS Native Open Force Field) format.

.. codeauthor:: John D. Chodera <john.chodera@choderalab.org>
.. codeauthor:: David L. Mobley <dmobley@mobleylab.org>
.. codeauthor:: Peter K. Eastman <peastman@stanford.edu>

.. todo ::

   * Speed up overall import time by putting non-global imports only where they are needed

"""

__all__ = [
    "get_available_force_fields",
    "MAX_SUPPORTED_VERSION",
    "ParameterHandlerRegistrationError",
    "SMIRNOFFVersionError",
    "SMIRNOFFAromaticityError",
    "SMIRNOFFParseError",
    "PartialChargeVirtualSitesError",
    "ForceField",
]

import copy
import logging
import os
import pathlib
import warnings
from collections import OrderedDict
from typing import TYPE_CHECKING, List

from openff.toolkit.topology.molecule import DEFAULT_AROMATICITY_MODEL
from openff.toolkit.typing.engines.smirnoff.io import ParameterIOHandler
from openff.toolkit.typing.engines.smirnoff.parameters import (
    IncompatibleParameterError,
    ParameterHandler,
)
from openff.toolkit.typing.engines.smirnoff.plugins import load_handler_plugins
from openff.toolkit.utils.exceptions import (
    ParameterHandlerRegistrationError,
    PartialChargeVirtualSitesError,
    SMIRNOFFAromaticityError,
    SMIRNOFFParseError,
    SMIRNOFFVersionError,
)
from openff.toolkit.utils.utils import (
    all_subclasses,
    convert_0_1_smirnoff_to_0_2,
    convert_0_2_smirnoff_to_0_3,
    convert_all_quantities_to_string,
    convert_all_strings_to_quantity,
    requires_package,
)

if TYPE_CHECKING:
    from openff.toolkit.topology import Topology

deprecated_names = ["ParseError"]


def __getattr__(name):
    if name in deprecated_names:
        warnings.filterwarnings("default", category=DeprecationWarning)
        warning_msg = f"{name} is DEPRECATED and will be removed in a future release of the OpenFF Toolkit."
        warnings.warn(warning_msg, DeprecationWarning)

        if name == "ParseError":
            from openff.toolkit.utils.exceptions import _DeprecatedParseError

            return _DeprecatedParseError

    raise AttributeError(f"module {__name__} has no attribute {name}")


# =============================================================================================
# CONFIGURE LOGGER
# =============================================================================================

logger = logging.getLogger(__name__)

# =============================================================================================
# PRIVATE METHODS
# =============================================================================================

# Directory paths used by ForceField to discover offxml files.
_installed_offxml_dir_paths: List[str] = []


def _get_installed_offxml_dir_paths() -> List[str]:
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
        for entry_point in iter_entry_points(
            group="openforcefield.smirnoff_forcefield_directory"
        ):
            _installed_offxml_dir_paths.extend(entry_point.load()())
    return _installed_offxml_dir_paths


def get_available_force_fields(full_paths=False):
    """
    Get the filenames of all available .offxml force field files.

    Availability is determined by what is discovered through the
    ``openforcefield.smirnoff_forcefield_directory`` entry point. If the
    ``openff-forcefields`` package is installed, this should include several
    .offxml files such as ``openff-1.0.0.offxml``\ .

     Parameters
     ----------
     full_paths : bool, default=False
         If False, return the name of each available \*.offxml file.
         If True, return the full path to each available \*.offxml file.

     Returns
     -------
     available_force_fields : List[str]
         List of available force field files

    """
    installed_paths = _get_installed_offxml_dir_paths()
    available_force_fields = []
    for installed_path in installed_paths:
        for globbed in pathlib.Path(installed_path).rglob("*.offxml"):
            if full_paths:
                available_force_fields.append(globbed.as_posix())
            else:
                available_force_fields.append(globbed.name)
    return available_force_fields


# TODO: Instead of having a global version number, alow each Force to have a separate version number
MAX_SUPPORTED_VERSION = (
    "1.0"  # maximum version of the SMIRNOFF spec supported by this SMIRNOFF force field
)


# =============================================================================================
# FORCEFIELD
# =============================================================================================

# QUESTION: How should we document private object fields?

# TODO: How do we serialize/deserialize `ForceField`'s object model? Can we rely on pickle?


class ForceField:
    """A factory that assigns SMIRNOFF parameters to a molecular system

    :class:`ForceField` is a factory that constructs an OpenMM :class:`openmm.System` object from a
    :class:`openff.toolkit.topology.Topology` object defining a (bio)molecular system containing one or more molecules.

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
        Registered list of :class:`ParameterHandler` classes that will handle different force field tags to create the parameter object model.
        ``parameter_object_handlers[tagname]`` is the :class:`ParameterHandler` that will be instantiated to process the force field definition section ``tagname``.
        :class:`ParameterHandler` classes are registered when the ForceField object is created, but can be manipulated afterwards.
    parameter_io_handlers : dict of str : ParameterIOHandler class
        Registered list of :class:`ParameterIOHandler` classes that will handle serializing/deserializing the parameter object model to string or file representations, such as XML.
        ``parameter_io_handlers[iotype]`` is the :class:`ParameterHandler` that will be instantiated to process the serialization scheme ``iotype``.
        :class:`ParameterIOHandler` classes are registered when the ForceField object is created, but can be manipulated afterwards.

    Examples
    --------

    Create a new ForceField containing the smirnoff99Frosst parameter set:

    >>> from openff.toolkit.typing.engines.smirnoff import ForceField
    >>> forcefield = ForceField('test_forcefields/test_forcefield.offxml')

    Create an OpenMM system from a :class:`openff.toolkit.topology.Topology` object:

    >>> from openff.toolkit.topology import Molecule, Topology
    >>> ethanol = Molecule.from_smiles('CCO')
    >>> topology = Topology.from_molecules(molecules=[ethanol])
    >>> system = forcefield.create_openmm_system(topology)

    Modify the long-range electrostatics method:

    >>> forcefield.get_parameter_handler('Electrostatics').method = 'PME'

    Inspect the first few vdW parameters:

    >>> low_precedence_parameters = forcefield.get_parameter_handler('vdW').parameters[0:3]

    Retrieve the vdW parameters by SMIRKS string and manipulate it:

    >>> parameter = forcefield.get_parameter_handler('vdW').parameters['[#1:1]-[#7]']
    >>> parameter.rmin_half += 0.1 * unit.angstroms
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

    >>> from openff.toolkit.typing.engines.smirnoff import vdWHandler
    >>> new_parameter = vdWHandler.vdWType(smirks='[*:1]', epsilon=0.0157*unit.kilocalories_per_mole, rmin_half=0.6000*unit.angstroms)
    >>> forcefield.get_parameter_handler('vdW').parameters.insert(0, new_parameter)

    .. warning ::

       We currently don't check whether removing a parameter could accidentally remove the root type, so it's possible to no longer type all molecules this way.

    """

    def __init__(
        self,
        *sources,
        aromaticity_model=DEFAULT_AROMATICITY_MODEL,
        parameter_handler_classes=None,
        parameter_io_handler_classes=None,
        disable_version_check=False,
        allow_cosmetic_attributes=False,
        load_plugins=False,
    ):
        """Create a new :class:`ForceField` object from one or more SMIRNOFF parameter definition files.

        Parameters
        ----------
        sources : string or file-like object or open file handle or URL (or iterable of these)
            A list of files defining the SMIRNOFF force field to be loaded.
            Currently, only `the SMIRNOFF XML format <https://openforcefield.github.io/standards/standards/smirnoff/>`_ is supported.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a ``read()`` method from which the force field XML data can be loaded.
            If multiple files are specified, any top-level tags that are repeated will be merged if they are compatible,
            with files appearing later in the sequence resulting in parameters that have higher precedence.
            Support for multiple files is primarily intended to allow solvent parameters to be specified by listing them last in the sequence.
        aromaticity_model : string, default='OEAroModel_MDL'
            The aromaticity model used by the force field. Currently, only 'OEAroModel_MDL' is supported
        parameter_handler_classes : iterable of ParameterHandler classes, optional, default=None
            If not None, the specified set of ParameterHandler classes will be instantiated to create the parameter object model.
            By default, all imported subclasses of ParameterHandler are automatically registered.
        parameter_io_handler_classes : iterable of ParameterIOHandler classes
            If not None, the specified set of ParameterIOHandler classes will be used to parse/generate serialized parameter sets.
            By default, all imported subclasses of ParameterIOHandler are automatically registered.
        disable_version_check : bool, optional, default=False
            If True, will disable checks against the current highest supported force field version.
            This option is primarily intended for force field development.
        allow_cosmetic_attributes : bool, optional. Default = False
            Whether to retain non-spec kwargs from data sources.
        load_plugins: bool, optional. Default = False
            Whether to load ``ParameterHandler`` classes which have been registered
            by installed plugins.

        Examples
        --------

        Load one SMIRNOFF parameter set in XML format (searching the package data directory by default, which includes some standard parameter sets):

        >>> forcefield = ForceField('test_forcefields/test_forcefield.offxml')

        Load multiple SMIRNOFF parameter sets:

        forcefield = ForceField('test_forcefields/test_forcefield.offxml', 'test_forcefields/tip3p.offxml')

        Load a parameter set from a string:

        >>> offxml = '<SMIRNOFF version="0.2" aromaticity_model="OEAroModel_MDL"/>'
        >>> forcefield = ForceField(offxml)

        """
        # Clear all object fields
        self._initialize()

        self.aromaticity_model = aromaticity_model
        # Store initialization options
        self.disable_version_check = disable_version_check
        # if True, we won't check which SMIRNOFF version number we're parsing

        # Register all ParameterHandler objects that will process SMIRNOFF force definitions
        # TODO: We need to change this to just find all ParameterHandler objects in this file;
        # otherwise, we can't define two different ParameterHandler subclasses to compare for a new type of energy term
        # since both will try to register themselves for the same XML tag and an Exception will be raised.
        if parameter_handler_classes is None:
            parameter_handler_classes = all_subclasses(ParameterHandler)
        if load_plugins:

            registered_handlers = load_handler_plugins()

            # Make sure the same handlers aren't added twice.
            parameter_handler_classes += [
                handler
                for handler in registered_handlers
                if handler not in parameter_handler_classes
            ]

        self._register_parameter_handler_classes(parameter_handler_classes)

        # Register all ParameterIOHandler objects that will process serialized parameter representations
        if parameter_io_handler_classes is None:
            parameter_io_handler_classes = all_subclasses(ParameterIOHandler)

        self._register_parameter_io_handler_classes(parameter_io_handler_classes)

        # Parse all sources containing SMIRNOFF parameter definitions
        self.parse_sources(sources, allow_cosmetic_attributes=allow_cosmetic_attributes)

    def _initialize(self):
        """
        Initialize all object fields.
        """
        self._MIN_SUPPORTED_SMIRNOFF_VERSION = 0.1
        self._MAX_SUPPORTED_SMIRNOFF_VERSION = 0.3
        self._disable_version_check = (
            False  # if True, will disable checking compatibility version
        )
        self._aromaticity_model = None
        self._parameter_handler_classes = (
            OrderedDict()
        )  # Parameter handler classes that _can_ be initialized if needed
        self._parameter_handlers = (
            OrderedDict()
        )  # ParameterHandler classes to be instantiated for each parameter type
        self._parameter_io_handler_classes = (
            OrderedDict()
        )  # ParameterIOHandler classes that _can_ be initialiazed if needed
        self._parameter_io_handlers = (
            OrderedDict()
        )  # ParameterIO classes to be used for each file type
        self._author = None
        self._date = None

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
            (
                packaging.version.parse(str(version))
                > packaging.version.parse(str(self._MAX_SUPPORTED_SMIRNOFF_VERSION))
            )
            or (
                packaging.version.parse(str(version))
                < packaging.version.parse(str(self._MIN_SUPPORTED_SMIRNOFF_VERSION))
            )
        ):
            raise SMIRNOFFVersionError(
                "SMIRNOFF offxml file was written with version {}, but this version of ForceField only supports "
                "version {} to version {}".format(
                    version,
                    self._MIN_SUPPORTED_SMIRNOFF_VERSION,
                    self._MAX_SUPPORTED_SMIRNOFF_VERSION,
                )
            )

    @property
    def aromaticity_model(self):
        """Returns the aromaticity model for this ForceField object.

        Returns
        -------
        aromaticity_model
            The aromaticity model for this force field.
        """
        return self._aromaticity_model

    @aromaticity_model.setter
    def aromaticity_model(self, aromaticity_model):
        """
        Register that this force field is using an aromaticity model. Will check for
        compatibility with other aromaticity model(s) already in use.

        Parameters
        ----------
        aromaticity_model : str
            The aromaticity model to register.

        Raises
        ------
        SMIRNOFFAromaticityError if an incompatible aromaticity model is passed in.

        .. notes ::
           * Currently, the only supported aromaticity model is 'OEAroModel_MDL'.

        """
        # Implement better logic here if we ever support another aromaticity model
        if aromaticity_model != "OEAroModel_MDL":
            raise SMIRNOFFAromaticityError(
                "Read aromaticity model {}. Currently only "
                "OEAroModel_MDL is supported.".format(aromaticity_model)
            )

        self._aromaticity_model = aromaticity_model

    def _add_author(self, author):
        """
        Add an author to this force field. If this functional is called multiple times, all provided authors
        will be concatenated with the string " AND ". No redundancy checking is performed by this function.

        Parameters
        ----------
        author : str
            The author to add to this ForceField object
        """
        if self._author is None:
            self._author = author
        else:
            self._author += " AND " + author

    def _add_date(self, date):
        """
        Add an date to this force field. If this functional is called multiple times, all provided dates
        will be concatenated with the string " AND ". No redundancy checking is performed by this function.

        Parameters
        ----------
        date : str
            The author to add to this ForceField object
        """
        if self._date is None:
            self._date = date
        else:
            self._date += " AND " + date

    @property
    def author(self):
        """Returns the author data for this ForceField object. If not defined in any loaded files, this will be None.

        Returns
        -------
        author : str
            The author data for this force field.
        """
        return self._author

    @author.setter
    def author(self, author):
        """Set the author data for this ForceField object. If not defined in any loaded files, this will be None.

        Parameters
        ----------
        author : str
            The author data to set for this force field.
        """
        self._author = author

    @property
    def date(self):
        """Returns the date data for this ForceField object. If not defined in any loaded files, this will be None.

        Returns
        -------
        date : str
            The date data for this force field.
        """
        return self._date

    @date.setter
    def date(self, date):
        """Set the author data for this ForceField object. If not defined in any loaded files, this will be None.

        Parameters
        ----------
        date : str
            The date data to set for this force field.
        """
        self._date = date

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
                            parameter_handler_class,
                            tagname,
                            self._parameter_handler_classes[tagname],
                        )
                    )
                self._parameter_handler_classes[tagname] = parameter_handler_class

    def _register_parameter_io_handler_classes(self, parameter_io_handler_classes):
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
                if serialization_format in self._parameter_io_handler_classes.keys():
                    raise Exception(
                        "Attempting to register ParameterIOHandler {}, which provides a IO parser for format "
                        "'{}', but ParameterIOHandler {} has already been registered to handle that tag.".format(
                            parameter_io_handler_class,
                            serialization_format,
                            self._parameter_io_handler_classes[serialization_format],
                        )
                    )
                self._parameter_io_handler_classes[
                    serialization_format
                ] = parameter_io_handler_class

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
                    parameter_handler, tagname, self._parameter_handlers[tagname]
                )
            )

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
                    parameter_io_handler,
                    io_format,
                    self._parameter_io_handlers[io_format],
                )
            )
        self._parameter_io_handlers[io_format] = parameter_io_handler

    @property
    def registered_parameter_handlers(self):
        """
        Return the list of registered parameter handlers by name

        .. warning :: This API is experimental and subject to change.

        Returns
        -------
            registered_parameter_handlers: iterable of names of ParameterHandler objects in this ForceField

        """
        return [*self._parameter_handlers.keys()]

    # TODO: Do we want to make this optional?

    @staticmethod
    def _check_for_missing_valence_terms(
        name, topology, assigned_terms, topological_terms
    ):
        """
        Check to ensure there are no missing valence terms in the given topology, identifying potential gaps in
        parameter coverage.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        name : str
            Name of the calling force Handler
        topology : openff.toolkit.topology.Topology
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
            topology_set = set(
                [
                    ordered_tuple(atom.index for atom in atomset)
                    for atomset in topological_terms
                ]
            )
            assigned_set = set(
                [
                    ordered_tuple(index for index in atomset)
                    for atomset in assigned_terms
                ]
            )
        except TypeError as te:
            topology_set = set([atom.index for atom in topological_terms])
            assigned_set = set([atomset[0] for atomset in assigned_terms])

        def render_atoms(atomsets):
            msg = ""
            for atomset in atomsets:
                msg += f"{atomset:30} :"
                try:
                    for atom_index in atomset:
                        atom = atoms[atom_index]
                        msg += f" {atom.residue.index:5} {atom.residue.name:3} {atom.name:3}"
                except TypeError as te:
                    atom = atoms[atomset]
                    msg += (
                        f" {atom.residue.index:5} {atom.residue.name:3} {atom.name:3}"
                    )

                msg += "\n"
            return msg

        if set(assigned_set) != set(topology_set):
            # Form informative error message
            msg = f"{name}: Mismatch between valence terms added and topological terms expected.\n"
            atoms = [atom for atom in topology.topology_atoms]
            if len(assigned_set.difference(topology_set)) > 0:
                msg += "Valence terms created that are not present in Topology:\n"
                msg += render_atoms(assigned_set.difference(topology_set))
            if len(topology_set.difference(assigned_set)) > 0:
                msg += "Topological atom sets not assigned parameters:\n"
                msg += render_atoms(topology_set.difference(assigned_set))
            msg += "topology_set:\n"
            msg += str(topology_set) + "\n"
            msg += "assigned_set:\n"
            msg += str(assigned_set) + "\n"
            raise Exception(
                msg
            )  # TODO: Should we raise a more specific exception here?

    def get_parameter_handler(
        self, tagname, handler_kwargs=None, allow_cosmetic_attributes=False
    ):
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
            Dict to be passed to the handler for construction or checking compatibility. If this is None and no
            existing ParameterHandler exists for the desired tag, a handler will be initialized with all default
            values. If this is None and a handler for the desired tag exists, the existing ParameterHandler will
            be returned.
        allow_cosmetic_attributes : bool, optional. Default = False
            Whether to permit non-spec kwargs in smirnoff_data.

        Returns
        -------
        handler : An openff.toolkit.engines.typing.smirnoff.ParameterHandler

        Raises
        ------
        KeyError if there is no ParameterHandler for the given tagname
        """
        # If there are no kwargs for the handler, initialize handler_kwargs as an empty dict
        skip_version_check = False
        if handler_kwargs is None:
            handler_kwargs = dict()
            skip_version_check = True

        # Ensure that the ForceField has a ParameterHandler class registered that can handle this tag
        ph_class = self._get_parameter_handler_class(tagname)

        if tagname in self._parameter_handlers:
            # If a handler of this class already exists, ensure that the two handlers encode compatible science
            old_handler = self._parameter_handlers[tagname]
            # If no handler kwargs were provided, skip the compatibility check
            if handler_kwargs != {}:
                # Initialize a new instance of this parameter handler class with the given kwargs
                new_handler = ph_class(
                    **handler_kwargs,
                    allow_cosmetic_attributes=allow_cosmetic_attributes,
                    skip_version_check=skip_version_check,
                )
                old_handler.check_handler_compatibility(new_handler)
            return_handler = old_handler
        elif tagname in self._parameter_handler_classes:
            # Otherwise, register this handler in the force field
            # Initialize a new instance of this parameter handler class with the given kwargs
            new_handler = ph_class(
                **handler_kwargs,
                allow_cosmetic_attributes=allow_cosmetic_attributes,
                skip_version_check=skip_version_check,
            )
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
        io_handler : An openff.toolkit.engines.typing.smirnoff.ParameterIOHandler

        Raises
        ------
        KeyError if there is no ParameterIOHandler for the given tagname
        """
        # Uppercase the format string to avoid case mismatches
        io_format = io_format.upper()

        # Remove "." (ParameterIOHandler tags do not include it)
        io_format = io_format.strip(".")

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
                io_format
            )
            msg += "Registered parameter IO handlers: {}\n".format(
                self._parameter_io_handlers.keys()
            )
            raise KeyError(msg)

        return io_handler

    def deregister_parameter_handler(self, handler):
        """
        Deregister a parameter handler specified by tag name, class, or instance.

        Parameters
        ----------
        handler: str, openff.toolkit.typing.engines.smirnoff.ParameterHandler-derived type or object
            The handler to deregister.
        """
        if isinstance(handler, ParameterHandler):
            tagname = handler.TAGNAME
        elif isinstance(
            handler, str
        ):  # Catch case of name (as str) before checking subclass
            tagname = handler
        elif issubclass(handler, ParameterHandler):
            tagname = handler._TAGNAME
        else:
            tagname = handler
        del self._parameter_handlers[tagname]

    def parse_sources(self, sources, allow_cosmetic_attributes=True):
        """Parse a SMIRNOFF force field definition.

        Parameters
        ----------
        sources : string or file-like object or open file handle or URL (or iterable of these)
            A list of files defining the SMIRNOFF force field to be loaded.
            Currently, only `the SMIRNOFF XML format <https://openforcefield.github.io/standards/standards/smirnoff/>`_ is supported.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a ``read()`` method from which the force field XML data can be loaded.
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

        # TODO: If a non-first source fails here, the force field might be partially modified
        for source in sources:
            smirnoff_data = self.parse_smirnoff_from_source(source)
            self._load_smirnoff_data(
                smirnoff_data, allow_cosmetic_attributes=allow_cosmetic_attributes
            )

    def _to_smirnoff_data(self, discard_cosmetic_attributes=False):
        """
        Convert this ForceField and all related ParameterHandlers to an OrderedDict representing a SMIRNOFF
        data object.

        Returns
        -------
        smirnoff_dict : OrderedDict
            A nested OrderedDict representing this ForceField as a SMIRNOFF data object.
        discard_cosmetic_attributes : bool, optional. Default=False
            Whether to discard any non-spec attributes stored in the ForceField.

        """
        l1_dict = OrderedDict()

        # Assume we will write out SMIRNOFF data in compliance with the max supported spec version
        l1_dict["version"] = self._MAX_SUPPORTED_SMIRNOFF_VERSION

        # Write out the aromaticity model used
        l1_dict["aromaticity_model"] = self._aromaticity_model

        # Write out author and date (if they have been set)
        if not (self._author is None):
            l1_dict["Author"] = self._author

        # Write out author and date (if they have been set)
        if not (self._date is None):
            l1_dict["Date"] = self._date

        for handler_format, parameter_handler in self._parameter_handlers.items():
            handler_tag = parameter_handler._TAGNAME
            l1_dict[handler_tag] = parameter_handler.to_dict(
                discard_cosmetic_attributes=discard_cosmetic_attributes
            )

        smirnoff_dict = OrderedDict()
        smirnoff_dict["SMIRNOFF"] = l1_dict
        smirnoff_dict = convert_all_quantities_to_string(smirnoff_dict)
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
        import packaging.version

        # Check that the SMIRNOFF version of this data structure is supported by this ForceField implementation

        if "SMIRNOFF" in smirnoff_data:
            version = smirnoff_data["SMIRNOFF"]["version"]
        elif "SMIRFF" in smirnoff_data:
            version = smirnoff_data["SMIRFF"]["version"]
        else:
            raise SMIRNOFFParseError(
                "'version' attribute must be specified in SMIRNOFF tag"
            )

        self._check_smirnoff_version_compatibility(str(version))
        # Convert 0.1 spec files to 0.3 SMIRNOFF data format by converting
        # from 0.1 spec to 0.2, then 0.2 to 0.3
        if packaging.version.parse(str(version)) == packaging.version.parse("0.1"):
            # NOTE: This will convert the top-level "SMIRFF" tag to "SMIRNOFF"
            smirnoff_data = convert_0_1_smirnoff_to_0_2(smirnoff_data)
            smirnoff_data = convert_0_2_smirnoff_to_0_3(smirnoff_data)

        # Convert 0.2 spec files to 0.3 SMIRNOFF data format by removing units
        # from section headers and adding them to quantity strings at all levels.
        elif packaging.version.parse(str(version)) == packaging.version.parse("0.2"):
            smirnoff_data = convert_0_2_smirnoff_to_0_3(smirnoff_data)

        # Ensure that SMIRNOFF is a top-level key of the dict
        if not ("SMIRNOFF" in smirnoff_data):
            raise SMIRNOFFParseError(
                "'SMIRNOFF' must be a top-level key in the SMIRNOFF object model"
            )

        # Check that the aromaticity model required by this parameter set is compatible with
        # others loaded by this ForceField
        if "aromaticity_model" in smirnoff_data["SMIRNOFF"]:
            aromaticity_model = smirnoff_data["SMIRNOFF"]["aromaticity_model"]
            self.aromaticity_model = aromaticity_model

        elif self._aromaticity_model is None:
            raise SMIRNOFFParseError(
                "'aromaticity_model' attribute must be specified in SMIRNOFF "
                "tag, or contained in a previously-loaded SMIRNOFF data source"
            )

        if "Author" in smirnoff_data["SMIRNOFF"]:
            self._add_author(smirnoff_data["SMIRNOFF"]["Author"])

        if "Date" in smirnoff_data["SMIRNOFF"]:
            self._add_date(smirnoff_data["SMIRNOFF"]["Date"])

        # Go through the whole SMIRNOFF data structure, trying to convert all strings to Quantity
        smirnoff_data = convert_all_strings_to_quantity(smirnoff_data)

        # Go through the subsections, delegating each to the proper ParameterHandler

        # Define keys which are expected from the spec, but are not parameter sections
        l1_spec_keys = ["Author", "Date", "version", "aromaticity_model"]
        # TODO: Throw SMIRNOFFSpecError for unrecognized keywords

        for parameter_name in smirnoff_data["SMIRNOFF"]:
            # Skip (for now) cosmetic l1 items. They're handled above
            if parameter_name in l1_spec_keys:
                continue
            # Handle cases where a parameter name has no info (eg. ToolkitAM1BCC)
            if smirnoff_data["SMIRNOFF"][parameter_name] is None:
                # "get"ting the parameter handler here will also initialize it.
                _ = self.get_parameter_handler(parameter_name, {})
                continue

            # Otherwise, we expect this l1_key to correspond to a ParameterHandler
            section_dict = smirnoff_data["SMIRNOFF"][parameter_name]

            # TODO: Implement a ParameterHandler.from_dict() that knows how to deserialize itself for extensibility.
            #       We could let it load the ParameterTypes from the dict and append them to the existing handler
            #       after verifying that they are compatible.

            # Get the parameter types serialization that is not passed to the ParameterHandler constructor.
            ph_class = self._get_parameter_handler_class(parameter_name)
            try:
                infotype = ph_class._INFOTYPE
                parameter_list_tagname = infotype._ELEMENT_NAME
            except AttributeError:
                # The ParameterHandler doesn't have ParameterTypes (e.g. ToolkitAM1BCCHandler).
                parameter_list_dict = {}
            else:
                parameter_list_dict = section_dict.pop(parameter_list_tagname, {})

            # Must be wrapped into its own tag.
            # Assumes that parameter_list_dict is always a list

            # If the parameter list isn't empty, it must be transferred into its own tag.
            # This is necessary for deserializing SMIRNOFF force field sections which may or may
            # not have any smirks-based elements (like an empty ChargeIncrementModel section)
            if parameter_list_dict != {}:
                parameter_list_dict = {parameter_list_tagname: parameter_list_dict}

            # Retrieve or create parameter handler, passing in section_dict to check for
            # compatibility if a handler for this parameter name already exists
            handler = self.get_parameter_handler(
                parameter_name,
                section_dict,
                allow_cosmetic_attributes=allow_cosmetic_attributes,
            )
            handler._add_parameters(
                parameter_list_dict, allow_cosmetic_attributes=allow_cosmetic_attributes
            )

    def parse_smirnoff_from_source(self, source):
        """
        Reads a SMIRNOFF data structure from a source, which can be one of many types.

        Parameters
        ----------
        source : str or bytes
            sources : string or file-like object or open file handle or URL (or iterable of these)
            A list of files defining the SMIRNOFF force field to be loaded
            Currently, only `the SMIRNOFF XML format <https://openforcefield.github.io/standards/standards/smirnoff/>`_ is supported.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a ``read()`` method from which the force field XML data can be loaded.

        Returns
        -------
        smirnoff_data : OrderedDict
            A representation of a SMIRNOFF-format data structure. Begins at top-level 'SMIRNOFF' key.

        """
        from openff.toolkit.utils import get_data_file_path

        # Check whether this could be a file path. It could also be a
        # file handler or a simple XML string.
        if isinstance(source, str):
            # Try first the simple path.
            searched_dirs_paths = [""]
            # Then try a relative file path w.r.t. an installed directory.
            searched_dirs_paths.extend(_get_installed_offxml_dir_paths())
            # Finally, search in openff/toolkit/data/.
            # TODO: Remove this when smirnoff99Frosst 1.0.9 will be released.
            searched_dirs_paths.append(get_data_file_path(""))
            searched_dirs_paths.append(get_data_file_path("test_forcefields"))
            searched_dirs_paths.append(get_data_file_path("test_forcefields/old"))

            # Determine the actual path of the file.
            # TODO: What is desired toolkit behavior if two files with the desired name are available?
            for dir_path in searched_dirs_paths:
                dir_path = pathlib.Path(dir_path)
                for file_path in dir_path.glob("*.offxml"):
                    if str(file_path).lower().endswith(source.lower()):
                        source = str(file_path.absolute())
                        break

        # Process all SMIRNOFF definition files or objects
        # QUESTION: Allow users to specify force field URLs so they can pull force field definitions from the web too?
        io_formats_to_try = self._parameter_io_handler_classes.keys()

        # Parse content depending on type
        for parameter_io_format in io_formats_to_try:
            parameter_io_handler = self.get_parameter_io_handler(parameter_io_format)

            # Try parsing as a force field file or file-like object
            try:
                smirnoff_data = parameter_io_handler.parse_file(source)
                return smirnoff_data
            except SMIRNOFFParseError as e:
                exception_msg = e.msg
            except (FileNotFoundError, OSError):
                # If this is not a file path or a file handle, attempt parsing as a string.
                try:
                    smirnoff_data = parameter_io_handler.parse_string(source)
                    return smirnoff_data
                except SMIRNOFFParseError as e:
                    exception_msg = e.args[0]

        # If we haven't returned by now, the parsing was unsuccessful
        valid_formats = [
            input_format for input_format in self._parameter_io_handlers.keys()
        ]
        msg = f"Source {source} could not be read. If this is a file, ensure that the path is correct.\n"
        msg += "If the file is present, ensure it is in a known SMIRNOFF encoding.\n"
        msg += f"Valid formats are: {valid_formats}\n"
        msg += f"Parsing failed with the following error:\n{exception_msg}\n"
        raise IOError(msg)

    def to_string(self, io_format="XML", discard_cosmetic_attributes=False):
        """
        Write this Forcefield and all its associated parameters to a string in a given format which
        complies with the SMIRNOFF spec.


        Parameters
        ----------
        io_format : str or ParameterIOHandler, optional. Default='XML'
            The serialization format to write to
        discard_cosmetic_attributes : bool, default=False
            Whether to discard any non-spec attributes stored in the ForceField.

        Returns
        -------
        forcefield_string : str
            The string representation of the serialized force field
        """
        # Resolve which IO handler to use
        if isinstance(io_format, ParameterIOHandler):
            io_handler = io_format
        else:
            io_handler = self.get_parameter_io_handler(io_format)

        smirnoff_data = self._to_smirnoff_data(
            discard_cosmetic_attributes=discard_cosmetic_attributes
        )
        string_data = io_handler.to_string(smirnoff_data)
        return string_data

    def to_file(self, filename, io_format=None, discard_cosmetic_attributes=False):
        """
        Write this Forcefield and all its associated parameters to a string in a given format which
        complies with the SMIRNOFF spec.


        Parameters
        ----------
        filename : str
            The filename to write to
        io_format : str or ParameterIOHandler, optional. Default=None
            The serialization format to write out. If None, will attempt to be inferred from the filename.
        discard_cosmetic_attributes : bool, default=False
            Whether to discard any non-spec attributes stored in the ForceField.

        Returns
        -------
        forcefield_string : str
            The string representation of the serialized force field
        """

        if io_format is None:
            basename, io_format = os.path.splitext(filename)

        # Resolve which IO handler to use
        if isinstance(io_format, ParameterIOHandler):
            io_handler = io_format
        else:
            # Handle the fact that .offxml is the same as .xml
            if io_format.lower() == "offxml" or io_format.lower() == ".offxml":
                io_format = "xml"
            io_handler = self.get_parameter_io_handler(io_format)

        # Write out the SMIRNOFF data to the IOHandler
        smirnoff_data = self._to_smirnoff_data(
            discard_cosmetic_attributes=discard_cosmetic_attributes
        )
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
        for tagname, parameter_handler in self._parameter_handlers.items():
            if parameter_handler._DEPENDENCIES is not None:
                for dependency in parameter_handler._DEPENDENCIES:
                    G.add_edge(dependency._TAGNAME, parameter_handler._TAGNAME)

        # Ensure there are no loops in handler order
        if not (nx.is_directed_acyclic_graph(G)):
            raise RuntimeError(
                "Unable to resolve order in which to run ParameterHandlers. Dependencies do not form "
                "a directed acyclic graph."
            )
        # Resolve order
        ordered_parameter_handlers = list()
        for tagname in nx.topological_sort(G):
            if tagname in self._parameter_handlers:
                ordered_parameter_handlers.append(self._parameter_handlers[tagname])
        return ordered_parameter_handlers

    # TODO: Should we add convenience methods to parameterize a Topology and export directly to AMBER, gromacs, CHARMM, etc.?
    #       Or should we create an "enhanced" OpenFF System object that knows how to convert to all of these formats?
    #       We could even create a universal applyParameters(format='AMBER') method that allows us to export to whatever system we want.

    # TODO: Should the Topology contain the default box vectors? Or should we require they be specified externally?

    # TODO: How do we know if the system is periodic or not?
    # TODO: Should we also accept a Molecule as an alternative to a Topology?

    # TODO: Fall back to old code path if Interchange not installed?
    @requires_package("openmm")
    def create_openmm_system(
        self,
        topology: "Topology",
        use_interchange: bool = False,
        **kwargs,
    ):
        """Create an OpenMM System from this ForceField and a Topology.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The ``Topology`` which is to be parameterized with this ``ForceField``.
        use_interchange : bool, optional, default=True
            Whether to use the Interchange module in creation of an OpenMM ``System``.

        """
        if use_interchange:
            return self.create_interchange(topology, **kwargs,).to_openmm(
                combine_nonbonded_forces=True,
            )
        else:
            return self._old_create_openmm_system(topology, **kwargs)

    @requires_package("openmm")
    def _old_create_openmm_system(self, topology, **kwargs):
        """Create an OpenMM System representing the interactions for the specified Topology with the current force field

        Parameters
        ----------
        topology : openff.toolkit.topology.Topology
            The ``Topology`` corresponding to the system to be parameterized
        charge_from_molecules : List[openff.toolkit.molecule.Molecule], optional. default =[]
            If specified, partial charges will be taken from the given molecules
            instead of being determined by the force field.
        partial_bond_orders_from_molecules : List[openff.toolkit.molecule.Molecule], optional. default=[]
            If specified, partial bond orders will be taken from the given molecules
            instead of being determined by the force field.
            **All** bonds on each molecule given must have ``fractional_bond_order`` specified.
            A `ValueError` will be raised if any bonds have ``fractional_bond_order=None``.
            Molecules in the topology not represented in this list will have fractional
            bond orders calculated using underlying toolkits as needed.
        return_topology : bool, optional. default=False
            If ``True``, return tuple of ``(system, topology)``, where
            ``topology`` is the processed topology. Default ``False``. This topology will have the
            final partial charges assigned on its reference_molecules attribute, as well as partial
            bond orders (if they were calculated).
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry, optional. default=GLOBAL_TOOLKIT_REGISTRY
            The toolkit registry to use for operations like conformer generation and
            partial charge assignment.

        Returns
        -------
        system : openmm.System
            The newly created OpenMM System corresponding to the specified ``topology``
        topology : openff.toolkit.topology.Topology, optional.
            If the `return_topology` keyword argument is used, this method will also return a Topology. This
            can be used to inspect the partial charges and partial bond orders assigned to the molecules
            during parameterization.

        """
        import openmm
        from openff.units.openmm import to_openmm

        return_topology = kwargs.pop("return_topology", False)

        # Make a deep copy of the topology so we don't accidentally modify it
        topology = copy.deepcopy(topology)

        # set all fractional_bond_orders in topology to None
        for ref_mol in topology.reference_molecules:
            for bond in ref_mol.bonds:
                bond.fractional_bond_order = None

        # Set the topology aromaticity model to that used by the current force field
        # TODO: See openff-toolkit issue #206 for proposed implementation of aromaticity
        # topology.set_aromaticity_model(self._aromaticity_model)

        # Create an empty OpenMM System
        system = openmm.System()

        # Set periodic boundary conditions if specified
        if topology.box_vectors is not None:
            system.setDefaultPeriodicBoxVectors(*to_openmm(topology.box_vectors))

        # Add atom particles with appropriate masses
        # Virtual site particle creation is handled in the parameter handler
        # create_force call
        # This means that even though virtual sites may have been created via
        # the molecule API, an empty VirtualSites tag must exist in the FF
        for atom in topology.topology_atoms:
            # addParticle(mass.m_as(unit.dalton)) would be safer but slower
            system.addParticle(atom.mass.m)

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
                unknown_kwargs
            )
            msg += "Known keyword arguments: {}".format(known_kwargs)
            raise ValueError(msg)

        # Add forces and parameters to the System
        for parameter_handler in parameter_handlers:
            parameter_handler.create_force(system, topology, **kwargs)

        # Let force Handlers do postprocessing
        for parameter_handler in parameter_handlers:
            parameter_handler.postprocess_system(system, topology, **kwargs)

        # Handle 1-4 scaling interactions here, instead of in handlers, since OpenMM
        # does things slightly differently than other engines may
        electrostatics_14 = self.get_parameter_handler(tagname="Electrostatics").scale14
        vdw_14 = self.get_parameter_handler(tagname="vdW").scale14

        # Create exceptions based on bonds.
        # QUESTION: Will we want to do this for *all* cases, or would we ever want flexibility here?
        bond_particle_indices = []

        bond_particle_indices = [
            (topology.particle_index(bond.atom1), topology.particle_index(bond.atom2))
            for bond in topology.bonds
        ]

        # TODO: Can we generalize this to allow for `CustomNonbondedForce` implementations too?
        forces = [system.getForce(i) for i in range(system.getNumForces())]
        nonbonded_force = [f for f in forces if type(f) == openmm.NonbondedForce][0]

        nonbonded_force.createExceptionsFromBonds(
            bond_particle_indices,
            electrostatics_14,
            vdw_14,
        )

        if hasattr(self.get_parameter_handler("Electrostatics"), "cutoff"):
            vdw_cutoff = self.get_parameter_handler("vdW").cutoff
            coul_cutoff = self.get_parameter_handler("Electrostatics").cutoff
            coul_method = self.get_parameter_handler("Electrostatics").method
            if vdw_cutoff != coul_cutoff:
                if coul_method == "PME":
                    nonbonded_force.setCutoffDistance(to_openmm(vdw_cutoff))
                else:
                    raise IncompatibleParameterError(
                        "In its current implementation of the OpenFF Toolkit, with "
                        f"With electrostatics method {coul_method}, the electrostatics "
                        f"cutoff must equal the vdW cutoff. Found vdw cutoff {vdw_cutoff} "
                        f"and {coul_cutoff}."
                    )

        if return_topology:
            return (system, topology)
        else:
            return system

    @requires_package("openff.interchange")
    def create_interchange(self, topology: "Topology", box=None):
        """
        Create an Interchange object from a ForceField, Topology, and (optionally) box vectors.

        WARNING: This API and functionality are experimental and not suitable for production.

        Parameters
        ----------
        topology : openff.toolkit.topology.Topology
            The topology to create this `Interchange` object from.
        box : array-like, optional
            The box vectors or lengths to use for the Interchange object.

        Returns
        -------
        interchange : openff.interchange.Interchange
            An `Interchange` object resulting from applying this `ForceField` to a `Topology`.

        """
        from openff.interchange.components.interchange import Interchange

        return Interchange.from_smirnoff(force_field=self, topology=topology, box=box)

    def label_molecules(self, topology):
        """Return labels for a list of molecules corresponding to parameters from this force field.
        For each molecule, a dictionary of force types is returned, and for each force type,
        each force term is provided with the atoms involved, the parameter id assigned, and the corresponding SMIRKS.

        Parameters
        ----------
        topology : openff.toolkit.topology.Topology
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
        from openff.toolkit.topology import Topology
        from openff.toolkit.typing.engines.smirnoff.parameters import VirtualSiteHandler

        # Loop over molecules and label
        molecule_labels = list()
        for molecule_idx, molecule in enumerate(topology.reference_molecules):
            top_mol = Topology.from_molecules([molecule])
            current_molecule_labels = dict()
            param_is_list = False
            for tag, parameter_handler in self._parameter_handlers.items():

                if type(parameter_handler) == VirtualSiteHandler:
                    param_is_list = True

                matches = parameter_handler.find_matches(top_mol)

                # Remove the chemical environment matches from the
                # matched results.

                # Because we sometimes need to enforce atom ordering,
                # `matches` here is not a normal `dict`, but rather
                # one that transforms keys in arbitrary ways. Thus,
                # we need to make a copy of its specific class here.

                parameter_matches = matches.__class__()

                # Now make parameter_matches into a dict mapping
                # match objects to ParameterTypes

                if param_is_list:
                    for match in matches:
                        parameter_matches[match] = [
                            m.parameter_type for m in matches[match]
                        ]
                else:
                    for match in matches:
                        parameter_matches[match] = matches[match].parameter_type

                current_molecule_labels[tag] = parameter_matches

            molecule_labels.append(current_molecule_labels)
        return molecule_labels

    def _get_parameter_handler_class(self, tagname):
        """Retrieve the ParameterHandler class associated to the tagname and throw a custom error if not found."""
        try:
            ph_class = self._parameter_handler_classes[tagname]
        except KeyError:
            msg = "Cannot find a registered parameter handler class for tag '{}'\n".format(
                tagname
            )
            msg += "Known parameter handler class tags are {}".format(
                self._parameter_handler_classes.keys()
            )
            raise KeyError(msg)
        return ph_class

    def get_partial_charges(self, molecule, **kwargs):
        """Generate the partial charges for the given molecule in this force field.

        Parameters
        ----------
        molecule : :class:`openff.toolkit.topology.Molecule`
            The ``Molecule`` corresponding to the system to be parameterized
        toolkit_registry : :class:`openff.toolkit.utils.toolkits.ToolkitRegistry`, default=GLOBAL_TOOLKIT_REGISTRY
            The toolkit registry to use for operations like conformer generation and
            partial charge assignment.

        Returns
        -------
        charges : ``openmm.unit.Quantity`` with shape ``(n_atoms,)`` and dimensions of charge
            The partial charges of the provided molecule in this force field.

        Raises
        ------
        PartialChargeVirtualSitesError
            If the ``ForceField`` applies virtual sites to the ``Molecule``.
            ``get_partial_charges`` cannot identify which virtual site charges
            may belong to which atoms in this case.

        Other exceptions
            As any ``ParameterHandler`` may in principle modify charges, the entire
            force field must be applied to the molecule to produce the charges.
            Calls to this method from incorrectly or incompletely specified ``ForceField``
            objects thus may raise an exception.

        Examples
        --------

        >>> from openff.toolkit.typing.engines.smirnoff import ForceField, Molecule
        >>> ethanol = Molecule.from_smiles('CCO')
        >>> force_field = ForceField('test_forcefields/test_forcefield.offxml')

        Assign partial charges to the molecule according to the force field:

        >>> ethanol.partial_charges = force_field.get_partial_charges(ethanol)

        Use the assigned partial charges when creating an OpenMM ``System``:

        >>> topology = ethanol.to_topology()
        >>> system = forcefield.create_openmm_system(
        ...    topology,
        ...    charge_from_molecules=[ethanol]
        ... )

        This is especially useful when you want to create multiple systems
        with the same molecule or molecules, as it allows the expensive
        charge calculation to be cached.

        """
        _, top_with_charges = self.create_openmm_system(
            molecule.to_topology(), return_topology=True, **kwargs
        )

        if top_with_charges.n_topology_virtual_sites != 0:
            raise PartialChargeVirtualSitesError(
                "get_partial_charges is not supported on molecules with virtual sites"
            )

        charges = [*top_with_charges.reference_molecules][0].partial_charges
        return charges

    def __getitem__(self, val):
        """
        Syntax sugar for lookikng up a ParameterHandler. Note that only
        string-based lookups are currently supported.
        """
        if isinstance(val, str):
            if val in self._parameter_handlers:
                return self.get_parameter_handler(val)
            else:
                raise KeyError(f"Parameter handler with name '{val}' not found.")
        elif isinstance(val, ParameterHandler) or issubclass(val, ParameterHandler):
            raise NotImplementedError

    def __hash__(self):
        """Deterministically hash a ForceField object

        Notable behavior:
          * `author` and `date` are stripped from the ForceField
          * `id` and `parent_id` are stripped from each ParameterType"""

        # Completely re-constructing the force field may be overkill
        # compared to deepcopying and modifying, but is not currently slow
        ff_copy = ForceField()
        ff_copy.date = None
        ff_copy.author = None

        param_attrs_to_strip = ["_id", "_parent_id"]

        for handler_name in self.registered_parameter_handlers:
            handler = copy.deepcopy(self.get_parameter_handler(handler_name))

            for param in handler._parameters:
                for attr in param_attrs_to_strip:
                    # param.__dict__.pop(attr, None) may be faster
                    # https://stackoverflow.com/a/42303681/4248961
                    if hasattr(param, attr):
                        delattr(param, attr)

            ff_copy.register_parameter_handler(handler)

        return hash(ff_copy.to_string(discard_cosmetic_attributes=True))
