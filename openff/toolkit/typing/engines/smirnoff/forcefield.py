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
from typing import TYPE_CHECKING, List, Optional, Tuple, Union

from packaging.version import Version

from openff.toolkit.typing.engines.smirnoff.io import ParameterIOHandler
from openff.toolkit.typing.engines.smirnoff.parameters import ParameterHandler
from openff.toolkit.typing.engines.smirnoff.plugins import load_handler_plugins
from openff.toolkit.utils.constants import (
    ALLOWED_AROMATICITY_MODELS,
    DEFAULT_AROMATICITY_MODEL,
)
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
    import openmm
    from openff.units import unit

    from openff.toolkit.topology import Molecule, Topology
    from openff.toolkit.utils.base_wrapper import ToolkitWrapper
    from openff.toolkit.utils.toolkit_registry import ToolkitRegistry

logger = logging.getLogger(__name__)

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
        from importlib_metadata import entry_points

        # Find all registered entry points that should return a list of
        # paths to directories where to search for offxml files.
        for entry_point in entry_points().select(
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
    .offxml files such as ``openff-1.0.0.offxml``.

    Parameters
    ----------
    full_paths : bool, default=False
        If False, return the name of each available *.offxml file.
        If True, return the full path to each available *.offxml file.

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

    Processing a :class:`Topology` object defining a chemical system will then call all :class:`ParameterHandler`
    objects in an order guaranteed to satisfy the declared processing order constraints of each
    :class:`ParameterHandler`.

    Examples
    --------

    Create a new ForceField object from the distributed OpenFF 2.0 ("Sage") file:

    >>> from openff.toolkit import ForceField
    >>> force_field = ForceField('openff-2.0.0.offxml')

    Create an OpenMM system from a :class:`openff.toolkit.topology.Topology` object:

    >>> from openff.toolkit import Molecule, Topology
    >>> ethanol = Molecule.from_smiles('CCO')
    >>> topology = Topology.from_molecules(molecules=[ethanol])
    >>> system = force_field.create_openmm_system(topology)

    Modify the long-range electrostatics method:

    >>> force_field.get_parameter_handler('Electrostatics').periodic_potential = 'PME'

    Inspect the first few vdW parameters:

    >>> low_precedence_parameters = force_field.get_parameter_handler('vdW').parameters[0:3]

    Retrieve the vdW parameters by SMIRKS string and manipulate it:

    >>> from openff.units import unit
    >>> parameter = force_field.get_parameter_handler('vdW').parameters['[#1:1]-[#7]']
    >>> parameter.rmin_half += 0.1 * unit.angstroms
    >>> parameter.epsilon *= 1.02

    Make a child vdW type more specific (checking modified SMIRKS for validity):

    >>> force_field.get_parameter_handler('vdW').parameters[-1].smirks += '~[#53]'

    .. warning ::

       While we check whether the modified SMIRKS is still valid and has the appropriate valence type,
       we currently don't check whether the typing remains hierarchical, which could result in some types
       no longer being assignable because more general types now come *below* them and preferentially match.

    Delete a parameter:

    >>> del force_field.get_parameter_handler('vdW').parameters['[#1:1]-[#6X4]']

    Insert a parameter at a specific point in the parameter tree:

    >>> from openff.toolkit.typing.engines.smirnoff import vdWHandler
    >>> new_parameter = vdWHandler.vdWType(
    ...     smirks='[*:1]',
    ...     epsilon=0.0157*unit.kilocalories_per_mole,
    ...     rmin_half=0.6000*unit.angstroms,
    ... )
    >>> force_field.get_parameter_handler('vdW').parameters.insert(0, new_parameter)

    .. warning ::

       We currently don't check whether removing a parameter could accidentally remove the root type, so it's possible
       to no longer type all molecules this way.

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
            Currently, only `the SMIRNOFF XML format <https://openforcefield.github.io/standards/standards/smirnoff/>`_
            is supported.  Each entry may be an absolute file path, a path relative to the current working directory, a
            path relative to this module's data subdirectory (for built in force fields), or an open file-like object
            with a ``read()`` method from which the force field XML data can be loaded.  If multiple files are
            specified, any top-level tags that are repeated will be merged if they are compatible, with files appearing
            later in the sequence resulting in parameters that have higher precedence.  Support for multiple files is
            primarily intended to allow solvent parameters to be specified by listing them last in the sequence.
        aromaticity_model : str, optional, default="OEAroModel_MDL"
            The aromaticity model to use. Only OEAroModel_MDL is supported.
        parameter_handler_classes : iterable of ParameterHandler classes, optional, default=None
            If not None, the specified set of ParameterHandler classes will be instantiated to create the parameter
            object model.  By default, all imported subclasses of ParameterHandler are automatically registered.
        parameter_io_handler_classes : iterable of ParameterIOHandler classes
            If not None, the specified set of ParameterIOHandler classes will be used to parse/generate serialized
            parameter sets.  By default, all imported subclasses of ParameterIOHandler are automatically registered.
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

        Load one SMIRNOFF parameter set in XML format (searching the package data directory by default, which includes
        some standard parameter sets):

        >>> forcefield = ForceField('openff-2.0.0.offxml')

        Load multiple SMIRNOFF parameter sets:

        >>> from openff.toolkit.tests.utils import get_data_file_path
        >>> forcefield = ForceField('openff-2.0.0.offxml', get_data_file_path('test_forcefields/tip3p.offxml'))

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
            plugin_classes = load_handler_plugins()

            for handler in plugin_classes:
                if handler not in parameter_handler_classes:
                    parameter_handler_classes.append(handler)
                    self._plugin_parameter_handler_classes.append(handler)

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
        self._MIN_SUPPORTED_SMIRNOFF_VERSION = Version("0.1")
        self._MAX_SUPPORTED_SMIRNOFF_VERSION = Version("0.3")
        self._disable_version_check = (
            False  # if True, will disable checking compatibility version
        )
        self._aromaticity_model = DEFAULT_AROMATICITY_MODEL
        # Parameter handler classes that _can_ be initialized if needed
        self._parameter_handler_classes = dict()
        # ParameterHandler classes to be instantiated for each parameter type
        self._parameter_handlers = dict()
        # classes of ParameterHandlers that were registered via the plugin interface
        self._plugin_parameter_handler_classes = list()
        # ParameterIOHandler classes that _can_ be initialiazed if needed
        self._parameter_io_handler_classes = dict()
        # ParameterIO classes to be used for each file type
        self._parameter_io_handlers = dict()
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
        SMIRNOFFVersionError
            If an incompatible version is passed in.

        """
        from packaging.version import parse

        # Use PEP-440 compliant version number comparison, if requested
        if self.disable_version_check:
            pass
        else:
            if (
                parse(str(version)) > parse(str(self._MAX_SUPPORTED_SMIRNOFF_VERSION))
            ) or (
                parse(str(version)) < parse(str(self._MIN_SUPPORTED_SMIRNOFF_VERSION))
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
        SMIRNOFFAromaticityError
            If an incompatible aromaticity model is passed in.

        Notes
        -----
           * Currently, the only supported aromaticity model is 'OEAroModel_MDL'.

        """
        # Implement better logic here if we ever support another aromaticity model
        if aromaticity_model not in ALLOWED_AROMATICITY_MODELS:
            raise SMIRNOFFAromaticityError(
                f"Read aromaticity model {aromaticity_model} which is not in the set of allowed aromaticity models: "
                f"{ALLOWED_AROMATICITY_MODELS}"
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
        Exception
            If two ParameterIOHandlers are attempted to be registered for the same file format.

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
        self._parameter_handler_classes[parameter_handler._TAGNAME] = type(
            parameter_handler
        )

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
    def registered_parameter_handlers(self) -> List[str]:
        """
        Return the list of registered parameter handlers by name

        .. warning :: This API is experimental and subject to change.

        Returns
        -------
            registered_parameter_handlers: iterable of names of ParameterHandler objects in this ForceField

        """
        return [*self._parameter_handlers.keys()]

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
        KeyError
            If there is no ParameterHandler for the given tagname
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
        KeyError
            If there is no ParameterIOHandler for the given tagname
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
            msg = f"Cannot find a registered parameter IO handler for format '{io_format}'\n"
            msg += f"Registered parameter IO handlers: {self._parameter_io_handlers.keys()}\n"
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
            Currently, only `the SMIRNOFF XML format <https://openforcefield.github.io/standards/standards/smirnoff/>`_
            is supported.  Each entry may be an absolute file path, a path relative to the current working directory, a
            path relative to this module's data subdirectory (for built in force fields), or an open file-like object
            with a ``read()`` method from which the force field XML data can be loaded.  If multiple files are
            specified, any top-level tags that are repeated will be merged if they are compatible, with files appearing
            later in the sequence resulting in parameters that have higher precedence.  Support for multiple files is
            primarily intended to allow solvent parameters to be specified by listing them last in the sequence.
        allow_cosmetic_attributes : bool, optional. Default = False
            Whether to permit non-spec kwargs present in the source.

        Notes
        -----

            * New SMIRNOFF sections are handled independently, as if they were specified in the same file.
            * If a SMIRNOFF section that has already been read appears again, its definitions are appended to the end
                of the previously-read definitions if the sections are configured with compatible attributes;
                otherwise, an ``IncompatibleTagException`` is raised.

        """
        # Ensure that we are working with an iterable
        try:
            sources = iter(sources)
        except TypeError:
            # Make iterable object
            sources = [sources]

        # TODO: If a non-first source fails here, the force field might be partially modified
        for source in sources:
            smirnoff_data = self.parse_smirnoff_from_source(source)
            self._load_smirnoff_data(
                smirnoff_data, allow_cosmetic_attributes=allow_cosmetic_attributes
            )

    def _to_smirnoff_data(self, discard_cosmetic_attributes=False) -> dict:
        """
        Convert this ForceField and all related ParameterHandlers to a dict representing a SMIRNOFF
        data object.

        Parameters
        ----------
        discard_cosmetic_attributes : bool, optional. Default=False
            Whether to discard any non-spec attributes stored in the ForceField.

        Returns
        -------
        smirnoff_data : dict
            A nested dict representing this ForceField as a SMIRNOFF data object.

        """
        l1_dict = dict()

        # Assume we will write out SMIRNOFF data in compliance with the max supported spec version
        l1_dict["version"] = str(self._MAX_SUPPORTED_SMIRNOFF_VERSION)

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

        smirnoff_data = dict()
        smirnoff_data["SMIRNOFF"] = l1_dict
        smirnoff_data = convert_all_quantities_to_string(smirnoff_data)
        return smirnoff_data

    # TODO: Should we call this "from_dict"?
    def _load_smirnoff_data(self, smirnoff_data: dict, allow_cosmetic_attributes=False):
        """
        Add parameters from a SMIRNOFF-format data structure to this ForceField.

        Parameters
        ----------
        smirnoff_data : dict
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
        smirnoff_data = convert_all_strings_to_quantity(
            smirnoff_data,
            ignore_keys=["smirks", "name"],
        )

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

    def parse_smirnoff_from_source(self, source) -> dict:
        """
        Reads a SMIRNOFF data structure from a source, which can be one of many types.

        Parameters
        ----------
        source : str or bytes or file-like object
            File defining the SMIRNOFF force field to be loaded Currently, only `the SMIRNOFF XML format
            <https://openforcefield.github.io/standards/standards/smirnoff/>`_ is supported.  The file may be an
            absolute file path, a path relative to the current working directory, a path relative to this module's data
            subdirectory (for built in force fields), or an open file-like object with a ``read()`` method from which
            the force field XML data can be loaded.

        Returns
        -------
        smirnoff_data : dict
            A representation of a SMIRNOFF-format data structure. Begins at top-level 'SMIRNOFF' key.

        """
        # First, see if a file exists with a name `source` in the current directory or in directories known to the
        # plugin system. It could also be a raw XML-like string ...
        if isinstance(source, str):
            # Try first the simple path.
            searched_dirs_paths: List[str] = [os.getcwd()]
            # Then try a relative file path w.r.t. an installed directory.
            searched_dirs_paths.extend(_get_installed_offxml_dir_paths())

            # Determine the actual path of the file.
            # TODO: What is desired toolkit behavior if two files with the desired name are available?
            for dir_path in searched_dirs_paths:
                for file_path in pathlib.Path(dir_path).glob("*.offxml"):
                    if str(file_path).lower().endswith(source.lower()):
                        source = str(file_path.absolute())
                        break

        else:
            # ... or a file-like object, in which case we shouldn't look through the plugin system.
            # if it's raw bytes, no need to search for paths, either.
            searched_dirs_paths = list()

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
                exception_type = type(e)
                exception_context = "while trying to parse source as an object"
                exception_msg = e.msg
            except (FileNotFoundError, OSError):
                # If this is not a file path or a file handle, attempt parsing as a string.
                # TODO: Do we actually support parsing bytes?
                try:
                    smirnoff_data = parameter_io_handler.parse_string(source)
                    return smirnoff_data
                except SMIRNOFFParseError as e:
                    exception_type = type(e)
                    exception_context = "while trying to parse source as a file"
                    exception_msg = e.args[0]

        # If we haven't returned by now, the parsing was unsuccessful
        # There is different parsing behavior for str and file-like objects, so raise errors separately
        if isinstance(source, str):
            pretty_searched_paths = "\n    ".join(searched_dirs_paths)
            msg = (
                f"Source '{source}' could not be read. If this is a file, ensure that the path is correct.\n"
                f"Looked in the following paths and found no files named '{source}':"
                f"\n    {pretty_searched_paths}\n"
                f"If '{source}' is present as a file, ensure it is in a known SMIRNOFF encoding.\n"
                f"Valid formats are: {[*self._parameter_io_handler_classes.keys()]}\n"
                f"Parsing failed {exception_context} with the following exception and "
                f"message:\n{exception_type}\n{exception_msg}\n"
            )

        else:
            msg = (
                f"Source '{source}' could not be read.\n"
                f"Parsing failed {exception_context} with the following exception and message:\n"
                f"{exception_type}\n{exception_msg}\n"
            )

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

    # TODO: Should we also accept a Molecule as an alternative to a Topology?
    @requires_package("openmm")
    def create_openmm_system(
        self,
        topology: "Topology",
        *,
        return_topology: bool = False,
        toolkit_registry: Optional[Union["ToolkitRegistry", "ToolkitWrapper"]] = None,
        charge_from_molecules: Optional[List["Molecule"]] = None,
        partial_bond_orders_from_molecules: Optional[List["Molecule"]] = None,
        allow_nonintegral_charges: bool = False,
    ) -> Union["openmm.System", Tuple["openmm.System", "Topology"]]:
        """Create an OpenMM System from this ForceField and a Topology.

        Parameters
        ----------
        topology
            The ``Topology`` which is to be parameterized with this
            ``ForceField``.
        toolkit_registry
            The toolkit registry to use for parametrization (eg, for calculating
            partial charges and partial bond orders)
        charge_from_molecules
            Take partial charges from the input topology rather than calculating
            them. This may be useful for avoiding recalculating charges, but
            take care to ensure that your charges are appropriate for the force
            field.
        partial_bond_orders_from_molecules
            Take partial bond orders from the input topology rather than
            calculating them. This may be useful for avoiding recalculating
            PBOs, but take to ensure that they are appropriate for the
            force field.
        allow_nonintegral_charges
            Allow charges that do not sum to an integer.
        return_topology
            .. deprecated:: 0.11.0
                The ``return_topology`` argument has been deprecated and will be
                removed in v0.12.0. Call :meth:`ForceField.create_interchange`
                and take the topology from `interchange.topology` instead.

            Return the Topology with any modifications needed to parametrize it
            in a tuple along with the OpenMM system.
        """
        interchange = self.create_interchange(
            topology,
            toolkit_registry,
            charge_from_molecules=charge_from_molecules,
            partial_bond_orders_from_molecules=partial_bond_orders_from_molecules,
            allow_nonintegral_charges=allow_nonintegral_charges,
        )

        openmm_system = interchange.to_openmm(combine_nonbonded_forces=True)

        if not return_topology:
            return openmm_system
        else:
            warning_msg = (
                "The `create_openmm_system` kwarg `return_topology` is DEPRECATED and will be "
                "removed in version 0.12.0 of the OpenFF Toolkit. "
                "Use `ForceField.create_interchange` followed by `Interchange.topology`, "
                "`Interchange.to_openmm_topology`, and `Interchange.to_openmm` "
                "for long-term replacements for `return_topology` functionality."
            )
            warnings.warn(warning_msg, DeprecationWarning)
            return openmm_system, copy.deepcopy(interchange.topology)

    @requires_package("openff.interchange")
    def create_interchange(
        self,
        topology: "Topology",
        toolkit_registry: Optional[Union["ToolkitRegistry", "ToolkitWrapper"]] = None,
        charge_from_molecules: Optional[List["Molecule"]] = None,
        partial_bond_orders_from_molecules: Optional[List["Molecule"]] = None,
        allow_nonintegral_charges: bool = False,
    ):
        """
        Create an Interchange object from a ForceField, Topology, and (optionally) box vectors.

        WARNING: This API and functionality are experimental and not suitable for production.

        Parameters
        ----------
        topology : openff.toolkit.topology.Topology
            The topology to create this `Interchange` object from.
        toolkit_registry
            The toolkit registry to use for parametrization (eg, for calculating
            partial charges and partial bond orders)
        charge_from_molecules
            Take charges from the input topology rather than calculating them.
            This may be useful for avoiding recalculating charges, but take care
            to ensure that your charges are appropriate for the force field.
        partial_bond_orders_from_molecules
            Take partial bond orders from the input topology rather than
            calculating them. This may be useful for avoiding recalculating
            PBOs, but take to ensure that they are appropriate for the
            force field.
        allow_nonintegral_charges
            Allow charges that do not sum to an integer.

        Returns
        -------
        interchange : openff.interchange.Interchange
            An `Interchange` object resulting from applying this `ForceField` to a `Topology`.

        """
        from openff.interchange import Interchange  # type: ignore[import]

        from openff.toolkit.utils.toolkit_registry import _toolkit_registry_manager

        if toolkit_registry is not None:
            used_registry = toolkit_registry
        else:
            from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY

            used_registry = GLOBAL_TOOLKIT_REGISTRY

        with _toolkit_registry_manager(used_registry):
            return Interchange.from_smirnoff(
                force_field=self,
                topology=topology,
                charge_from_molecules=charge_from_molecules,
                partial_bond_orders_from_molecules=partial_bond_orders_from_molecules,
                allow_nonintegral_charges=allow_nonintegral_charges,
            )

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
           Or should we label all interactions in a :class:`Topology` instead of just labeling its
            ``unique_molecules``?

        """
        from openff.toolkit import Topology
        from openff.toolkit.typing.engines.smirnoff.parameters import VirtualSiteHandler

        # Loop over molecules and label
        molecule_labels = list()

        # TODO: This was previously ... enumerate(topology.reference_molecules). It's currently
        # unclear if this should be topology.unique_molecules instead, since that might be faster
        # (if also modifying this to label _all_ duplicates of each unique molecule)
        for molecule_idx, molecule in enumerate(topology.molecules):
            top_mol = Topology.from_molecules([molecule])
            current_molecule_labels = dict()
            for tag, parameter_handler in self._parameter_handlers.items():
                param_is_list = False

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
            msg = f"Cannot find a registered parameter handler class for tag '{tagname}'\n"
            msg += f"Known parameter handler class tags are {self._parameter_handler_classes.keys()}"
            raise KeyError(msg)
        return ph_class

    def get_partial_charges(self, molecule: "Molecule", **kwargs) -> "unit.Quantity":
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

        >>> from openff.toolkit import ForceField, Molecule
        >>> ethanol = Molecule.from_smiles('CCO')
        >>> force_field = ForceField('openff-2.0.0.offxml')

        Assign partial charges to the molecule according to the force field:

        >>> ethanol.partial_charges = force_field.get_partial_charges(ethanol)

        Use the assigned partial charges when creating an OpenMM ``System``:

        >>> topology = ethanol.to_topology()
        >>> system = force_field.create_openmm_system(
        ...    topology,
        ...    charge_from_molecules=[ethanol]
        ... )

        This is especially useful when you want to create multiple systems
        with the same molecule or molecules, as it allows the expensive
        charge calculation to be cached.

        """
        from openff.toolkit.topology.molecule import Molecule

        if not isinstance(molecule, Molecule):
            raise ValueError(
                "`molecule` argument must be a `Molecule` or subclass object. Found type "
                f"{type(molecule)}"
            )

        _, top_with_charges = self.create_openmm_system(
            molecule.to_topology(), return_topology=True, **kwargs
        )

        assert top_with_charges.n_molecules == 1, (
            "Expected a single molecule in the topology produced by Interchange. "
            f"Found {len(top_with_charges.n_molecules)} molecules."
        )

        for molecule in top_with_charges.molecules:
            return molecule.partial_charges

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
