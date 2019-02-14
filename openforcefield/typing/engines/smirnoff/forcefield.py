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

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import copy
import logging

from collections import OrderedDict

import packaging

from simtk import openmm, unit

from openforcefield.utils import all_subclasses
from openforcefield.topology import DEFAULT_AROMATICITY_MODEL
from openforcefield.typing.engines.smirnoff.parameters import ParameterList, ParameterHandler
from openforcefield.typing.engines.smirnoff.io import ParameterIOHandler

#=============================================================================================
# CONFIGURE LOGGER
#=============================================================================================

logger = logging.getLogger(__name__)

#=============================================================================================
# PRIVATE METHODS
#=============================================================================================

# TODO: Instead of having a global version number, alow each Force to have a separate version number
MAX_SUPPORTED_VERSION = '1.0'  # maximum version of the SMIRNOFF spec supported by this SMIRNOFF forcefield


class ParameterHandlerRegistrationError(Exception):
    """
    Exception for errors in ParameterHandler registration errors
    """

    def __init__(self, msg):
        super().__init__(self, msg)
        self.msg = msg


#=============================================================================================
# FORCEFIELD
#=============================================================================================

# QUESTION: How should we document private object fields?

# TODO: How do we serialize/deserialize `ForceField`'s object model? Can we rely on pickle?


class ForceField(object):
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
    >>> forcefield = ForceField('smirnoff99Frosst.offxml')

    Create an OpenMM system from a :class:`openforcefield.topology.Topology` object:

    >>> from openforcefield.topology.testsystems import WaterBox
    >>> system = forcefield.create_system(WaterBox.topology)

    Modify the long-range electrostatics method:

    >>> forcefield.forces['Electrostatics'].method = 'PME'

    Inspect the first few vdW parameters:

    >>> print(forcefield.forces['vdW'].parameters[0:3])

    Retrieve the vdW parameters by SMIRKS string and manipulate it:

    >>> parameter = forcefield.forces['vdW'].parameters['[#1:1]-[#7]']
    >>> parameter.sigma += 0.1 * unit.angstroms
    >>> parameter.epsilon *= 1.02

    Make a child vdW type more specific (checking modified SMIRKS for validity):

    >>> forcefield.forces['vdW'].parameters[-1].smirks += '$(*~[#53])'

    .. warning ::

       While we check whether the modified SMIRKS is still valid and has the appropriate valence type,
       we currently don't check whether the typing remains hierarchical, which could result in some types
       no longer being assignable because more general types now come *below* them and preferentially match.

    Delete a parameter:

    >>> del forcefield.forces['vdW'].parameters['[#1:1]-[#6X4]']

    .. warning ::

       We currently don't check whether removing a parameter could accidentally remove the root type, so it's possible to no longer type all molecules this way.

    Insert a parameter at a specific point in the parameter tree:

    >>> new_parameter = vdWType(smirks='[*:1]', epsilon=0.0157*unit.kilocalories_per_mole, rmin_half=0.6000*unit.angstroms)
    >>> forcefield.forces['vdW'].parameters.insert(0, new_parameter)

    .. warning ::

       We currently don't check whether removing a parameter could accidentally remove the root type, so it's possible to no longer type all molecules this way.

    """

    def __init__(self,
                 *sources,
                 parameter_handler_classes=None,
                 parameter_io_handler_classes=None,
                 disable_version_check=False):
        """Create a new :class:`ForceField` object from one or more SMIRNOFF parameter definition files.

        Parameters
        ----------
        sources : string or file-like object or open file handle or URL (or iterable of these)
            A list of files defining the SMIRNOFF force field to be loaded
            Currently, only `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_ is supported.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a ``read()`` method from which the forcefield XML data can be loaded.
            If multiple files are specified, any top-level tags that are repeated will be merged if they are compatible,
            with files appearing later in the sequence resulting in parameters that have higher precedence.
            Support for multiple files is primarily intended to allow solvent parameters to be specified by listing them last in the sequence.
        parameter_handler_classes : iterable of ParameterHandler classes, optional, default=None
            If not None, the specified set of ParameterHandler classes will be instantiated to create the parameter object model
            By default, all imported subclasses of ParameterHandler are automatically registered
        parameter_io_handler_classes : iterable of ParameterIOHandler classes
            If not None, the specified set of ParameterIOHandler classes will be used to parse/generate serialized parameter sets
            By default, all imported subclasses of ParameterIOHandler are automatically registered
        disable_version_check : bool, optional, default=False
            If True, will disable checks against the current highest supported forcefield version.
            This option is primarily intended for forcefield development.

        Examples
        --------

        Load one SMIRNOFF parameter set in XML format (searching the package data directory by default, which includes some standard parameter sets):

        >>> forcefield = ForceField('smirnoff99Frosst.offxml')

        Load multiple SMIRNOFF parameter sets:

        >>> forcefield = ForceField(['smirnoff99Frosst.offxml', 'tip3p.offxml'])

        Load a parameter set from a URL:

        >>> forcefield = ForceField('https://raw.githubusercontent.com/openforcefield/openforcefield/master/openforcefield/data/forcefield/')

        Load a parameter set from a string:

        >>> offxml = '<SMIRNOFF version=1.0/>'
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
        self.parse(sources)

    def _initialize(self):
        """
        Initialize all object fields.
        """
        self._disable_version_check = False  # if True, will disable checking compatibility version
        self._aromaticity_model = DEFAULT_AROMATICITY_MODEL  # aromaticity model
        self._parameter_handler_classes = OrderedDict(
        )  # Parameter handler classes that _can_ be initialized if needed
        self._parameter_handlers = OrderedDict(
        )  # ParameterHandler classes to be instantiated for each parameter type
        self._parameter_io_handler_classes = OrderedDict(
        )  # ParameterIOHandler classes that _can_ be initialiazed if needed
        self._parameter_io_handlers = OrderedDict(
        )  # ParameterIO classes to be used for each file type
        self._parameters = ParameterList(
        )  # ParameterHandler objects instantiated for each parameter type

    # TODO : This is unused and shouldn't be here -- Move to parameteriohandler
    def _parse_version(self, root):
        """
        Parse the forcefield version number and make sure it is supported.

        Parameters
        ----------
        root : etree.Element
            The document root

        """
        if 'version' in root.attrib:
            version = root.attrib['version']
            # Use PEP-440 compliant version number comparison, if requested
            if (not self.disable_version_check) and (
                    packaging.version.parse(version) >
                    packaging.version.parse(MAX_SUPPORTED_VERION)):
                self._raise_parsing_exception(
                    root,
                    'SMIRNOFF offxml file was written with version %s, but this version of ForceField only supports up to version %s'
                    % (self.version, MAX_SUPPORTED_VERSION))
        else:
            self._raise_parsing_exception(
                root, "'version' attribute must be specified in SMIRNOFF tag")

    def _register_parameter_handler_classes(self, parameter_handler_classes):
        """
        Register multiple ParameterHandler classes, ensuring they specify unique tags to process

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        parameter_handler_classes : iterable of ParameterHandler subclasses
            List of ParameterHandler classes to register for this ForceField.
        """
        #parsers = dict()
        for parameter_handler_class in parameter_handler_classes:
            #tagname = subclass._TAGNAME
            tagname = parameter_handler_class._TAGNAME
            if tagname is not None:
                if tagname in self._parameter_handler_classes:
                    raise Exception(
                        "ParameterHandler {} provides a parser for tag '{}', but ParameterHandler {} has "
                        "already been registered to handle that tag.".format(
                            parameter_handler_class, tagname,
                            self._parameter_handler_classes[tagname]))
                #parsers[tagname] = parameter_handler
                self._parameter_handler_classes[
                    tagname] = parameter_handler_class
        #return parsers

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
        #parsers = dict()
        for parameter_io_handler_class in parameter_io_handler_classes:
            serialization_format = parameter_io_handler_class._FORMAT
            if serialization_format is not None:
                if serialization_format in self._parameter_io_handler_classes.keys(
                ):
                    raise Exception(
                        "ParameterIOHandler {} provides a IO parser for format '{}', but ParameterIOHandler {} has "
                        "already been registered to handle that tag.".format(
                            parameter_io_handler_class, serialization_format,
                            self._parameter_io_handler_classes[
                                serialization_format]))
                self._parameter_io_handler_classes[
                    serialization_format] = parameter_io_handler_class
        #return parsers

    def register_parameter_handler(self, parameter_handler_class,
                                   parameter_handler_kwargs):
        """
        Register a new ParameterHandler from a specified class, instantiating the ParameterHandler object and making it
        available for lookup in the ForceField.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        parameter_handler_class : A ParameterHandler-derived object
            The ParameterHandler to register

        Returns
        -------
        new_handler : an openforcefield.engines.typing.smirnoff.ParameterHandler-derived object.
            The newly-created ParameterHandler
        """
        tagname = parameter_handler_class._TAGNAME
        if tagname in self._parameter_handlers.keys():
            raise ParameterHandlerRegistrationError(
                "Tried to register parameter handler '{}' for tag '{}', but "
                "tag is already registered to {}".format(
                    parameter_handler_class, tagname,
                    self._parameter_handlers[tagname]))
        new_handler = parameter_handler_class(self, **parameter_handler_kwargs)

        self._parameter_handlers[new_handler._TAGNAME] = new_handler
        return new_handler

    def register_parameter_io_handler(self, parameter_io_handler_class):
        """
        Register a new ParameterIOHandler from a specified class, instantiating the ParameterIOHandler object and making
        it available for lookup in the ForceField.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        parameter_io_handler_class : A subclass of ParameterIOHandler

        """
        io_format = parameter_io_handler_class._FORMAT
        if io_format in self._parameter_io_handlers.keys():
            raise ParameterHandlerRegistrationError(
                "Tried to register parameter IO handler '{}' for tag '{}', but "
                "tag is already registered to {}".format(
                    parameter_io_handler_class, io_format,
                    self._parameter_io_handlers[io_format]))
        new_io_handler = parameter_io_handler_class(self)

        self._parameter_io_handlers[io_format] = new_io_handler
        return new_io_handler

    # TODO: Do we want to make this optional?

    @staticmethod
    def _check_for_missing_valence_terms(name, topology, assigned_terms,
                                         topological_terms):
        """
        Check to ensure there are no missing valence terms in the given topology, identifying potential gaps in parameter coverage.

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

    def get_handler(self, tagname, handler_kwargs):
        """Retrieve the parameter handlers associated with the provided tagname.

        If the parameter handler has not yet been instantiated, it will be created.
        If a parameter handler object already exists, it will be checked for compatibility
        and an Exception raised if it is incompatible with the provided kwargs.

        Parameters
        ----------
        tagame : str
            The name of the parameter to be handled.
        handler_kwargs : dict
            Dict to be passed to the handler for construction or checking compatibility.

        Returns
        -------
        handler : An openforcefield.engines.typing.smirnoff.ParameterHandler

        Raises
        ------
        KeyError if there is no ParameterHandler for the given tagname
        """
        handler = None
        if tagname in self._parameter_handlers:
            # If handler already exists, make sure it is compatible
            handler = self._parameter_handlers[tagname]
            handler.check_handler_compatibility(handler_kwargs)
        elif tagname in self._parameter_handler_classes:
            new_ph_class = self._parameter_handler_classes[tagname]
            handler = self.register_parameter_handler(new_ph_class,
                                                      handler_kwargs)

        if handler is None:
            msg = "Cannot find a registered parameter handler for tag '{}'\n".format(
                tagname)
            msg += "Registered parameter handlers: {}\n".format(
                self._parameter_handlers.keys())
            raise KeyError(msg)

        return handler

    def get_io_handler(self, io_format):
        """Retrieve the parameter handlers associated with the provided tagname.

        If the parameter handler has not yet been instantiated, it will be created.
        If a parameter handler object already exists, it will be checked for compatibility
        and an Exception raised if it is incompatible with the provided kwargs.

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
        io_handler = None
        if io_format in self._parameter_io_handlers.keys():
            io_handler = self._parameter_io_handlers[io_format]
        elif io_format in self._parameter_io_handler_classes.keys():
            new_ph_class = self._parameter_io_handler_classes[io_format]
            io_handler = self.register_parameter_io_handler(new_ph_class)
        if io_handler is None:
            msg = "Cannot find a registered parameter IO handler for format '{}'\n".format(
                io_format)
            msg += "Registered parameter IO handlers: {}\n".format(
                self._parameter_io_handlers.keys())
            raise KeyError(msg)

        return io_handler

    # TODO: Delegate this to the XML handler

    def parse(self, sources, input_format=None):
        """Parse a SMIRNOFF force field definition.

        Parameters
        ----------
        sources : string or file-like object or open file handle or URL (or iterable of these)
            A list of files defining the SMIRNOFF force field to be loaded
            Currently, only `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_ is supported.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a ``read()`` method from which the forcefield XML data can be loaded.
            If multiple files are specified, any top-level tags that are repeated will be merged if they are compatible,
            with files appearing later in the sequence resulting in parameters that have higher precedence.
            Support for multiple files is primarily intended to allow solvent parameters to be specified by listing them last in the sequence.
        input_format : string, optional. Default = None
            The format of the input file(s). If not specified, all parsers will be tried.
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

        # Process all SMIRNOFF definition files or objects
        # QUESTION: Allow users to specify forcefield URLs so they can pull forcefield definitions from the web too?
        #trees = list()
        if input_format is None:
            # If we don't know, try all known formats
            io_formats_to_try = self._parameter_io_handler_classes.keys()
        else:
            io_formats_to_try = [input_format]

        for source in sources:
            # TODO: Add support for other sources (other serialized formats, URLs, etc)
            # Parse content depending on type
            parse_successful = False

            for parameter_io_format in io_formats_to_try:
                parameter_io_handler = self.get_io_handler(parameter_io_format)
                if type(source) == bytes:
                    parameter_io_handler.parse_string(source)
                    parse_successful = True
                elif type(source) == str:
                    #try: #TODO: What if it's the text of a XML file?
                    #    byte_source = str.encode(source)
                    #    parameter_io_handler.parse_string(byte_source)
                    #    parse_successful = True
                    #except:
                    parameter_io_handler.parse_file(source)
                    parse_successful = True
                else:
                    parameter_io_handler.parse_file(source)
                    parse_successful = True
                    # TODO: The try/except stuff here is dangerous. The io handler modifies the forcefield as it reads
                    # the file. If it crashes midway through, the forcefield will be partially modified.
                #except ParseError:
                #    pass
            if not (parse_successful):
                valid_formats = [
                    iohandler.format
                    for iohandler in self._parameter_io_handlers
                ]
                msg = "Source {} does not appear to be in a known SMIRNOFF encoding.\n".format(
                    source)
                msg += "Valid formats are: {}".format(valid_formats)
                raise Exception(msg)

    # TODO : Move to ParameterIOHandler
    def _parse_aromaticity_model(self, root):
        """Parse the aromaticity model, make sure it is supported, and make sure it does not contradict previously-specified aromaticity models.

        Parameters
        ----------
        root : etree.Element
            The document root

        """
        if not 'aromaticity_model' in root.attrib:
            self._raise_parsing_exception(
                root,
                "'aromaticity_model' attribute must be specified in top-level tag"
            )

        aromaticity_model = root.attrib['aromaticity_model']

        if aromaticity_model not in topology.ALLOWED_AROMATICITY_MODELS:
            self._raise_parsing_exception(
                root,
                "'aromaticity_model' (%s) must be one of the supported models: "
                % (aromaticity_model, topology.ALLOWED_AROMATICITY_MODELS))

        if (self._aromaticity_model is not None) and (self._aromaticity_model
                                                      != aromaticity_model):
            self._raise_parsing_exception(
                root,
                "'aromaticity_model' (%s) does not match earlier read 'aromaticity_model' (%s)"
                % (aromaticity_model, self._aromaticity_model))

    def _resolve_parameter_handler_order(self):
        """Resolve the order in which ParameterHandler objects should execute to satisfy constraints.

        Returns
        -------
        Iterable of ParameterHandlers
            The ParameterHandlers in this ForceField, in the order that they should be called to satisfy constraints.
        """
        ordered_parameter_handlers = list()
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
                    self.get_handler(tagname, {}))
        #ordered_parameter_handlers = [ self.get_handler(tagname, {}) for tagname in nx.topological_sort(G) ]
        return ordered_parameter_handlers

    # TODO: Should we add convenience methods to parameterize a Topology and export directly to AMBER, gromacs, CHARMM, etc.?
    #       Or should we create an "enhanced" openforcefield System object that knows how to convert to all of these formats?
    #       We could even create a universal applyParameters(format='AMBER') method that allows us to export to whatever system we want.

    # TODO: Should the Topology contain the default box vectors? Or should we require they be specified externally?

    # TODO: How do we know if the system is periodic or not?
    # TODO: Should we also accept a Molecule as an alternative to a Topology?

    def create_openmm_system(self,
                             topology,
                             default_box_vectors=None,
                             **kwargs):
        """Create an OpenMM System representing the interactions for the specified Topology with the current force field

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The ``Topology`` corresponding to the system to be parameterized
        default_box_vectors : simtk.unit.Quanity of shape [3,3] with units compatible with nanometers, optional, default=None
            Default box vectors to use.
            If not specified, default box vectors will be set to 1.0 nm edges.
            Note that, for periodic systems, after creating a Context, box vectors *must* be set to the appropriate dimensions.
        verbose : bool
            If True, verbose output will be printed.

        Returns
        -------
        system : simtk.openmm.System
            The newly created OpenMM System corresponding to the specified ``topology``

        """
        # Make a deep copy of the topology so we don't accidentally modify it
        topology = copy.deepcopy(topology)

        # Set the topology aromaticity model to that used by the current forcefield
        #topology.set_aromaticity_model(self._aromaticity_model)

        # Create an empty OpenMM System
        system = openmm.System()

        # Set periodic boundary conditions if specified
        if default_box_vectors is not None:
            system.setDefaultPeriodicBoxVectors(default_box_vectors)

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
        # TODO: Delete kwargs?
        for parameter_handler in parameter_handlers:
            parameter_handler.create_force(system, topology, **kwargs)

        # Let force Handlers do postprocessing
        # TODO: Delete kwargs?
        for parameter_handler in parameter_handlers:
            parameter_handler.postprocess_system(system, topology, **kwargs)

        return system

    def create_parmed_structure(self,
                                topology,
                                positions,
                                default_box_vectors=None,
                                **kwargs):
        """Create a ParmEd Structure object representing the interactions for the specified Topology with the current force field

        This method creates a `ParmEd <http://github.com/parmed/parmed>`_ ``Structure`` object containing a topology, positions, and parameters.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The ``Topology`` corresponding to the ``System`` object to be created.
        positions : simtk.unit.Quantity of dimension (natoms,3) with units compatible with angstroms
            The positions corresponding to the ``System`` object to be created
        default_box_vectors : simtk.unit.Quanity of shape [3,3] with units compatible with nanometers, optional, default=None
            Default box vectors to use.
            If not specified, default box vectors will be set to 1.0 nm edges.
            Note that, for periodic systems, after creating a Context, box vectors *must* be set to the appropriate dimensions.
        verbose : bool
            If True, verbose output will be printed.

        Returns
        -------
        structure : parmed.Structure
            The newly created ``parmed.Structure`` object

        """
        raise NotImplementedError
        import parmed
        # TODO: Automagically handle expansion of virtual sites? Or is Topology supposed to do that?

        # Create OpenMM System
        system = self.create_openmm_system(
            topology, default_box_vectors=default_box_vectors, **kwargs)

        # Create a ParmEd Structure object
        structure = parmed.openmm.topsystem.load_topology(
            topology.to_openmm(), system, positions)

        return structure

    def label_molecules(self, topology, verbose=False):
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
            ``[ ( [ atom1, ..., atomN], parameter_id, SMIRKS), ... ]``

        .. todo ::

           What is the most useful API for this method?
           Should we instead accept :class:`Molecule` objects as input and individually return labels?
           Should we attach the labels to the :class:`Molecule` object?
           Or should we label all interactions in a :class:`Topology` instead of just labeling its ``unique_molecules``?

        """
        # Loop over molecules and label
        molecule_labels = list()
        for molecule in enumerate(topology.unique_molecules):
            current_molecule_labels = dict()
            for force in self.forces:
                matches = force.get_matches(molecule)
                molecule_labels[idx][force.name] = matches
            molecule_labels.append(current_molecule_labels)
        return molecule_labels
