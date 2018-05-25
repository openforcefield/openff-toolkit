#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Parser for the SMIRNOFF (SMIRKS Native Open Force Field) format.

.. codeauthor:: John D. Chodera <john.chodera@choderalab.org>
.. codeauthor:: David L. Mobley <dmobley@mobleylab.org>
.. codeauthor:: Peter K. Eastman <peastman@stanford.edu>

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import re
import sys
import math
import copy
import time
import string
import random
import logger
import itertools
import collections

import numpy
import packaging

import lxml.etree as etree

from simtk import openmm, unit
from simtk.openmm.app import element as elem

from openforcefield.utils import get_data_filename, all_subclasses
from openforcefield.topology import Topology, ValenceDict, ImproperDict
from openforcefield.topology import DEFAULT_AROMATICITY_MODEL
from openforcefield.typing.chemistry import ChemicalEnvironment, SMIRKSParsingError

#=============================================================================================
# CONFIGURE LOGGER
#=============================================================================================

logger = logging.getLogger(__name__)

#=============================================================================================
# PRIVATE METHODS
#=============================================================================================

MAX_SUPPORTED_VERION = '1.0' # maximum version of the SMIRNOFF spec supported by this SMIRNOFF forcefield

#=============================================================================================
# FORCEFIELD
#=============================================================================================

# QUESTION: Should we process the XML files only when the ForceField is created, or should we be able to read more XML files later?

# QUESTION: How should we document private object fields?

# TODO: How do we serialize/deserialize this object?

# TODO: Should we add methods to retrieve string XML representations?

# QUESTION: How should we support other non-XML representations in future? Should we integrate these representations throughout, or collect all the XML reading/writing in one place?

class ForceField(object):
    """A factory initialized from a SMIRNOFF force field definition that assigns parameters to a molecular system.

    Specifically, :class:`ForceField` is a factory that constructs an OpenMM :class:`simtk.openmm.System` object from a :class:`openforcefield.topology.Topology` object.

    When a :class:`ForceField` object is created from one or more specified XML files, all :class:`ParameterHandler` subclasses
    currently imported are identified and registered to handle different sections of the SMIRNOFF force field definition file(s).
    The force field definition is processed by these Handlers to populate internal parameter definition data structures that
    can be manipulated via the API.

    Processing a :class:`Topology` object defining a chemical system will then call all :class`ParameterHandler` objects in an order
    guaranteed to satisfy the declared processing order constraints of each :class`ParameterHandler`.

    Attributes
    ----------
    parsers : dict of str : openforcefield.typing.engines.smirnoff.ParameterHandler
        Registered list of parsers that will handle forcefield tags.
        parsers[tagname] is the ``ParameterHandler`` that will be called to process the force field definition section ``tagname``
        Parsers are registered when the ForceField object is created, but can be manipulated afterwards.

    """

    def __init__(self, *sources, parameter_handlers=None, disable_version_check=False):
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
        parameter_handlers : iterable of ParameterHandler classes, optional, default=None
            If not None, the specified set of ParameterHandler objects will be used to parse the XML tags.
            By default, all imported subclasses of ParameterHandler are automatically registered to parse XML tags.
        disable_version_check : bool, optional, default=False
            If True, will disable checks against the current highest supported forcefield version.
            This option is primarily intended for forcefield development.

        """
        # Clear all object fields
        self._initialize()

        # Store initialization options
        self.disable_version_check = disable_version_check # if True, we won't check which SMIRNOFF version number we're parsing

        # Register all ParameterHandler objects that will process SMIRNOFF force definitions
        # TODO: Should force Handlers be specified as classes or objects?
        # TODO: Should we call them something besides parameter_handlers?
        if parameter_handlers is None:
            # Find all imported subclasses of ParameterHandler
            parameter_handlers = self._find_parameter_handlers()
        self._register_parsers(parameter_handlers)

        # Parse all sources containing SMIRNOFF parameter definitions
        self.parse(sources)

    def _initialize(self):
        """Initialize all object fields.
        """
        self.disable_version_check = False # if True, will disable checking compatibility version
        self.aromaticity_model = DEFAULT_AROMATICITY_MODEL # aromaticity model
        self.parameter_handlers = OrderedDict() # ParameterHandler classes to handle parameter types
        self.parameters = OrderedDict() # ParameterHandler objects instantiated fro parameter type

    # TODO: We need to change this to just find all ParameterHandler objects in this file;
    # otherwise, we can't define two different ParameterHandler subclasses to compare for a new type of energy term
    # since both will try to register themselves for the same XML tag and an Exception will be raised.
    @staticmethod
    def _find_parameter_handlers():
        """Identify all imported subclasses of ParameterHandler.

        Returns
        -------
        parameter_handlers : list of ParameterHandler subclasses
            List of ParameterHandler subclasses (not objects)
        """
        return all_subclasses(ParameterHandler)

    @staticmethod
    def _register_parsers(parameter_handlers):
        """Register all force Handlers with their ``_TAGNAME`` tag namesself.

        Parameters
        ----------
        parameter_handlers : list of ParameterHandler subclasses
            The force Handlers to register.

        """
        parsers = dict()
        for parameter_handler in parameter_handlers:
            tagname = subclass._TAGNAME
            if tagname is not None:
                if tagname in parsers.keys():
                    raise Exception("ParameterHandler {} provides a parser for tag '{}', but ParameterHandler {} has already been registered to handle that tag.".format(subclass, tagname, self.parsers[tagname]))
                parsers[tagname] = subclass
        return parsers

    @staticmethod
    def _raise_parsing_exception(node=None, msg=None):
        """Raise a ValueError during XML file parsing, printing the source line number to aid debugging if the XML parsing library supports itself.

        Parameters
        ----------
        node : lxml.etree.Node like, optional, default=None
            XML DOM node object, which may optionally have a ``sourceline`` field
        msg : str, optional, default=None
            Exception message to include in ValueError

        """
        if (node is not None):
            if hasattr(node, 'sourceline'):
                raise ValueError("Line %s : %s\n%s" % (node.sourceline, str(node), msg))
            else:
                raise ValueError("%s\n%s" % (str(node), msg))
        else:
            raise Exception(msg)

    @staticmethod
    def _validate_smarts(smarts, node=None, ensure_valence_type=None):
        """Validate the specified SMARTS string can be used to assign forcefield parameters.

        This checks to ensure the SMARTS string
        * is a valid SMARTS string
        * the tagged atoms form a fully connected subset of atoms
        * if ``ensure_valence_type`` is specified, ensure the tagged atoms specify the appropriate valence type

        Parameters
        ----------
        smarts : str
           The SMARTS string to be validated
        node : xml.etree.ElementTree.Element, optional, default=None
           Node of etree, used only for reporting errors
        ensure_valence_type : str, optional, default=None
           If not ``None``, will check to ensure tagged atoms specify appropriate valence types
           Supported ChemicalEnvironment getType() types: ['Atom', 'Bond', 'Angle', 'ProperTorsion', 'ImproperTorsion']
           If ``None``, will ensure that it is any one of the above valid valence types.

        """
        # Create a chemical environment to see if this is a valid SMARTS string
        try:
            chemenv = ChemicalEnvironment(smarts)
        except Exception as e:
            self._raise_parsing_exception(node, "Error parsing SMARTS '%s' : %s" % (smarts, str(e)))

        # Check type, if specified
        ensure_valence_type = chemenv.getType()
        if ensure_valence_type:
            if valence_type != ensure_valence_type:
                self._raise_parsing_exception(node, "Tagged atoms in SMARTS string '%s' specifies valence type '%s', expected '%s'." % (smarts, valence_type, ensure_valence_type))
        else:
            if valence_type is None:
                self._raise_parsing_exception(node, "Tagged atoms in SMARTS string '%s' did not tag atoms in a way that correctly specifies a valence type." % smarts)

    @staticmethod
    def _extract_quantity_from_xml_element(node, parent, name, unit_name=None, default=None):
        """
        Form a (potentially unit-bearing) quantity from the specified attribute name.

        node : xml.etree.ElementTree.Element
           Node of etree corresponding to force type entry.
        parent : xml.etree.ElementTree.Element
           Node of etree corresponding to parent Force.
        name : str
           Name of parameter to extract from attributes.
        unit_name : str, optional, default=None
           If specified, use this attribute name of 'parent' to look up units
        default : optional, default=None
           If not None, the value in ``default`` will be returned if ``name`` is not found
           instead of raising an exception.

        """
        # Check for expected attributes
        if (name not in node.attrib):
            self._raise_parsing_exception(node, "Expected XML attribute '%s' not found" % (name))

        # Most attributes will be converted to floats, but some are strings
        string_names = ['parent_id', 'id']
        # Handle case where this is a normal quantity
        if name not in string_names:
            quantity = float(node.attrib[name])
        # Handle label or string
        else:
            quantity = node.attrib[name]
            return quantity

        if unit_name is None:
            unit_name = name + '_unit'

        if unit_name in parent.attrib:
            # TODO: This is very dangerous. Replace it with safer scheme from YANK.
            string = '(%s * %s).value_in_unit_system(md_unit_system)' % (node.attrib[name], parent.attrib[unit_name])
            quantity = eval(string, unit.__dict__)

        return quantity

    @staticmethod
    def _check_for_missing_valence_terms(name, topology, assigned_terms, topological_terms):
        """
        Check to ensure there are no missing valence terms in the given topology, identifying potential gaps in parameter coverage.

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
        assigned_terms = [ item for item in assigned_terms ]
        topological_terms = [ item for item in topological_terms ]

        def ordered_tuple(atoms):
            atoms = list(atoms)
            if atoms[0] < atoms[-1]:
                return tuple(atoms)
            else:
                return tuple(reversed(atoms))
        try:
            topology_set = set([ ordered_tuple( atom.index for atom in atomset ) for atomset in topological_terms ])
            assigned_set = set([ ordered_tuple( index for index in atomset ) for atomset in assigned_terms ])
        except TypeError as te:
            topology_set = set([ atom.index for atom in topological_terms ])
            assigned_set = set([ atomset[0] for atomset in assigned_terms ])

        def render_atoms(atomsets):
            msg = ""
            for atomset in atomsets:
                msg += '%30s :' % str(atomset)
                try:
                    for atom_index in atomset:
                        atom = atoms[atom_index]
                        msg += ' %5s %3s %3s' % (atom.residue.index, atom.residue.name, atom.name)
                except TypeError as te:
                    atom = atoms[atomset]
                    msg += ' %5s %3s %3s' % (atom.residue.index, atom.residue.name, atom.name)

                msg += '\n'
            return msg

        if set(assigned_set) != set(topology_set):
            # Form informative error message
            msg = '%s: Mismatch between valence terms added and topological terms expected.\n' % name
            atoms = [ atom for atom in topology.atoms ]
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
            raise Exception(msg) # TODO: Should we raise a more specific exception here?

    def parse(self, sources):
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

        .. notes ::

           * New SMIRNOFF sections are handled independently, as if they were specified in the same file.
           * If a SMIRNOFF section that has already been read appears again, its definitions are appended to the end of the previously-read
             definitions if the sections are configured with compatible attributes; otherwise, an ``IncompatibleTagException`` is raised.

        """
        # Ensure that we are working with an iterable
        try:
            some_object_iterator = iter(files)
        except TypeError as te:
            # Make iterable object
            files = [files]

        # Process all SMIRNOFF definition files or objects
        # QUESTION: Allow users to specify forcefield URLs so they can pull forcefield definitions from the web too?
        trees = list()
        for source in sources:
            # TODO: Load content from source

            # Parse content depending on type
            try:
                self.parse_xml(source_string)
            except ParseError:
                msg = "Source {} does not appear to be in a known SMIRNOFF encoding.\n".format(source)
                msg += "Valid encodings are: ['XML']"
                raise Exception(msg)

    def parse_xml(self, xml):
        """Parse a SMIRNOFF force field definition in XML format.

        A ``ParseError`` is raised if the XML cannot be processed.

        Parameters
        ----------
        xml : string
            A SMIRNOFF force field definition in `the SMIRNOFF XML format <https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md>`_.

        .. notes ::

           * New SMIRNOFF sections are handled independently, as if they were specified in the same file.
           * If a SMIRNOFF section that has already been read appears again, its definitions are appended to the end of the previously-read
             definitions if the sections are configured with compatible attributes; otherwise, an ``IncompatibleTagException`` is raised.
        """

        parser = etree.XMLParser(remove_blank_text=True) # For pretty print on write
        try:
            # this handles either filenames or open file-like objects
            tree = etree.parse(file, parser)
        except IOError:
            temp_file = get_data_filename(file)
            tree = etree.parse(temp_file, parser)
        except Exception as e:
            # Fail with an error message about which file could not be read.
            # TODO: Also handle case where fallback to 'data' directory encounters problems,
            # but this is much less worrisome because we control those files.
            msg  = str(e) + '\n'
            if hasattr(file, 'name'):
                filename = file.name
            else:
                filename = str(file)
            msg += "ForceField.loadFile() encountered an error reading file '%s'\n" % filename
            raise Exception(msg)

        # Run through parsers

    def save(self, filename):
        """Write the current forcefield parameter set to a file, autodetecting the type from the extension.

        Parameters
        ----------
        filename : str
            The name of the file to be written.
            The `.offxml` file extension is auto-detected.

        """
        (basename, extension) = os.path.splitext(filename)
        if extension == '.offxml':
            xml = self.to_xml()
            with open(filename, 'w') as f:
                f.write(xml)
        else:
            msg = 'Export to extension {} not implemented yet.\n'.format(extension)
            msg += "Supported choices are: ['.offxml']"
            raise NotImplementedError(msg)

    @staticmethod
    def from_xml(xml):
        """Construct a ForceField object from a SMIRNOFF XML file.
        """
        return ForceField(xml)

    def to_lxml(self):
        """Render the forcefield parameter set to an lxml.etree

        Returns
        -------
        root : lxml.etree
            Root node
        """
        # Create the XML tree
        root = etree.XML('<SMIRNOFF/>')
        # Populate top-level attributes
        root['version'] = MAX_SUPPORTED_VERSION # TODO: Should we instead store a version?
        root['aromaticity_mode'] = self.aromaticity_model
        # Render all force Handlers
        for tag, parameter_handler in self.parameter_handlers.items():
            subtree = parameter_handler.to_lxml()
            root.insert(len(root), subtree)

        return root

    def to_xml(self):
        """Render the forcefield parameter set to XML.

        Returns
        -------
        xml : str
            The SMIRNOFF parameter set rendered as XML.
        """
        return etree.tostring(self.to_lxml())

    # TODO: Do we need this? Should we use built-in dict-based serialization?
    def __getstate__(self):
        """Serialize to XML.
        """
        return self.to_xml()

    # TODO: Do we need this? Should we use built-in dict-based serialization?
    def __setstate__(self, state):
        """Deserialize from XML.
        """
        self._initialize()
        self.parse_xml(state)

    def _parse_version(self, root):
        """Parse the forcefield version number and make sure it is supported.

        Parameters
        ----------
        root : etree.Element
            The document root

        """
        if 'version' in root.attrib:
            version = root.attrib['version']
            # Use PEP-440 compliant version number comparison, if requested
            if (not self.disable_version_check) and (packaging.version.parse(version) > packaging.version.parse(MAX_SUPPORTED_VERION)):
                self._raise_parsing_exception(root, 'SMIRNOFF offxml file was written with version %s, but this version of ForceField only supports up to version %s' % (self.version, MAX_SUPPORTED_VERSION))
        else:
            self._raise_parsing_exception(root, "'version' attribute must be specified in SMIRNOFF tag")

    def _parse_aromaticity_model(self, root):
        """Parse the aromaticity model, make sure it is supported, and make sure it does not contradict previously-specified aromaticity models.

        Parameters
        ----------
        root : etree.Element
            The document root

        """
        if not 'aromaticity_model' in root.attrib:
            self._raise_parsing_exception(root, "'aromaticity_model' attribute must be specified in top-level tag")

        aromaticity_model = root.attrib['aromaticity_model']

        if aromaticity_model not in topology.ALLOWED_AROMATICITY_MODELS:
            self._raise_parsing_exception(root, "'aromaticity_model' (%s) must be one of the supported models: " % (aromaticity_model, topology.ALLOWED_AROMATICITY_MODELS))

        if (self._aromaticity_model is not None) and (self._aromaticity_model != aromaticity_model):
            self._raise_parsing_exception(root, "'aromaticity_model' (%s) does not match earlier read 'aromaticity_model' (%s)" % (aromaticity_model, self._aromaticity_model))

    def parseXMLTrees(self):
        """
        Initialize all registered ParameterHandler objects by parsing loaded XML files.

        .. notes ::

           * New SMIRNOFF tags are handled independently, as if they were specified in the same file.
           * If a SMIRNOFF tag that has already been read appears again, its definitions are appended to the end of the previously-read
             definitions if the tags are configured with compatible attributes; otherwise, an ``IncompatibleTagException`` is raised.

        """
        trees = self._XMLTrees

        # Create all ParameterHandler objects from scratch
        self._forces = list()

        # Load force definitions
        for tree in trees:
            root = tree.getroot()

            # Formerly known as SMIRFF, now SMIRNOFF
            if not root.tag in ['SMIRFF',  'SMIRNOFF']:
                self._raise_parsing_exception(root, "Error: ForceField parses a SMIRNOFF forcefield, but this does not appear to be one as the root tag is %s." % root.tag)

            # Parse attributes from top-level tag
            self._parse_version(root)
            self._parse_aromaticity_model(root)

            # Load force parameters using registered parsers
            for child in root:
                if child.tag in parsers:
                    self.parsers[child.tag].parseElement(child.tag, child, self)
                else:
                    msg = "There is no registered parser for the tag '{}'. Parsers are registered for the following tags: {}".format(tagname, self.parsers.keys())
                    raise Exception(msg)

        # Reset flag tracking whether stored XML has been modified since the last parse
        self._XMLModified = False

    def _reparse_xml_if_needed(self):
        """
        Ensure the XML trees have been parsed if they have been changed.

        """
        if self._XMLModified:
            if verbose: print("Re-parsing XML because it was modified.")
            self.parseXMLTrees()

    def getHandlers(self):
        """Get the list of all registered Handlers."""
        return self._forces

    def registerHandler(self, Handler):
        """Register a new Handler."""
        # Special case: ConstraintHandler has to come before BondHandler and AngleHandler.
        # TODO: Figure out a more general way to allow Handlers to specify enforced orderings via dependencies.
        if isinstance(Handler, ConstraintHandler):
            self._forces.insert(0, Handler)
        else:
            self._forces.append(Handler)

    # TODO: Rework the API for writing ffxml parameters, since writing to multiple files is currently ambiguous
    def writeFile(self, files):
        """Write forcefield trees out to specified files.

        Parameters
        ----------
        files : str or tuple of str
            File or files to write XML trees to.

        """

        # Ensure that we are working with a tuple of files.
        if not isinstance(files, tuple):
            files = (files,)

        for idx, filenm in enumerate(files):
            tree = self._XMLTrees[idx]
            tree.write(filenm, xml_declaration=True, pretty_print=True)

    def _resolve_parameter_handler_order(self):
        """Resolve the order in which ParameterHandler objects should execute to satisfy constraints.
        """
        ordered_parameter_handlers = list()
        # Create a DAG expressing dependencies
        import networkx as nx
        G = nx.DiGraph()
        for parameter_handler in self.parsers.items():
            G.add_node(parameter_handler._TAGNAME)
            if parameter_handler._DEPENDENCIES is not None:
                for dependency in parameter_handler._DEPENDENCIES:
                    G.add_edge(dependency._TAGNAME, parameter_handler._TAGNAME)
        # TODO: Check to make sure DAG isn't cyclic
        # Resolve order
        ordered_parameter_handlers = [ self.parsers[tagname] for tagname in nx.topological_sort(G) ]
        return ordered_parameter_handlers

    def createSystem(self, topology, box_vectors=None, verbose=False, **kwargs):
        """Construct an OpenMM System representing a Topology with this force field. XML will be re-parsed if it is modified prior to system creation.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The ``Topology`` corresponding to the ``System`` object to be created.
        default_box_vectors : simtk.unit.Quanity of shape [3,3] with units compatible with nanometers, optional, default=None
            Default box vectors to use.
            If not specified, default box vectors will be set to 1.0 nm edges.
            Note that, for periodic systems, after creating a Context, box vectors *must* be set to the appropriate dimensions.
        verbose : bool
            If True, verbose output will be printed.
        kwargs
            Arbitrary additional keyword arguments may also be specified.
            This allows extra parameters to be specified that are specific to particular force fields.
            A ``ValueError`` will be raised if parameters are specified that are not used by any force Handlers.

        Returns
        -------
        system : simtk.openmm.System
            The newly created System

        """
        # TODO: Have `verbose` flag set whether logging info is displayed or not.

        # Re-parse XML by ParameterHandler objects if loaded XML files have been added or modified
        self._reparse_xml_if_needed()

        # Make a deep copy of the topology so we don't accidentaly modify it
        topology = copy.deepcopy(topology)

        # Set aromaticity model to that used by forcefield
        topology.set_aromaticity_model(self._aromaticity_model)

        # Create the System
        system = openmm.System()

        # Add particles
        # TODO: Optionally allow SMARTS-specified masses
        for atom in topology.atoms:
            system.addParticle(atom.element.mass)

        # Set periodic boundary conditions if specified
        if default_box_vectors is not None:
            system.setDefaultPeriodicBoxVectors(default_box_vectors)

        # Determine the order in which to process ParameterHandler objects in order to satisfy dependencies
        parameter_handlers = self._resolve_parameter_handler_order()

        # Check if any kwargs have been provided that aren't handled by force Handlers
        known_kwargs = set()
        for parameter_handler in parameter_handlers:
            known_args.update(parameter_handler.known_kwargs)
        unknown_kwargs = set(kwargs.keys()).difference(known_kwargs)
        if len(unknown_kwargs) > 0:
            msg = "The following keyword arguments to createSystem() are not used by any registered force Handler: {}\n".format(unknown_kwargs)
            msg += "Known keyword arguments: {}".format(known_kwargs)
            raise ValueError(msg)

        # Add forces to the System
        for parameter_handler in parameter_handlers:
            parameter_handler.createForce(system, topology, **kwargs)

        # Let force Handlers do postprocessing
        for parameter_handler in parameter_handlers:
            parameter_handler.postprocessSystem(system, topology, **kwargs)

        return system

    # TODO: Rework this
    def labelMolecules(self, molecules, verbose=False):
        """Return labels for a list of molecules corresponding to parameters from this force field.
        For each molecule, a dictionary of force types is returned, and for each force type,
        each force term is provided with the atoms involved, the parameter id assigned, and the corresponding SMIRKS.

        Parameters
        ----------
        molecules : list of openforcefield.topology.Molecule
            The molecules to be labeled.

        Returns
        -------
        molecule_labels : list
            List of labels for molecules. Each entry in the list corresponds to
            one molecule from the provided list of molecules and is a dictionary
            keyed by force type, i.e., ``molecule_labels[0]['HarmonicBondForce']``
            gives details for the harmonic bond parameters for the first
            molecule. Each element is a list of the form:
            ``[ ( [ atom1, ..., atomN], parameter_id, SMIRKS), ... ]``

        """
        self._ensure_xml_has_been_parsed()

        # Loop over molecules and label
        molecule_labels = list()
        for (idx, mol) in enumerate(molecules):
            molecule_labels.append(dict())
            for force in self._forces:
                matches = force.getMatches(mol)
                # TODO Reformat matches
                # QUESTION: Do we want to store labeled force terms by OpenMM force type (e.g., `HarmonicBondForce`) or by SMIRNOFF tag name (`BondForce`)
                molecule_labels[idx][force._TAGNAME] = matches

        return molecule_labels

#=============================================================================================
# The following classes are Handlers that know how to create Force subclasses and add them to a System that is being
# created.  Each Handler class must define three methods: 1) a static method that takes an etree Element and a ForceField,
# and returns the corresponding Handler object; 2) a createForce() method that constructs the Force object and adds it
# to the System; and 3) a labelForce() method that provides access to which
# terms are applied to which atoms in specified oemols.
# The static method should be added to the parsers map.
#=============================================================================================

#=============================================================================================
# Force Handlers
#=============================================================================================

## @private
class ParameterHandler(object):
    """Base class for parameter type handlers.
    """
    _TAGNAME = None # str of section type handled by this ParameterHandler (XML element name for SMIRNOFF XML representation)
    _VALENCE_TYPE = None # ChemicalEnvironment valence type string expected by SMARTS string for this Handler
    _INFOTYPE = None # container class with type information that will be stored in self._types
    _OPENMMTYPE = None # OpenMM Force class (or None if no equivalent)
    _DEPENDENCIES = None # list of ParameterHandler classes that must precede this, or None
    _DEFAULTS = {} # dict of attributes and their default values at tag-level
    _KWARGS = [] # list of keyword arguments accepted by the force Handler on initialization
    _SMIRNOFF_VERSION_INTRODUCED = 0.0 # the earliest version of SMIRNOFF spec that supports this ParameterHandler
    _SMIRNOFF_VERSION_DEPRECATED = None # if deprecated, the first SMIRNOFF version number it is no longer used
    _REQUIRE_UNITS = None # list of parameters that require units to be defined

    def __init__(self, forcefield):
        self._forcefield = forcefield # the ForceField object that this ParameterHandler is registered with
        self.parameters = list() # list of ParmaeterType objects

    @property
    def known_kwargs(self):
        """List of kwargs that can be parsed by the function.
        """
        # TODO: Should we use introspection to inspect the function signature instead?
        return set(self._KWARGS)

    def getMatches(self, topology):
        """Retrieve all force terms for a chemical entity.

        Parameters
        ----------
        entity : openforcefield.topology.ChemicalEntity
            Chemical entity for which constraints are to be enumerated

        Returns
        ---------
        matches : ValenceDict
            matches[atoms] is the ParameterType object corresponding to the tuple of Atom objects ``Atoms``

        """
        logger.info(self.__class__.__name__)
        matches = ValenceDict()
        for force_type in self._types:
            matches_for_this_type = { atoms : force_type for atoms in topology.chemical_environment_matches(force_type.smirks) }
            matches.update(matches_for_this_type)
            logger.info('{:64} : {:8} matches'.format(force_type.smirks, len(matches_for_this_type)))

        logger.info('{} matches identified'.format(len(matches)))
        return matches

    # TODO: Handle appending tags
    @classmethod
    def parseElement(cls, tag, element, ff):
        """
        Parse the XML tag/section this ParameterHandler is registered for.

        SMIRNOFF sections may be split across multiple files or otherwise appear multiple times,
        so we need to be able to handle multiple calls to parseElement().

        Parameters
        ----------
        tag : str
        element : lxml.Element
        ff : openforcefield.typing.engines.smirnoff.ForceField
            The ForceField object parsing the element

        """
        existing = [f for f in ff._forces if isinstance(f, cls)]
        if len(existing) == 0:
            Handler = cls(ff)
            ff.registerHandler(Handler)
        else:
            Handler = existing[0]

        # Extract all tag-level attributes (or their defaults)
        for attribute in _DEFAULTS.keys():
            setattr(self, attribute, _extract_quantity_from_xml_element(node, parent, attribute, default=DEFAULTS[attribute]))

        # Register all SMIRNOFF definitions for all occurrences of its registered tag
        for section in element.findall(tag):
            Handler.registerType(section, element)

        # TODO: Check that all required units are defined in the top-level tag attributes

    def registerType(self, node, parent):
        """Register a SMIRNOFF constraint type definition."""
        if self._INFOTYPE:
            self._types.append(self._INFOTYPE(node, parent))

    def createForce(self, system, topology, **kwargs):
        force = None
        if self._OPENMMTYPE:
            # Find existing force or create new one.
            existing = [ system.getForce(i) for i in range(system.getNumForces()) ]
            existing = [ f for f in existing if type(f) == openmm.HarmonicBondForce ]
            if len(existing) == 0:
                force = _OPENMMTYPE()
                system.addForce(force)
            else:
                force = existing[0]

        return force

    def postprocessSystem(self, system, topology, **args):
        pass

# TODO: Rename to better reflect role as parameter base class?
class ParameterType(object):
    """
    Base class for SMIRNOFF parameter types.

    """
    def __init__(self, smirks, **kwargs):
        """
        Create a ParameterType

        Parameters
        ----------

        """
        # Store SMIRKS
        self.smirks = smirks # TODO: Validate SMIRKS is appropriate for this type?
        # Store any other attributes provided
        for k, v in kwargs.items():
            setattr(self, k, v)

    @staticmethod
    def from_xml(node, parent):
        smirks = _validate_smarts(node.attrib['smirks'], node=node, valence_type=_VALENCE_TYPE)
        # TODO: Generally process other node attributes?
        kwargs = dict()
        if 'id' in node.attrib:
            kwargs['id'] = _extract_quantity_from_xml_element(node, parent, 'id')
        return ParameterType(smirks, **kwargs)

#=============================================================================================

## @private
class ConstraintHandler(ParameterHandler):
    """Handle SMIRNOFF ``<Constraints>`` tags

    ``ConstraintHandler`` must be applied before ``BondHandler`` and ``AngleHandler``,
    since those classes add constraints for which equilibrium geometries are needed from those tags.
    """
    class ConstraintType(ParameterType):
        """A SMIRNOFF constraint type"""
        def __init__(self, node, parent):
            super(ConstraintType, self).__init__(node, parent) # Base class handles ``smirks`` and ``id`` fields
            if 'distance' in node.attrib:
                self.distance = _extract_quantity_from_xml_element(node, parent, 'distance') # Constraint with specified distance will be added by ConstraintHandler
            else:
                self.distance = True # Constraint to equilibrium bond length will be added by HarmonicBondHandler

    _TAGNAME = 'Constraint'
    _VALENCE_TYPE = 'Bond' # ChemicalEnvironment valence type expected for SMARTS # TODO: Do we support more exotic types as well?
    _INFOTYPE = ConstraintType
    _OPENMMTYPE = None # don't create a corresponding OpenMM Force class
    _REQUIRED_UNITS = ['distance']

    def __init__(self, forcefield):
        super(ConstraintHandler, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        constraints = self.getMatches(topology)
        for (atoms, constraint) in constraints.items():
            # Update constrained atom pairs in topology
            topology.add_constraint(*atoms, constraint.distance)
            # If a distance is specified (constraint.distance != True), add the constraint here.
            # Otherwise, the equilibrium bond length will be used to constrain the atoms in HarmonicBondHandler
            if constraint.distance is not True:
                system.addConstraint(*atoms, constraint.distance)

#=============================================================================================

## @private
class BondHandler(ParameterHandler):
    """Handle SMIRNOFF ``<BondForce>`` tags"""

    class BondType(ParameterType):
        """A SMIRNOFF Bond parameter type"""
        def __init__(self, node, parent):
            super(ConstraintType, self).__init__(node, parent) # Base class handles ``smirks`` and ``id`` fields

            # Determine if we are using fractional bond orders for this bond
            # First, check if this force uses fractional bond orders
            if 'fractional_bondorder_method' in parent.attrib:
                # If it does, see if this parameter line provides fractional bond order parameters
                if 'length_bondorder1' in node.attrib and 'k_bondorder1' in node.attrib:
                    # Store what interpolation scheme we're using
                    self.fractional_bondorder = parent.attrib['fractional_bondorder']
                    # Store bondorder1 and bondorder2 parameters
                    self.k = list()
                    self.length = list()
                    for ct in range(1,3):
                        self.length.append( _extract_quantity_from_xml_element(node, parent, 'length_bondorder%s' % ct, unit_name = 'length_unit') )
                        self.k.append( _extract_quantity_from_xml_element(node, parent, 'k_bondorder%s' % ct, unit_name = 'k_unit') )
                else:
                    self.fractional_bondorder = None
            else:
                self.fractional_bondorder = None

            # If no fractional bond orders, just get normal length and k
            if self.fractional_bondorder is None:
                self.length = _extract_quantity_from_xml_element(node, parent, 'length')
                self.k = _extract_quantity_from_xml_element(node, parent, 'k')

    _TAGNAME = 'BondForce' # SMIRNOFF tag name to process
    _VALENCE_TYPE = 'Bond' # ChemicalEnvironment valence type expected for SMARTS
    _INFOTYPE = BondType # class to hold force type info
    _OPENMMTYPE = openmm.HarmonicBondForce # OpenMM force class to create
    _DEPENDENCIES = [ConstraintHandler] # ConstraintHandler must be executed first

    def __init__(self, forcefield):
        super(HarmonicBondHandler, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        # Create or retrieve existing OpenMM Force object
        force = super(BondHandler, self).createForce(system, topology, **kwargs)

        # Add all bonds to the system.
        bonds = self.getMatches(topology)
        skipped_constrained_bonds = 0 # keep track of how many bonds were constrained (and hence skipped)
        for (atoms, bond) in bonds.items():
            # Get corresponding particle indices in Topology
            particle_indices = tuple([ atom.particle_index for atom in atoms ])

            # Ensure atoms are actually bonded correct pattern in Topology
            topology.assert_bonded(atoms[0], atoms[1])

            # Compute equilibrium bond length and spring constant.
            if bond.fractional_bondorder is None:
                [k, length] = [bond.k, bond.length]
            else:
                # Interpolate using fractional bond orders
                # TODO: Do we really want to allow per-bond specification of interpolation schemes?
                order = topology.get_fractional_bond_order(*atoms)
                if bond.fractional_bondorder_interpolation == 'interpolate-linear':
                    k = bond.k[0] + (bond.k[1]-bond.k[0])*(order-1.)
                    length = bond.length[0] + (bond.length[1]-bond.length[0])*(order-1.)
                else:
                    raise Exception("Partial bondorder treatment {} is not implemented.".format(bond.fractional_bondorder))

            # Handle constraints.
            if topology.atom_pair_is_constrained(*atoms):
                # Atom pair is constrained; we don't need to add a bond term.
                skipped_constrained_bonds += 1
                # Check if we need to add the constraint here to the equilibrium bond length.
                if topology.atom_pair_is_constrained(*atoms) is True:
                    # Mark that we have now assigned a specific constraint distance to this constraint.
                    topology.add_constraint(*atoms, length)
                    # Add the constraint to the System.
                    system.addConstraint(*particle_indices, length)
                continue

            # Add harmonic bond to HarmonicBondForce
            force.addBond(*particle_indices, length, k)

        logger.info('{} bonds added ({} skipped due to constraints)'.format(len(bonds) - skipped_constrained_bonds, skipped_constrained_bonds))

        # Check that no topological bonds are missing force parameters
        _check_for_missing_valence_terms('BondForce', topology, bonds.keys(), topology.bonds)

#=============================================================================================

## @private
class AngleHandler(ParameterHandler):
    """Handle SMIRNOFF ``<AngleForce>`` tags"""

    class AngleType(ParameterType):
        """A SMIRNOFF angle type."""
        def __init__(self, node, parent):
            super(AngleType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.angle = _extract_quantity_from_xml_element(node, parent, 'angle')
            self.k = _extract_quantity_from_xml_element(node, parent, 'k')
            if 'fractional_bondorder' in parent.attrib:
                self.fractional_bondorder = parent.attrib['fractional_bondorder']
            else:
                self.fractional_bondorder = None

    _TAGNAME = 'AngleForce' # SMIRNOFF tag name to process
    _VALENCE_TYPE = 'Angle' # ChemicalEnvironment valence type expected for SMARTS
    _INFOTYPE = AngleType # class to hold force type info
    _OPENMMTYPE = openmm.HarmonicAngleForce # OpenMM force class to create

    def __init__(self, forcefield):
        super(AngleHandler, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        force = super(AngleHandler, self).createForce(system, topology, **kwargs)

        # Add all angles to the system.
        angles = self.getMatches(topology)
        skipped_constrained_angles = 0 # keep track of how many angles were constrained (and hence skipped)
        for (atoms, angle) in angles.items():
            # Get corresponding particle indices in Topology
            particle_indices = tuple([ atom.particle_index for atom in atoms ])

            # Ensure atoms are actually bonded correct pattern in Topology
            for (i,j) in [ (0,1), (1,2) ]:
                topology.assert_bonded(atoms[i], atoms[j])

            if topology.is_constrained(atoms[0], atoms[1]) and topology.is_constrained(atoms[1], atoms[2]) and topology.is_constrained(atoms[0], atoms[2]):
                # Angle is constrained; we don't need to add an angle term.
                skipped_constrained_angles += 1
                continue

            force.addAngle(*particle_indices, angle.angle, angle.k)

        logger.info('{} angles added ({} skipped due to constraints)'.format(len(angles) - skipped_constrained_angles, skipped_constrained_angles))

        # Check that no topological angles are missing force parameters
        _check_for_missing_valence_terms('AngleForce', topology, angles.keys(), topology.angles())

#=============================================================================================

## @private
class ProperTorsionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ProperTorsionForce>`` tags"""

    class ProperTorsionType(ParameterType):
        """A SMIRNOFF torsion type for proper torsions."""
        def __init__(self, node, parent):
            super(ProperTorsionType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.periodicity = list()
            self.phase = list()
            self.k = list()

            # Check that the SMIRKS pattern matches the type it's supposed to
            try:
                chemenv = ChemicalEnvironment(self.smirks)
                thistype = chemenv.getType()
                if thistype != 'ProperTorsion':
                    raise Exception("Error: SMIRKS pattern %s (parameter %s) does not specify a %s torsion, but it is supposed to." % (self.smirks, self.pid, 'Proper'))
            except SMIRKSParsingError:
                print("Warning: Could not confirm whether smirks pattern %s is a valid %s torsion." % (self.smirks, self.torsiontype))


            # TODO: Fractional bond orders should be processed on the per-force basis instead of per-bond basis
            if 'fractional_bondorder_method' in parent.attrib:
                self.fractional_bondorder = parent.attrib['fractional_bondorder']
            else:
                self.fractional_bondorder = None

            # Store parameters.
            index = 1
            while 'phase%d'%index in node.attrib:
                self.periodicity.append( int(_extract_quantity_from_xml_element(node, parent, 'periodicity%d' % index)) )
                self.phase.append( _extract_quantity_from_xml_element(node, parent, 'phase%d' % index, unit_name='phase_unit') )
                self.k.append( _extract_quantity_from_xml_element(node, parent, 'k%d' % index, unit_name='k_unit') )
                # Optionally handle 'idivf', which divides the periodicity by the specified value
                if ('idivf%d' % index) in node.attrib:
                    idivf = _extract_quantity_from_xml_element(node, parent, 'idivf%d' % index)
                    self.k[-1] /= float(idivf)
                index += 1

            # Check for errors, i.e. 'phase' instead of 'phase1'
            # TODO: Can we raise a more useful error if there is no ``id``?
            if len(self.phase)==0:
               raise Exception("Error: Torsion with id %s has no parseable phase entries." % self.pid)

    _TAGNAME = 'ProperTorsionForce' # SMIRNOFF tag name to process
    _VALENCE_TYPE = 'ProperTorsion' # ChemicalEnvironment valence type expected for SMARTS
    _INFOTYPE = ProperTorsionType # info type to store
    _OPENMMTYPE = openmm.PeriodicTorsionForce # OpenMM force class to create

    def __init__(self, forcefield):
        super(ProperTorsionHandler, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        force = super(ProperTorsionHandler, self).createForce(system, topology, **kwargs)

        # Add all proper torsions to the system.
        torsions = self.getMatches(topology)
        for (atoms, torsion) in torsions.items():
            # Ensure atoms are actually bonded correct pattern in Topology
            for (i,j) in [ (0,1), (1,2), (2,3) ]:
                topology.assert_bonded(atoms[i], atoms[j])

            for (periodicity, phase, k) in zip(torsion.periodicity, torsion.phase, torsion.k):
                force.addTorsion(atom_indices[0], atom_indices[1], atom_indices[2], atom_indices[3], periodicity, phase, k)

        logger.info('{} torsions added'.format(len(torsions)))

        # Check that no topological torsions are missing force parameters
        _check_for_missing_valence_terms('ProperTorsionForce', topology, torsions.keys(), topology.torsions())

## @private
class ImproperTorsionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<ImproperTorsionForce>`` tags"""

    class ImproperTorsionType(ParameterType):
        """A SMIRNOFF torsion type for improper torsions."""
        def __init__(self, node, parent):
            super(ImproperTorsionType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.periodicity = list()
            self.phase = list()
            self.k = list()

            # Check that the SMIRKS pattern matches the type it's supposed to
            try:
                chemenv = ChemicalEnvironment(self.smirks)
                thistype = chemenv.getType()
                if thistype != 'ImproperTorsion':
                    raise Exception("Error: SMIRKS pattern %s (parameter %s) does not specify a %s torsion, but it is supposed to." % (self.smirks, self.pid, 'Improper'))
            except SMIRKSParsingError:
                print("Warning: Could not confirm whether smirks pattern %s is a valid %s torsion." % (self.smirks, self.torsiontype))

            if 'fractional_bondorder' in parent.attrib:
                self.fractional_bondorder = parent.attrib['fractional_bondorder']
            else:
                self.fractional_bondorder = None

            # Store parameters.
            index = 1
            while 'phase%d'%index in node.attrib:
                self.periodicity.append( int(_extract_quantity_from_xml_element(node, parent, 'periodicity%d' % index)) )
                self.phase.append( _extract_quantity_from_xml_element(node, parent, 'phase%d' % index, unit_name='phase_unit') )
                self.k.append( _extract_quantity_from_xml_element(node, parent, 'k%d' % index, unit_name='k_unit') )
                # Optionally handle 'idivf', which divides the periodicity by the specified value
                if ('idivf%d' % index) in node.attrib:
                    idivf = _extract_quantity_from_xml_element(node, parent, 'idivf%d' % index)
                    self.k[-1] /= float(idivf)
                index += 1
                # SMIRNOFF applies trefoil (three-fold, because of right-hand rule) impropers unlike AMBER
                # If it's an improper, divide by the factor of three internally
                if node.tag=='Improper':
                    self.k[-1] /= 3.
            # Check for errors, i.e. 'phase' instead of 'phase1'
            # TODO: What can we do if there is no ``id``?
            if len(self.phase)==0:
               raise Exception("Error: Torsion with id %s has no parseable phase entries." % self.pid)

    _TAGNAME = 'ImproperTorsionForce' # SMIRNOFF tag name to process
    _VALENCE_TYPE = 'ImproperTorsion' # ChemicalEnvironment valence type expected for SMARTS
    _INFOTYPE = ImproperTorsionType # info type to store
    _OPENMMTYPE = openmm.PeriodicTorsionForce # OpenMM force class to create

    def __init__(self, forcefield):
        super(ImproperTorsionHandler, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        force = super(ImproperTorsionHandler, self).createForce(system, topology, **kwargs)

        # Add all improper torsions to the system
        torsions = self.getMatches(topology)
        for (atom_indices, improper) in impropers.items():
            # Ensure atoms are actually bonded correct pattern in Topology
            # For impropers, central atom is atom 1
            for (i,j) in [ (0,1), (1,2), (1,3) ]:
                topology.assert_bonded(atoms[i], atoms[j])

            # Impropers are applied in three paths around the trefoil having the same handedness
            for (periodicity, phase, k) in zip(improper.periodicity, improper.phase, improper.k):
                # Permute non-central atoms
                others = [ atom_indices[0], atom_indices[2], atom_indices[3] ]
                for p in [ (others[i], others[j], others[k]) for (i,j,k) in [(0,1,2), (1,2,0), (2,0,1)] ]
                    force.addTorsion(atom_indices[1], p[0], p[1], p[2], periodicity, phase, k)

        logger.info('{} impropers added, each applied in a six-fold trefoil' % (len(impropers)))

        # Check that no topological torsions are missing force parameters
        _check_for_missing_valence_terms('ImproperTorsionForce', topology, torsions.keys(), topology.impropers())

## @private
class vdWHandler(ParameterHandler):
    """Handle SMIRNOFF ``<vdWForce>`` tags"""

    # TODO: Is this necessary
    SCALETOL = 1e-5

    # DEFAULTS
    _DEFAULTS = {
        'potential' : 'Lennard-Jones-12-6',
        'scale12' : 0.0,
        'scale13' : 0.0,
        'scale14' : 0.5,
        'scale15' : 1.0,
    }

    class vdWType(ParameterType):
        """A SMIRNOFF vdWForce type."""
        def __init__(self, node, parent):
            # NOTE: Currently we support radius definition via 'sigma' or 'rmin_half'.
            super(StericsType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields

            # Make sure we don't have BOTH rmin_half AND sigma
            try:
                a = _extract_quantity_from_xml_element(node, parent, 'sigma')
                a = _extract_quantity_from_xml_element(node, parent, 'rmin_half')
                raise Exception("Error: BOTH sigma and rmin_half cannot be specified simultaneously in the .offxml file.")
            except:
                pass

            # Handle Lennard-Jones sigma
            try:
                self.sigma = _extract_quantity_from_xml_element(node, parent, 'sigma')
            #Handle rmin_half, AMBER-style
            except:
                rmin_half = _extract_quantity_from_xml_element(node, parent, 'rmin_half', unit_name='sigma_unit')
                self.sigma = 2.*rmin_half/(2.**(1./6.))
            self.epsilon = _extract_quantity_from_xml_element(node, parent, 'epsilon')

    _TAGNAME = 'vdWForce' # SMIRNOFF tag name to process
    _INFOTYPE = vdWType # info type to store
    _OPENMMTYPE = openmm.NonbondedForce # OpenMM force class to create

    def __init__(self, forcefield):
        super(NonbondedParameterHandler, self).__init__(forcefield)


    # TODO: Handle the case where multiple <NonbondedForce> tags are found
    # if abs(Handler.coulomb14scale - float(element.attrib['coulomb14scale'])) > NonbondedHandler.SCALETOL or \
    #                 abs(Handler.lj14scale - float(element.attrib['lj14scale'])) > NonbondedHandler.SCALETOL:
    #            raise ValueError('Found multiple NonbondedForce tags with different 1-4 scales')
    #    for atom in element.findall('Atom'):
    #        Handler.registerAtom(atom, element)

    # TODO: nonbondedMethod and nonbondedCutoff should now be specified by StericsForce attributes
    def createForce(self, system, topology, nonbondedMethod=NoCutoff, nonbondedCutoff=0.9, **args):
        force = super(NonbondedParameterHandler, self).createForce(system, topology)

        methodMap = {NoCutoff:openmm.NonbondedForce.NoCutoff,
                     CutoffNonPeriodic:openmm.NonbondedForce.CutoffNonPeriodic,
                     CutoffPeriodic:openmm.NonbondedForce.CutoffPeriodic,
                     Ewald:openmm.NonbondedForce.Ewald,
                     PME:openmm.NonbondedForce.PME}
        if nonbondedMethod not in methodMap:
            raise ValueError('Illegal nonbonded method for NonbondedForce; method given was %s' % nonbondedMethod)
        force = openmm.NonbondedForce()
        force.setNonbondedMethod(methodMap[nonbondedMethod])
        force.setCutoffDistance(nonbondedCutoff)
        if 'ewaldErrorTolerance' in args:
            force.setEwaldErrorTolerance(args['ewaldErrorTolerance'])
        if 'useDispersionCorrection' in args:
            force.setUseDispersionCorrection(bool(args['useDispersionCorrection']))
        system.addForce(force)

        # Iterate over all defined Lennard-Jones types, allowing later matches to override earlier ones.
        atoms = self.getMatches(topology)

        # Create all particles.
        for particle in topology.particles:
            force.addParticle(0.0, 1.0, 0.0)

        # Set the particle Lennard-Jones terms.
        for (atoms, ljtype) in atoms.items():
            force.setParticleParameters(atoms[0].particle_index, 0.0, ljtype.sigma, ljtype.epsilon)

        # Check that no atoms are missing force parameters
        # QUESTION: Don't we want to allow atoms without force parameters? Or perhaps just *particles* without force parameters, but not atoms?
        _check_for_missing_valence_terms('NonbondedForce Lennard-Jones parameters', topology, atoms.keys(), topology.atoms)

        # Set the partial charges
        # TODO: We need to make sure we have already assigned partial charges to the Topology reference molecules
        for atom in topology.atoms:
            # Retrieve nonbonded parameters for reference atom (charge not set yet)
            _, sigma, epsilon = force.getParticleParameters(atom.particle_index)
            # Set the nonbonded force with the partial charge
            force.setParticleParameters(atom.particle_index, atom.charge, sigma, epsilon)

    def postprocessSystem(self, system, topology, **args):
        # Create exceptions based on bonds.
        # QUESTION: Will we want to do this for *all* cases, or would we ever want flexibility here?
        bond_particle_indices = [ (atom1.particle_index, atom2.particle_index) for (atom1, atom2) in topology.bonds ]
        for force in system.getForces():
            # TODO: Should we just store which `Force` object we are adding to and use that instead,
            # to prevent interference with other kinds of forces in the future?
            # TODO: Can we generalize this to allow for `CustomNonbondedForce` implementations too?
            if isinstance(force, openmm.NonbondedForce):
                nonbonded.createExceptionsFromBonds(bond_particle_indices, self.coulomb14scale, self.lj14scale)

## @private
class BondChargeCorrectionHandler(ParameterHandler):
    """Handle SMIRNOFF ``<BondChargeCorrection>`` tags"""

    class BondChargeCorrectionType(ParameterType):
        """A SMIRNOFF bond charge correction type."""
        def __init__(self, node, parent):
            super(BondChargeCorrectionHandler, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.increment = _extract_quantity_from_xml_element(node, parent, 'increment')
            # If no units are specified, assume elementary charge
            if type(self.increment) == float:
                self.increment *= unit.elementary_charge

    _TAGNAME = 'BondChargeCorrection' # SMIRNOFF tag name to process
    _INFOTYPE = BondChargeCorrectionType # info type to store
    _OPENMMTYPE = openmm.NonbondedForce # OpenMM force class to create or utilize

    def __init__(self, forcefield):
        super(BondChargeCorrectionHandler, self).__init__(forcefield)

    #if element.attrib['method'] != existing[0]._initialChargeMethod:
    #raise Exception("Existing BondChargeCorrectionHandler uses initial charging method '%s' while new BondChargeCorrectionHandler requests '%s'" % (existing[0]._initialChargeMethod, element.attrib['method']))

    def createForce(self, system, topology, **args):
        # No forces are created by this Handler.
        pass

    # TODO: Move chargeModel and library residue charges to SMIRNOFF spec
    def postprocessSystem(self, system, topology, **kwargs):
        bonds = self.getMatches(topology)

        # Apply bond charge increments to all appropriate force groups
        # QUESTION: Should we instead apply this to the Topology in a preprocessing step, prior to spreading out charge onto virtual sites?
        for force in system.getForces():
            if force.__class__.__name__ in ['NonbondedForce']: # TODO: We need to apply this to all Force types that involve charges
                for (atoms, bond) in bonds.items():
                    # Get corresponding particle indices in Topology
                    particle_indices = tuple([ atom.particle_index for atom in atoms ])
                    # Retrieve parameters
                    [charge0, sigma0, epsilon0] = force.getParticleParameters(particle_indices[0])
                    [charge1, sigma1, epsilon1] = force.getParticleParameters(particle_indices[1])
                    # Apply bond charge increment
                    charge0 -= bond.increment
                    charge1 += bond.increment
                    # Update charges
                    force.setParticleParameters(particle_indices[0], charge0, sigma0, epsilon0)
                    force.setParticleParameters(particle_indices[1], charge1, sigma1, epsilon1)

## @private
class GBSAParameterHandler(ParameterHandler):
    """Handle SMIRNOFF ``<GBSAParameterHandler>`` tags"""
    # TODO: Differentiate between global and per-particle parameters for each model.

    # Global parameters for surface area (SA) component of model
    SA_expected_parameters = {
        'ACE' : ['surface_area_penalty', 'solvent_radius'],
        None : [],
    }

    # Per-particle parameters for generalized Born (GB) model
    GB_expected_parameters = {
        'HCT' : ['radius', 'scale'],
        'OBC1' : ['radius', 'scale'],
        'OBC2' : ['radius', 'scale'],
    }

    class GBSAType(ParameterType):
        """A SMIRNOFF GBSA type."""
        def __init__(self, node, parent):
            super(GBSAType, self).__init__(node, parent)

            # Store model parameters.
            gb_model = parent.attrib['gb_model']
            expected_parameters = GBSAParameterHandler.GB_expected_parameters[gb_model]
            provided_parameters = list()
            missing_parameters = list()
            for name in expected_parameters:
                if name in node.attrib:
                    provided_parameters.append(name)
                    value = _extract_quantity_from_xml_element(node, parent, name)
                    setattr(self, name, value)
                else:
                    missing_parameters.append(name)
            if len(missing_parameters) > 0:
                msg  = 'GBSAForce: missing per-atom parameters for tag %s' % str(node)
                msg += 'model "%s" requires specification of per-atom parameters %s\n' % (gb_model, str(expected_parameters))
                msg += 'provided parameters : %s\n' % str(provided_parameters)
                msg += 'missing parameters: %s' % str(missing_parameters)
                raise Exception(msg)

    def __init__(self, forcefield):
        super(GBSAParameterHandler, self).__init__(forcefield)

    # TODO: Fix this
    def parseElement(self):
        # Initialize GB model
        gb_model = element.attrib['gb_model']
        valid_GB_models = GBSAParameterHandler.GB_expected_parameters.keys()
        if not gb_model in valid_GB_models:
            raise Exception('Specified GBSAForce model "%s" not one of valid models: %s' % (gb_model, valid_GB_models))
        self.gb_model = gb_model

        # Initialize SA model
        sa_model = element.attrib['sa_model']
        valid_SA_models = GBSAParameterHandler.SA_expected_parameters.keys()
        if not sa_model in valid_SA_models:
            raise Exception('Specified GBSAForce SA_model "%s" not one of valid models: %s' % (sa_model, valid_SA_models))
        self.sa_model = sa_model

        # Store parameters for GB and SA models
        # TODO: Deep copy?
        self.parameters = element.attrib

    # TODO: Generalize this to allow forces to know when their OpenMM Force objects can be combined
    def checkCompatibility(self, Handler):
        """
        Check compatibility of this Handler with another Handlers.
        """
        Handler = existing[0]
        if (Handler.gb_model != self.gb_model):
            raise ValueError('Found multiple GBSAForce tags with different GB model specifications')
        if (Handler.sa_model != self.sa_model):
            raise ValueError('Found multiple GBSAForce tags with different SA model specifications')
        # TODO: Check other attributes (parameters of GB and SA models) automatically?

    def createForce(self, system, topology, **args):
        # TODO: Rework this
        from openforcefield.typing.engines.smirnoff import gbsaforces
        force_class = getattr(gbsaforces, self.gb_model)
        force = force_class(**self.parameters)
        system.addForce(force)

        # Add all GBSA terms to the system.
        expected_parameters = GBSAParameterHandler.GB_expected_parameters[self.gb_model]

        # Create all particles with parameters set to zero
        atoms = self.getMatches(topology)
        nparams = 1 + len(expected_parameters) # charge + GBSA parameters
        params = [ 0.0 for i in range(nparams) ]
        for particle in topology.particles():
            force.addParticle(params)
        # Set the GBSA parameters (keeping charges at zero for now)
        for (atoms, gbsa_type) in atoms.items():
            atom = atoms[0]
            # Set per-particle parameters for assigned parameters
            params = [atom.charge] + [ getattr(gbsa_type, name) for name in expected_parameters ]
            force.setParticleParameters(atom.particle_index, params)
