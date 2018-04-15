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

# QUESTION: How should we support other non-XML representations in future?

class ForceField(object):
    """A factory initialized from a SMIRNOFF force field that constructs OpenMM System objects from corresponding openforcefield Topology objects.

    When a ``ForceField`` object is created from one or more specified XML files, all ``ForceGenerator`` subclasses currently imported are identified
    and registered to handle different sections of the SMIRNOFF forcefield definition file(s). The files are then processed by the generators to populate
    internal data structures.

    Processing a ``Topology`` object defining a chemical system will then call all ``ForceGenerator`` objects in an order guaranteed to satisfy the
    declared processing order constraints of each ``ForceGenerator``.

    A programmatic API can be used to modify loaded parameters or write out the new parameter set.

    Attributes
    ----------
    parsers : dict of str : openforcefield.typing.engines.smirnoff.ForceGenerator
        Registered list of parsers that will handle forcefield tags.
        parsers[tagname] is the ``ForceGenerator`` that will be called to process XML block ``tagname``
        Parsers are registered when the ForceField object is created.

    """

    def __init__(self, *files, force_generators=None, disable_version_check=False):
        """Create a new ForceField object from one or more SMIRNOFF XML parameter definition files.

        Parameters
        ----------
        files : list
            A list of XML files defining the SMIRNOFF force field.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory
            (for built in force fields), or an open file-like object with a read() method from which the forcefield XML data can be loaded.
            If multiple files are specified, any top-level tags that are repeated will be merged if they are compatible,
            with files appearing later in the sequence resulting in parameters that have higher precedence.
            Support for multiple files is primarily intended to allow solvent parameters to be specified by listing them last in the sequence.
        force_generators : iterable of ForceGenerator classes, optional, default=None
            If not None, the specified set of ForceGenerator objects will be used to parse the XML tags.
            By default, all imported subclasses of ForceGenerator are automatically registered to parse XML tags.
        disable_version_check : bool, optional, default=False
            If True, will disable checks against the current highest supported forcefield version.

        """
        self._aromaticity_model = None # the aromaticity model to use for parsing; denote uninitialized

        self.disable_version_check = disable_version_check # if True, we won't check which SMIRNOFF version number we're parsing

        # Load all XML files containing parameter definitions
        # TODO: We may not store XML files, but instead parse them right away
        self.load_file(files)

        # Register all ForceGenerator objects that will handle SMIRNOFF tags in processing XML files
        if force_generators is None:
            # Find all imported subclasses of ForceGenerator
            force_generators = self._find_force_generators()
        self._register_parsers(force_generators)

    @property
    def parsers(self):
        """
        Retrieve a read-only list of ForceGenerators listed as parsers.
        """
        return copy.deepcopy(self._parsers)

    # TODO: We need to change this to just find all ForceGenerator objects in this file;
    # otherwise, we can't define two different ForceGenerator subclasses to compare for a new type of energy term
    # since both will try to register themselves for the same XML tag and an Exception will be raised.
    @staticmethod
    def _find_force_generators():
        """Identify all imported subclasses of ForceGenerator.

        Returns
        -------
        force_generators : list of ForceGenerator subclasses
            List of ForceGenerator subclasses (not objects)
        """
        return all_subclasses(ForceGenerator)

    @staticmethod
    def _register_parsers(force_generators):
        """Register all force generators with their ``_TAGNAME`` tag namesself.

        Parameters
        ----------
        force_generators : list of ForceGenerator subclasses
            The force generators to register.

        """
        parsers = dict()
        for force_generator in force_generators:
            tagname = subclass._TAGNAME
            if tagname is not None:
                if tagname in parsers.keys():
                    raise Exception("ForceGenerator {} provides a parser for tag '{}', but ForceGemerator {} has already been registered to handle that tag.".format(subclass, tagname, self.parsers[tagname]))
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
            Name of the calling force generator
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

    def load_file(self, files):
        """Load a SMIRNOFF XML file and add the definitions from it to this ForceField.

        Parameters
        ----------
        files : string or file or tuple
            An XML file or tuple of XML files containing SMIRNOFF force field definitions.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in force fields), or an open file-like object with a read() method from which the forcefield XML data can be loaded.

        .. notes ::

           * New SMIRNOFF tags are handled independently, as if they were specified in the same file.
           * If a SMIRNOFF tag that has already been read appears again, its definitions are appended to the end of the previously-read
             definitions if the tags are configured with compatible attributes; otherwise, an ``IncompatibleTagException`` is raised.

        """
        # Ensure that we are working with a tuple of files.
        if not isinstance(files, tuple):
            files = (files,)

        # Load in all XML trees.
        # QUESTION: Should we allow users to specify forcefield URLs so they can pull forcefield definitions from the web too?
        trees = list()
        for file in files:
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

            trees.append(tree)

        # Retain XML trees internally
        self._XMLTrees = trees

        # Parse XML, get force definitions
        # QUESTION: Should we use lazy instantiation to only parse XML trees when we need to so that users will manipulate parameters only via XML tags?
        # QUESTION: If we only parse XML trees when needed, we may run into the problem where a call to ForceField(*files) is successful
        # and we *think* the files are valid, but we don't actually parse the file and find there is a problem until we actually use it.
        # Would this create problems for debugging? Maybe we can have an optional `parse_immediately=True` flag?
        self.parseXMLTrees()

    def save_file(self, filename, format="XML"):
        """Write the current forcefield parameter set to a file.

        Parameters
        ----------
        filename : str
            The name of the file to be written.
        format : str, optional, default="XML"
            The format of the file to be written.
            One of ['XML']

        """
        raise NotImplementedError('Feature implemented yet.')

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
        Initialize all registered ForceGenerator objects by parsing loaded XML files.

        .. notes ::

           * New SMIRNOFF tags are handled independently, as if they were specified in the same file.
           * If a SMIRNOFF tag that has already been read appears again, its definitions are appended to the end of the previously-read
             definitions if the tags are configured with compatible attributes; otherwise, an ``IncompatibleTagException`` is raised.

        """
        trees = self._XMLTrees

        # Create all ForceGenerator objects from scratch
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

    def getGenerators(self):
        """Get the list of all registered generators."""
        return self._forces

    def registerGenerator(self, generator):
        """Register a new generator."""
        # Special case: ConstraintGenerator has to come before BondGenerator and AngleGenerator.
        # TODO: Figure out a more general way to allow generators to specify enforced orderings via dependencies.
        if isinstance(generator, ConstraintGenerator):
            self._forces.insert(0, generator)
        else:
            self._forces.append(generator)

    # TODO: Rework the API for getting, setting, and adding parameters
    def getParameter(self, smirks=None, paramID=None, force_type='Implied'):
        """Get info associated with a particular parameter as specified by SMIRKS or parameter ID, and optionally force term.

        Parameters
        ----------
        smirks (optional) : str
            Default None. If specified, will pull parameters on line containing this `smirks`.
        paramID : str
            Default None. If specified, will pull parameters on line with this `id`
        force_type : str
            Default "Implied". Optionally, specify a particular force type such as
            "HarmonicBondForce" or "HarmonicAngleForce" etc. to search for a
            matching ID or SMIRKS.

        Returns
        -------
        params : dict
            Dictionary of attributes (parameters and their descriptions) from XML

        Usage notes: SMIRKS or parameter ID must be specified.

        .. todo::
           *  Update behavior of "Implied" force_type so it raises an exception if the parameter is not uniquely identified by the provided info.

        """
        # Check for valid input
        if smirks and paramID:
            raise ValueError("Error: Specify SMIRKS OR parameter ID but not both.")
        if (smirks is None) and (paramID is None):
            raise ValueError("Error: Must specify SMIRKS or parameter ID.")

        trees = self._XMLTrees
        # Loop over XML files we read
        for tree in trees:
            # Loop over tree
            for child in tree.getroot():
                # Check a particular section?
                checksection = True
                if force_type is not 'Implied':
                    # See whether this has the tag we want to check
                    checksection = (child.tag==force_type)

                if checksection:
                    #Loop over descendants
                    for elem in child.iterdescendants(tag=etree.Element):
                        if (smirks and elem.attrib['smirks']==smirks) or (paramID and elem.attrib['id']==paramID):
                            return copy.deepcopy(elem.attrib)

    # TODO: Rework the API for getting, setting, and adding parameters
    def setParameter(self, params, smirks=None, paramID=None, force_type="Implied"):
        """Get info associated with a particular parameter as specified by SMIRKS or parameter ID, and optionally force term.

        Parameters
        ----------
        params : dict
            Dictionary of attributes (parameters and their descriptions) for XML,
            i.e. as output by getParameter.
        smirks (optional) : str
            Default None. If specified, will set parameters on line containing this `smirks`.
        paramID (optional) : str
            Default None. If specified, will set parameters on line with this `id`
        force_type (optional) : str
            Default "Implied". Optionally, specify a particular force type such as
            "HarmonicBondForce" or "HarmonicAngleForce" etc. to search for a
            matching ID or SMIRKS.

        Returns
        -------
        success : bool
            True/False as to whether that parameter was found and successfully set

        Usage notes: SMIRKS or parameter ID must be specified.

        .. todo::
            * Update set/get parameter API.
            * Update behavior of "Implied" force_type so it raises an exception if the parameter is not uniquely identified by the provided info.

        """
        # Check for valid input
        if smirks and paramID:
            raise ValueError("Error: Specify SMIRKS OR parameter ID but not both.")
        if (smirks is None) and (paramID is None):
            raise ValueError("Error: Must specify SMIRKS or parameter ID.")
        if not params:
            raise ValueError("Error, parameters must be specified.")

        # Below, we should do due dilegence that we're working on a parameter line which has
        # roughly the same types of parameters (though for torsions (or bonds, if we have partial bond orders),
        # the number of terms might differ), so define a utility function to give
        # back the basic names of parameters (i.e. k) without suffixes
        def get_param_names( param_keys ):
            names = set()
            for param in param_keys:
                ct = 1
                while param[-ct].isdigit():
                    ct+=1
                if ct > 1:
                    names.add( param[:-(ct-1)])
                else:
                    names.add( param )
            return names

        # Set parameters
        trees = self._XMLTrees
        success = False
        # Loop over XML files we read
        for tree in trees:
            # Loop over tree
            for child in tree.getroot():
                # Check a particular section?
                checksection = True
                if force_type is not 'Implied':
                    # See whether this has the tag we want to check
                    checksection= (child.tag==force_type)

                if checksection:
                    #Loop over descendants
                    for elem in child.iterdescendants(tag=etree.Element):
                        if (smirks and elem.attrib['smirks']==smirks) or (paramID and elem.attrib['id']==paramID):
                            # Try to set parameters
                            old_params=elem.attrib
                            if get_param_names(old_params.keys()) != get_param_names(params.keys()):
                                raise ValueError('Error: Provided parameters have different keys (%s) than existing parameters (%s).' % (', '.join(old_params.keys()), ', '.join(params.keys())))

                            # Loop over attributes, change values
                            for tag in params.keys():
                                elem.set( tag, params[tag])

                            # Found parameters and set, so update success flag
                            success = True

        # If we made any changes to XML, set flag so it will be reprocessed prior to system creation
        if success:
            self._XMLModified = True

        return success

    # TODO: Rework the API for getting, setting, and adding parameters
    def addParameter(self, params, smirks, force_type, tag):
        """Add specified SMIRKS/parameter in the section under the specified force type.

        Parameters
        ----------
        params : dict
            Dictionary of attributes (parameters and their descriptions) for XML,
            i.e. as output by getParameter.
        smirks : str
            SMIRKS pattern to associate with this parameter
        force_type : str
            Specify a particular force type such as "HarmonicBondForce" or "HarmonicAngleForce" in which to add this parameter
        tag : str
            Tag to use identifying this parameter, i.e. 'Bond' for a HarmonicBondForce, etc.

        Returns
        -------
        success : bool
            Returns ``True`` on success, or ``False`` on failure.

        .. todo::
            * Update set/get parameter API.

        """

        trees = self._XMLTrees
        # Loop over XML files we read
        success = False
        for tree in trees:
            # Loop over tree
            for child in tree.getroot():
                if child.tag==force_type:
                    success = True
                    addl = etree.Element( tag, smirks=smirks, attrib = params)
                    child.append(addl)

        # If we made any changes to XML, set flag so it will be reprocessed prior to system creation
        if success:
            self._XMLModified = True

        return success

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

    def _resolve_force_generator_order(self):
        """Resolve the order in which ForceGenerator objects should execute to satisfy constraints.
        """
        ordered_force_generators = list()
        # Create a DAG expressing dependencies
        import networkx as nx
        G = nx.DiGraph()
        for force_generator in self.parsers.items():
            G.add_node(force_generator._TAGNAME)
            if force_generator._DEPENDENCIES is not None:
                for dependency in force_generator._DEPENDENCIES:
                    G.add_edge(dependency._TAGNAME, force_generator._TAGNAME)
        # TODO: Check to make sure DAG isn't cyclic
        # Resolve order
        ordered_force_generators = [ self.parsers[tagname] for tagname in nx.topological_sort(G) ]
        return ordered_force_generators

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
            A ``ValueError`` will be raised if parameters are specified that are not used by any force generators.

        Returns
        -------
        system : simtk.openmm.System
            The newly created System

        """
        # TODO: Have `verbose` flag set whether logging info is displayed or not.

        # Re-parse XML by ForceGenerator objects if loaded XML files have been added or modified
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

        # Determine the order in which to process ForceGenerator objects in order to satisfy dependencies
        force_generators = self._resolve_force_generator_order()

        # Check if any kwargs have been provided that aren't handled by force generators
        known_kwargs = set()
        for force_generator in force_generators:
            known_args.update(force_generator.known_kwargs)
        unknown_kwargs = set(kwargs.keys()).difference(known_kwargs)
        if len(unknown_kwargs) > 0:
            msg = "The following keyword arguments to createSystem() are not used by any registered force generator: {}\n".format(unknown_kwargs)
            msg += "Known keyword arguments: {}".format(known_kwargs)
            raise ValueError(msg)

        # Add forces to the System
        for force_generator in force_generators:
            force_generator.createForce(system, topology, **kwargs)

        # Let force generators do postprocessing
        for force_generator in force_generators:
            force_generator.postprocessSystem(system, topology, **kwargs)

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
# The following classes are generators that know how to create Force subclasses and add them to a System that is being
# created.  Each generator class must define three methods: 1) a static method that takes an etree Element and a ForceField,
# and returns the corresponding generator object; 2) a createForce() method that constructs the Force object and adds it
# to the System; and 3) a labelForce() method that provides access to which
# terms are applied to which atoms in specified oemols.
# The static method should be added to the parsers map.
#=============================================================================================

#=============================================================================================
# Force generators
#=============================================================================================

## @private
class ForceGenerator(object):
    """Base class for force generators.
    """
    _TAGNAME = None # str of XML element name handled by this ForceGenerator
    _VALENCE_TYPE = None # ChemicalEnvironment valence type string expected by SMARTS string for this Generator
    _INFOTYPE = None # container class with type information that will be stored in self._types
    _OPENMMTYPE = None # OpenMM Force class (or None if no equivalent)
    _DEPENDENCIES = None # list of ForceGenerator classes that must precede this, or None
    _DEFAULTS = {} # dict of attributes and their default values at tag-level
    _KWARGS = [] # list of keyword arguments accepted by the force generator on initialization
    _SMIRNOFF_VERSION_INTRODUCED = 0.0 # the earliest version of SMIRNOFF spec that supports this ForceGenerator
    _SMIRNOFF_VERSION_DEPRECATED = None # if deprecated, the first SMIRNOFF version number it is no longer used
    _REQUIRE_UNITS = None # list of parameters that require units to be defined

    def __init__(self, forcefield):
        self.ff = forcefield # the ForceField object that this ForceGenerator is registered with
        self._types = list() # list of ForceType objects of type cls._INFOTYPE

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
            matches[atoms] is the ForceType object corresponding to the tuple of Atom objects ``Atoms``

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
        Parse the XML tag/section this ForceGenerator is registered for.

        SMIRNOFF sections may be split across multiple files or otherwise appear multiple times,
        so we need to be able to handle multiple calls to parseElement().

        Parameters
        ----------


        """
        existing = [f for f in ff._forces if isinstance(f, cls)]
        if len(existing) == 0:
            generator = cls(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]

        # Extract all tag-level attributes (or their defaults)
        for attribute in _DEFAULTS.keys():
            setattr(self, attribute, _extract_quantity_from_xml_element(node, parent, attribute, default=DEFAULTS[attribute]))

        # Register all SMIRNOFF definitions for all occurrences of its registered tag
        for section in element.findall(tag):
            generator.registerType(section, element)

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

## @private
class ForceType(object):
    """
    Base class for force types.
    """
    def __init__(self, node, parent):
        """
        Create a ForceType that stores ``smirks`` and ``id`` attributes.
        """
        self.smirks = _validate_smarts(node.attrib['smirks'], node=node, valence_type=_VALENCE_TYPE)
        if 'id' in node.attrib:
            self.pid = _extract_quantity_from_xml_element(node, parent, 'id')

#=============================================================================================

## @private
class ConstraintGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<Constraints>`` tags

    ``ConstraintGenerator`` must be applied before ``HarmonicBondGenerator`` and ``HarmonicAngleGenerator``,
    since those classes add constraints for which equilibrium geometries are needed from those tags.
    """
    class ConstraintType(ForceType):
        """A SMIRNOFF constraint type"""
        def __init__(self, node, parent):
            super(ConstraintType, self).__init__(node, parent) # Base class handles ``smirks`` and ``id`` fields
            if 'distance' in node.attrib:
                self.distance = _extract_quantity_from_xml_element(node, parent, 'distance') # Constraint with specified distance will be added by ConstraintGenerator
            else:
                self.distance = True # Constraint to equilibrium bond length will be added by HarmonicBondGenerator

    _TAGNAME = 'Constraint'
    _VALENCE_TYPE = 'Bond' # ChemicalEnvironment valence type expected for SMARTS # TODO: Do we support more exotic types as well?
    _INFOTYPE = ConstraintType
    _OPENMMTYPE = None # don't create a corresponding OpenMM Force class
    _REQUIRED_UNITS = ['distance']

    def __init__(self, forcefield):
        super(ConstraintGenerator, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        constraints = self.getMatches(topology)
        for (atoms, constraint) in constraints.items():
            # Update constrained atom pairs in topology
            topology.add_constraint(*atoms, constraint.distance)
            # If a distance is specified (constraint.distance != True), add the constraint here.
            # Otherwise, the equilibrium bond length will be used to constrain the atoms in HarmonicBondGenerator
            if constraint.distance is not True:
                system.addConstraint(*atoms, constraint.distance)

#=============================================================================================

## @private
class BondGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<BondForce>`` tags"""

    class BondType(ForceType):
        """A SMIRNOFF bond type"""
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
    _DEPENDENCIES = [ConstraintGenerator] # ConstraintGenerator must be executed first

    def __init__(self, forcefield):
        super(HarmonicBondGenerator, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        # Create or retrieve existing OpenMM Force object
        force = super(BondGenerator, self).createForce(system, topology, **kwargs)

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
class AngleGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<AngleForce>`` tags"""

    class AngleType(ForceType):
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
        super(AngleGenerator, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        force = super(AngleGenerator, self).createForce(system, topology, **kwargs)

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
class ProperTorsionGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<ProperTorsionForce>`` tags"""

    class ProperTorsionType(ForceType):
        """A SMIRNOFF torsion type for proper torsions."""
        def __init__(self, node, parent):
            super(ProperTorsionType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.periodicity = list()
            self.phase = list()
            self.k = list()

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
        super(ProperTorsionGenerator, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        force = super(ProperTorsionGenerator, self).createForce(system, topology, **kwargs)

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
class ImproperTorsionGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<ImproperTorsionForce>`` tags"""

    class ImproperTorsionType(ForceType):
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
                if thistype != 'Improper':
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
                # SMIRNOFF applies trefoil (six-fold) impropers unlike AMBER
                # If it's an improper, divide by the factor of six internally
                if node.tag=='Improper':
                    self.k[-1] /= 6.

            # Check for errors, i.e. 'phase' instead of 'phase1'
            # TODO: What can we do if there is no ``id``?
            if len(self.phase)==0:
               raise Exception("Error: Torsion with id %s has no parseable phase entries." % self.pid)

    _TAGNAME = 'ImproperTorsionForce' # SMIRNOFF tag name to process
    _VALENCE_TYPE = 'ImproperTorsion' # ChemicalEnvironment valence type expected for SMARTS
    _INFOTYPE = ImproperTorsionType # info type to store
    _OPENMMTYPE = openmm.PeriodicTorsionForce # OpenMM force class to create

    def __init__(self, forcefield):
        super(ImproperTorsionGenerator, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        force = super(ImproperTorsionGenerator, self).createForce(system, topology, **kwargs)

        # Add all improper torsions to the system
        torsions = self.getMatches(topology)
        for (atom_indices, improper) in impropers.items():
            # Ensure atoms are actually bonded correct pattern in Topology
            # For impropers, central atom is atom 1
            for (i,j) in [ (0,1), (1,2), (1,3) ]:
                topology.assert_bonded(atoms[i], atoms[j])

            # Impropers are applied to all six paths around the trefoil
            for (periodicity, phase, k) in zip(improper.periodicity, improper.phase, improper.k):
                # Permute non-central atoms
                others = [ atom_indices[0], atom_indices[2], atom_indices[3] ]
                for p in itertools.permutations( others ):
                    force.addTorsion(p[0], atom_indices[1], p[1], p[2], periodicity, phase, k)

        logger.info('{} impropers added, each applied in a six-fold trefoil' % (len(impropers)))

        # Check that no topological torsions are missing force parameters
        _check_for_missing_valence_terms('ImproperTorsionForce', topology, torsions.keys(), topology.impropers())

## @private
class vdWGenerator(ForceGenerator):
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

    class vdWType(ForceType):
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
        super(NonbondedForceGenerator, self).__init__(forcefield)


    # TODO: Handle the case where multiple <NonbondedForce> tags are found
    # if abs(generator.coulomb14scale - float(element.attrib['coulomb14scale'])) > NonbondedGenerator.SCALETOL or \
    #                 abs(generator.lj14scale - float(element.attrib['lj14scale'])) > NonbondedGenerator.SCALETOL:
    #            raise ValueError('Found multiple NonbondedForce tags with different 1-4 scales')
    #    for atom in element.findall('Atom'):
    #        generator.registerAtom(atom, element)

    # TODO: nonbondedMethod and nonbondedCutoff should now be specified by StericsForce attributes
    def createForce(self, system, topology, nonbondedMethod=NoCutoff, nonbondedCutoff=0.9, **args):
        force = super(NonbondedForceGenerator, self).createForce(system, topology)

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
class BondChargeCorrectionGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<BondChargeCorrection>`` tags"""

    class BondChargeCorrectionType(ForceType):
        """A SMIRNOFF bond charge correction type."""
        def __init__(self, node, parent):
            super(BondChargeCorrectionGenerator, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.increment = _extract_quantity_from_xml_element(node, parent, 'increment')
            # If no units are specified, assume elementary charge
            if type(self.increment) == float:
                self.increment *= unit.elementary_charge

    _TAGNAME = 'BondChargeCorrection' # SMIRNOFF tag name to process
    _INFOTYPE = BondChargeCorrectionType # info type to store
    _OPENMMTYPE = openmm.NonbondedForce # OpenMM force class to create or utilize

    def __init__(self, forcefield):
        super(BondChargeCorrectionGenerator, self).__init__(forcefield)

    #if element.attrib['method'] != existing[0]._initialChargeMethod:
    #raise Exception("Existing BondChargeCorrectionGenerator uses initial charging method '%s' while new BondChargeCorrectionGenerator requests '%s'" % (existing[0]._initialChargeMethod, element.attrib['method']))

    def createForce(self, system, topology, **args):
        # No forces are created by this generator.
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
class GBSAForceGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<GBSAForceGenerator>`` tags"""
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

    class GBSAType(ForceType):
        """A SMIRNOFF GBSA type."""
        def __init__(self, node, parent):
            super(GBSAType, self).__init__(node, parent)

            # Store model parameters.
            gb_model = parent.attrib['gb_model']
            expected_parameters = GBSAForceGenerator.GB_expected_parameters[gb_model]
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
        super(GBSAForceGenerator, self).__init__(forcefield)

    # TODO: Fix this
    def parseElement(self):
        # Initialize GB model
        gb_model = element.attrib['gb_model']
        valid_GB_models = GBSAForceGenerator.GB_expected_parameters.keys()
        if not gb_model in valid_GB_models:
            raise Exception('Specified GBSAForce model "%s" not one of valid models: %s' % (gb_model, valid_GB_models))
        self.gb_model = gb_model

        # Initialize SA model
        sa_model = element.attrib['sa_model']
        valid_SA_models = GBSAForceGenerator.SA_expected_parameters.keys()
        if not sa_model in valid_SA_models:
            raise Exception('Specified GBSAForce SA_model "%s" not one of valid models: %s' % (sa_model, valid_SA_models))
        self.sa_model = sa_model

        # Store parameters for GB and SA models
        # TODO: Deep copy?
        self.parameters = element.attrib

    # TODO: Generalize this to allow forces to know when their OpenMM Force objects can be combined
    def checkCompatibility(self, generator):
        """
        Check compatibility of this generator with another generators.
        """
        generator = existing[0]
        if (generator.gb_model != self.gb_model):
            raise ValueError('Found multiple GBSAForce tags with different GB model specifications')
        if (generator.sa_model != self.sa_model):
            raise ValueError('Found multiple GBSAForce tags with different SA model specifications')
        # TODO: Check other attributes (parameters of GB and SA models) automatically?

    def createForce(self, system, topology, **args):
        # TODO: Rework this
        from openforcefield.typing.engines.smirnoff import gbsaforces
        force_class = getattr(gbsaforces, self.gb_model)
        force = force_class(**self.parameters)
        system.addForce(force)

        # Add all GBSA terms to the system.
        expected_parameters = GBSAForceGenerator.GB_expected_parameters[self.gb_model]

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
