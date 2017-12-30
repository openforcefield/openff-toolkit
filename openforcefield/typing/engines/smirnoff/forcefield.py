#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Parser for the SMIRNOFF (SMIRKS Native Open Force Field) format.

.. codeauthor:: John D. Chodera <john.chodera@choderalab.org>
.. codeauthor:: David L. Mobley <dmobley@mobleylab.org>
.. codeauthor:: Peter K. Eastman <peastman@stanford.edu>

.. todo::
    * Use xml parser with 'sourceline' node attributes to aid debugging
    http://stackoverflow.com/questions/6949395/is-there-a-way-to-get-a-line-number-from-an-elementtree-element
    * Use logger instead of debug printing

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

import numpy

import lxml.etree as etree

from simtk import openmm, unit
from simtk.openmm.app import element as elem

from openforcefield.utils import get_data_filename
from openforcefield.topology import Topology
from openforcefield.typing.chemistry import ChemicalEnvironment, SMIRKSParsingError # TODO: Get rid of this?

#=============================================================================================
# CONFIGURE LOGGER
#=============================================================================================

logger = logging.getLogger(__name__)

#=============================================================================================
# ENUMERATED TYPES
#=============================================================================================

# Enumerated values for nonbonded method

class NoCutoff(object):
    def __repr__(self):
        return 'NoCutoff'
NoCutoff = NoCutoff()

class CutoffNonPeriodic(object):
    def __repr__(self):
        return 'CutoffNonPeriodic'
CutoffNonPeriodic = CutoffNonPeriodic()

class CutoffPeriodic(object):
    def __repr__(self):
        return 'CutoffPeriodic'
CutoffPeriodic = CutoffPeriodic()

class Ewald(object):
    def __repr__(self):
        return 'Ewald'
Ewald = Ewald()

class PME(object):
    def __repr__(self):
        return 'PME'
PME = PME()

# Enumerated values for constraint type

class HBonds(object):
    def __repr__(self):
        return 'HBonds'
HBonds = HBonds()

class AllBonds(object):
    def __repr__(self):
        return 'AllBonds'
AllBonds = AllBonds()

class HAngles(object):
    def __repr__(self):
        return 'HAngles'
HAngles = HAngles()

#=============================================================================================
# PRIVATE METHODS
#=============================================================================================

def _validateSMIRKS(smirks, node=None):
    """Validate the specified SMIRKS string.

    Parameters
    ----------
    smirks : str
       The SMIRKS string to be validated
    node : xml.etree.ElementTree.Element
       Node of etree with 'sourceline' attribute.

    .. todo::

        * Can we make this not dependent on OEChem?
        * Should we move this into ``openforcefield.topology.Topology`` or ``openforcefield.chemistry``?

    """
    from openeye import oechem

    qmol = oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smirks):
        if (node is not None) and ('sourceline' in node.attrib):
            raise Exception("Line %s: Error parsing SMIRKS '%s'" % (node.attrib['sourceline'], node.attrib['smirks']))
        else:
            raise Exception("Error parsing SMIRKS '%s'" % (node.attrib['smirks']))

    return smirks

def _extractQuantity(node, parent, name, unit_name=None, default=None):
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
        if default is not None:
            if 'sourceline' in node.attrib:
                raise Exception("Line %d : Expected XML attribute '%s' not found" % (node.attrib['sourceline'], name))
            else:
                raise Exception("Expected XML attribute '%s' not found" % (name))
        else:
            return default

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
        # TODO: This is very dangerous.
        string = '(%s * %s).value_in_unit_system(md_unit_system)' % (node.attrib[name], parent.attrib[unit_name])
        quantity = eval(string, unit.__dict__)

    return quantity

def _check_for_missing_valence_terms(name, topology, assigned_terms, topological_terms):
    """
    Check to ensure there are no missing valence terms.

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
    # Convert to lists
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
        msg = '%s: Mismatch between valence terms added and topological terms expected.\n' % name
        atoms = [ atom for atom in topology.atoms() ]
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
        raise Exception(msg)

def assert_bonded(topology, atom1, atom2):
    """
    Raise an exception if the specified atoms are not bonded in the topology.

    Parameters
    ----------
    topology : openforcefield.topology.Topology
        The Topology object to check for bonded atoms
    atom1, atom2 : openforcefield.topology.Atom
        The atoms to check to ensure they are bonded
    """
    assert topology.is_bonded(atom1, atom2), 'Atoms {} and {} are not bonded in topology'.format(atom1, atom2)

#=============================================================================================
# FORCEFIELD
#=============================================================================================

# A map of functions to parse elements of the XML file.
parsers = {}

class ForceField(object):
    """A ForceField constructs OpenMM System objects based on a Topology.
    """

    def __init__(self, *files):
        """Load one or more XML parameter definition files and create a SMIRNOFF ForceField object based on them.

        Parameters
        ----------
        files : list
            A list of XML files defining the SMIRNOFF force field.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in force fields), or an open file-like object with a read() method from which the forcefield XML data can be loaded.

        """
        self._forces = []
        self.loadFile(files)

    def loadFile(self, files):
        """Load a SMIRNOFF XML file and add the definitions from it to this ForceField.

        Parameters
        ----------
        files : string or file or tuple
            An XML file or tuple of XML files containing SMIRNOFF force field definitions.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in force fields), or an open file-like object with a read() method from which the forcefield XML data can be loaded.
        """

        # Ensure that we are working with a tuple of files.
        if not isinstance(files, tuple):
            files = (files,)

        # Load in all XML trees.
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
        self.parseXMLTrees()

    def parseXMLTrees(self):
        """Parse XML trees, load force definitions."""

        trees = self._XMLTrees

        # We'll be creating all forces again from scratch by re-parsing
        self._forces = []

        # Store forcefield version info and, if present, aromaticity model
        root = trees[0].getroot()
        # Formerly known as SMIRFF, now SMIRNOFF
        if root.tag=='SMIRFF' or root.tag=='SMIRNOFF':
            if 'version' in root.attrib:
                #TO DO: Should this be a float, a string, or something else?
                self.version = float(root.attrib['version'])
            else:
                self.version = 0.0
            if 'aromaticity_model' in root.attrib:
                self._aromaticity_model = root.attrib['aromaticity_model']
            else:
                self._aromaticity_model = None

            # QUESTION: Do we really want to specify whether or not to use fractional bond order at the root level when we can easily
            # determine whether individual Force components use it? Perhaps just specify a ``fractional_bondorder_method`` or ``model`` at root level?
            if 'use_fractional_bondorder' in root.attrib:
                self._use_fractional_bondorder = root.attrib['use_fractional_bondorder'] == 'True'
            else:
                self._use_fractional_bondorder = False
        else:
            raise ValueError("Error: ForceField parses a SMIRNOFF forcefield, but this does not appear to be one as the root tag is %s." % root.tag)

        # Load force definitions
        for tree in trees:
            root = tree.getroot()

            # Before loading, do some error checking/consistency checking.
            # Warn if version number is not consistent
            if 'version' in root.attrib:
                if float(root.attrib['version']) != self.version:
                    print("Warning: Inconsistent version number in parsed FFXML files.")
            # Throw an exception if aromaticity model is not consistent
            if 'aromaticity_model' in root.attrib:
                if root.attrib['aromaticity_model'] != self._aromaticity_model:
                    raise ValueError("Error: Aromaticity model specified in FFXML files is inconsistent.")

            # Now actually load
            for child in root:
                if child.tag in parsers:
                    parsers[child.tag](child.tag, child, self)

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
        # Special case: ConstraintGenerator has to come before HarmonicBondGenerator and HarmonicAngleGenerator.
        # TODO: Figure out a more general way to allow generators to specify enforced orderings.
        if isinstance(generator, ConstraintGenerator):
            self._forces.insert(0, generator)
        else:
            self._forces.append(generator)

    # TODO: Rework the API
    def getParameter(self, smirks = None, paramID=None, force_type='Implied'):
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

    # TODO: Rework the API
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

    # TODO: Rework the API
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

    def writeFile(self, files):
        """Write forcefield trees out to specified files.

        .. todo::
            * Overhaul API

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

    def createSystem(self, topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*unit.nanometer,
                     constraints=None, rigidWater=True, removeCMMotion=True, hydrogenMass=None, residueTemplates=dict(),
                     chargeMethod=None, verbose=False, **kwargs):
        """Construct an OpenMM System representing a Topology with this force field. XML will be re-parsed if it is modified prior to system creation.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The ``Topology`` corresponding to the ``System`` object to be created.
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : boolean=True
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep
            their total mass the same.
        residueTemplates : dict=dict()
           Key: Topology Residue object
           Value: string, name of _TemplateData residue template object to use for
                  (Key) residue
           This allows user to specify which template to apply to particular Residues
           in the event that multiple matching templates are available (e.g Fe2+ and Fe3+
           templates in the ForceField for a monoatomic iron ion in the topology).
        chargeMethod : str, optional, default=None
           If 'BCC' is specified, bond charge corrections defined the `ForceField` will be applied to AM1-derived charges, otherwise charges from provided `molecules` will be used. (DEFAULT)
           If one of the `openeye.oequacpac.OECharges_` options is specified as a string (e.g. 'OECharges_AM1BCCSym'), this will be used and no bond charge corrections will be applied.
           If `None`, charges from the provided `molecules` will be used and no bond charge corrections will be applied.
        verbose : bool
           If True, verbose output will be printed.
        kwargs
             Arbitrary additional keyword arguments may also be specified.
             This allows extra parameters to be specified that are specific to
             particular force fields.

        Returns
        -------
        system : simtk.openmm.System
            The newly created System

        .. todo::

            * Replace ``verbose`` with log level selection.

            * Should the ``Topology`` be modified by ``createSystem``?
              There are a number of optional argument that we may want to make only defined by the SMIRNOFF forcefield:

                * ``chargeMethod``
                * ``hydrogenMass``
                * ``constraints``
                * ``rigidWater``

        """
        # Re-parse XML by generators if XML was modified
        self._reparse_xml_if_needed()

        # Make a deep copy of the topology so we don't accidentaly modify it
        topology = copy.deepcopy(topology)

        # Set aromaticity model to that used by forcefield
        topology.set_aromaticity_model(self._aromaticity_model)

        # Charge molecules if a valid charging method was identified
        # TODO: chargeMethod should be part of the SMIRNOFF spec.
        if chargeMethod is not None:
            for molecule in topology.molecules:
                topology.assign_partial_charges(method=chargeMethod)

        # Update bond orders stored in the topology if needed
        if self.bondOrderMethod is not None:
            for molecule in topology.molecules:
                topology.assign_fractional_bond_orders(method=self.bondOrderMethod)

        # Create the System and add atoms
        system = openmm.System()
        for atom in topology.atoms():
            # Add the particle to the OpenMM system.
            # QUESTION: Do we need an option to have SMARTS-specified fractional masses for compatibility with other forcefields?
            system.addParticle(atom.element.mass)

        # Adjust hydrogen masses if requested.
        # QUESTION: Do we need the capability to modify hydrogen masses via a keyword argument?
        if hydrogenMass is not None:
            if not unit.is_quantity(hydrogenMass):
                hydrogenMass *= unit.dalton
            for atom1, atom2 in topology.bonds():
                if atom1.element == elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == elem.hydrogen and atom1.element not in (elem.hydrogen, None):
                    transferMass = hydrogenMass-system.getParticleMass(atom2.index)
                    system.setParticleMass(atom2.index, hydrogenMass)
                    system.setParticleMass(atom1.index, system.getParticleMass(atom1.index)-transferMass)

        # Set periodic boundary conditions.
        boxVectors = topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            system.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
        elif nonbondedMethod not in [NoCutoff, CutoffNonPeriodic]:
            raise ValueError('Requested periodic boundary conditions for a Topology that does not specify periodic box dimensions')

        # Apply keyword-specified constraints.
        # QUESTION: Do we need the capability to enforce constraints specified via a keyword argument in addition to SMIRNOFF-specified constraints??
        if constraints != None:
            raise Exception("Constraints are not implemented yet.")

        # Allow the nonbonded method to be overridden
        # QUESTION: Do we want to make the choice of nonbonded method part of the forcefield spec?
        kwargs['nonbondedMethod'] = nonbondedMethod
        kwargs['nonbondedCutoff'] = nonbondedCutoff

        # Set user-specified charge method
        # QUESTION: If the charge method is part of the forcefield spec, do we really want users to be able to override it?
        kwargs['chargeMethod'] = chargeMethod

        # Add forces to the System
        for force in self._forces:
            if 'createForce' in dir(force):
                force.createForce(system, topology, **kwargs)

        # Add center-of-mass motion removal, if requested
        if removeCMMotion:
            system.addForce(openmm.CMMotionRemover())

        # Let force generators do postprocessing
        for force in self._forces:
            if 'postprocessSystem' in dir(force):
                force.postprocessSystem(system, topology, **kwargs)

        return system

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
                # QUESTION: Do we want to store these by OpenMM force type (e.g., `HarmonicBondForce`) or by SMIRNOFF tag name (`BondForce`)
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

# TODO: Move these collections to a utility file like ``openforcefield.utils`` or ``openforcefield.collections``?

import collections
class TransformedDict(collections.MutableMapping):
    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key

class ValenceDict(TransformedDict):
    """Enforce uniqueness in atom indices"""
    def __keytransform__(self, key):
        """Reverse tuple if first element is larger than last element."""
        # Ensure key is a tuple.
        key = tuple(key)
        # Reverse the key if the first element is bigger than the last.
        if key[0] > key[-1]:
            key = tuple(reversed(key))
        return key

class ImproperDict(TransformedDict):
    """Symmetrize improper torsions"""
    def __keytransform__(self,key):
        """Reorder tuple in numerical order except for element[1] which is the central atom; it retains its position."""
        # Ensure key is a tuple
        key = tuple(key)
        # Retrieve connected atoms
        connectedatoms = [key[0], key[2], key[3]]
        # Sort connected atoms
        connectedatoms.sort()
        # Re-store connected atoms
        key = tuple( [connectedatoms[0], key[1], connectedatoms[1], connectedatoms[2]])
        return(key)

#=============================================================================================
# Force generators
#=============================================================================================

## @private
class ForceGenerator(object):
    """Base class for force generators.
    """
    _TAGNAME = None # str of XML element name handled by this ForceGenerator
    _INFOTYPE = None # container class with type information that will be stored in self._types
    _OPENMMTYPE = None # OpenMM Force class (or None if no equivalent)

    def __init__(self, forcefield):
        self.ff = forcefield # the ForceField that this ForceGenerator is registered with
        self._types = list() # list of ForceType objects of type cls._INFOTYPE

    def getMatches(self, topology, **kwargs):
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

    @classmethod
    def parseElement(cls, tag, element, ff):
        # Find existing force generator or create new one.
        existing = [f for f in ff._forces if isinstance(f, cls)]
        if len(existing) == 0:
            generator = cls(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]

        # Register all SMIRNOFF definitions for all occurrences of its registered tag
        for section in element.findall(tag):
            generator.registerType(section, element)

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
        self.smirks = _validateSMIRKS(node.attrib['smirks'], node=node)
        if 'id' in node.attrib:
            self.pid = _extractQuantity(node, parent, 'id')

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
                self.distance = _extractQuantity(node, parent, 'distance') # Constraint with specified distance will be added by ConstraintGenerator
            else:
                self.distance = True # Constraint to equilibrium bond length will be added by HarmonicBondGenerator

    _TAGNAME = 'Constraint'
    _INFOTYPE = ConstraintType
    _OPENMMTYPE = None # don't create a corresponding OpenMM Force class

    def __init__(self, forcefield):
        super(ConstraintGenerator, self).__init__(forcefield)

    def createForce(self, system, topology, **kwargs):
        constraints = self.getMatches(topology)
        for (atoms, constraint) in constraints.items():
            # Update constrained atom pairs in topology
            topology.constrain_atom_pair(*atoms, constraint.distance)
            # If a distance is specified (constraint.distance != True), add the constraint here.
            # Otherwise, the equilibrium bond length will be used to constrain the atoms in HarmonicBondGenerator
            if constraint.distance is not True:
                system.addConstraint(*atoms, constraint.distance)

# QUESTION: Would we want all subclasses of ForceGenerator automatically register themselves?
parsers[ConstraintGenerator._TAGNAME] = ConstraintGenerator.parseElement

## @private
class BondGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<BondForce>`` tags"""

    class BondType(ForceType):
        """A SMIRNOFF bond type"""
        def __init__(self, node, parent):
            super(ConstraintType, self).__init__(node, parent) # Base class handles ``smirks`` and ``id`` fields

            # Determine if we are using fractional bond orders for this bond
            # First, check if this force uses fractional bond orders
            if 'fractional_bondorder' in parent.attrib:
                # If it does, see if this parameter line provides fractional bond order parameters
                if 'length_bondorder1' in node.attrib and 'k_bondorder1' in node.attrib:
                    # Store what interpolation scheme we're using
                    self.fractional_bondorder = parent.attrib['fractional_bondorder']
                    # Store bondorder1 and bondorder2 parameters
                    self.k = list()
                    self.length = list()
                    for ct in range(1,3):
                        self.length.append( _extractQuantity(node, parent, 'length_bondorder%s' % ct, unit_name = 'length_unit') )
                        self.k.append( _extractQuantity(node, parent, 'k_bondorder%s' % ct, unit_name = 'k_unit') )
                else:
                    self.fractional_bondorder = None
            else:
                self.fractional_bondorder = None

            # If no fractional bond orders, just get normal length and k
            if self.fractional_bondorder is None:
                self.length = _extractQuantity(node, parent, 'length')
                self.k = _extractQuantity(node, parent, 'k')

    _TAGNAME = 'BondForce' # SMIRNOFF tag name to process
    _INFOTYPE = BondType # class to hold force type info
    _OPENMMTYPE = openmm.HarmonicBondForce # OpenMM force class to create

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
            assert_bonded(atoms[0], atoms[1])

            # Compute equilibrium bond length and spring constant.
            if bond.fractional_bondorder is None:
                [k, length] = [bond.k, bond.length]
            else:
                # This bond uses partial bond orders
                # Make sure forcefield asks for fractional bond orders
                # TODO: Why shouldn't fractional bond orders just be automatically supported if used? This seems like a problem with setting this flag.
                if not self.ff._use_fractional_bondorder:
                    raise ValueError("Error: your forcefield file does not request to use fractional bond orders in its header, but a harmonic bond attempts to use them.")
                # Proceed to do interpolation
                order = topology.get_fractional_bond_order(*atoms)
                if bond.fractional_bondorder == 'interpolate-linear':
                    k = bond.k[0] + (bond.k[1]-bond.k[0])*(order-1.)
                    length = bond.length[0] + (bond.length[1]-bond.length[0])*(order-1.)
                else:
                    raise Exception("Partial bondorder treatment %s is not implemented." % bond.fractional_bondorder)

            # Handle constraints.
            if topology.atom_pair_is_constrained(*atoms):
                # Atom pair is constrained; we don't need to add a bond term.
                skipped_constrained_bonds += 1
                # Check if we need to add the constraint here to the equilibrium bond length.
                if topology.atom_pair_is_constrained(*atoms) is True:
                    # Mark that we have now assigned a specific constraint distance to this constraint.
                    topology.constrain_atom_pair(*atoms, length)
                    # Add the constraint to the System.
                    system.addConstraint(*particle_indices, length)
                continue

            # Add harmonic bond to HarmonicBondForce
            force.addBond(*particle_indices, length, k)

        logger.info('{} bonds added ({} skipped due to constraints)'.format(len(bonds) - skipped_constrained_bonds, skipped_constrained_bonds))

        # Check that no topological bonds are missing force parameters
        _check_for_missing_valence_terms('BondForce', topology, bonds.keys(), topology.bonds())

parsers[BondGenerator._TAGNAME] = BondGenerator.parseElement

#=============================================================================================

## @private
class AngleGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<AngleForce>`` tags"""

    class AngleType(ForceType):
        """A SMIRNOFF angle type."""
        def __init__(self, node, parent):
            super(AngleType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.angle = _extractQuantity(node, parent, 'angle')
            self.k = _extractQuantity(node, parent, 'k')
            if 'fractional_bondorder' in parent.attrib:
                self.fractional_bondorder = parent.attrib['fractional_bondorder']
            else:
                self.fractional_bondorder = None

    _TAGNAME = 'AngleForce' # SMIRNOFF tag name to process
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
                assert_bonded(topology, atoms[i], atoms[j])

            if topology.is_constrained(atoms[0], atoms[1]) and topology.is_constrained(atoms[1], atoms[2]) and topology.is_constrained(atoms[0], atoms[2]):
                # Angle is constrained; we don't need to add an angle term.
                skipped_constrained_angles += 1
                continue

            force.addAngle(*particle_indices, angle.angle, angle.k)

        logger.info('{} angles added ({} skipped due to constraints)'.format(len(angles) - skipped_constrained_angles, skipped_constrained_angles))

        # Check that no topological angles are missing force parameters
        _check_for_missing_valence_terms('AngleForce', topology, angles.keys(), topology.angles())

parsers[AngleForce._TAGNAME] = AngleGenerator.parseElement

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

            # Check that the SMIRKS pattern matches the type it's supposed to
            # TODO: Do we need this check?
            try:
                chemenv = ChemicalEnvironment(self.smirks)
                thistype = chemenv.getType()
                if thistype != 'Torsion':
                    raise Exception("Error: SMIRKS pattern %s (parameter %s) does not specify a %s torsion, but it is supposed to." % (self.smirks, self.pid, 'Proper'))
            except SMIRKSParsingError:
                print("Warning: Could not confirm whether smirks pattern %s is a valid %s torsion." % (self.smirks, self.torsiontype))

            if 'fractional_bondorder' in parent.attrib:
                self.fractional_bondorder = parent.attrib['fractional_bondorder']
            else:
                self.fractional_bondorder = None

            # Store parameters.
            index = 1
            while 'phase%d'%index in node.attrib:
                self.periodicity.append( int(_extractQuantity(node, parent, 'periodicity%d' % index)) )
                self.phase.append( _extractQuantity(node, parent, 'phase%d' % index, unit_name='phase_unit') )
                self.k.append( _extractQuantity(node, parent, 'k%d' % index, unit_name='k_unit') )
                # Optionally handle 'idivf', which divides the periodicity by the specified value
                if ('idivf%d' % index) in node.attrib:
                    idivf = _extractQuantity(node, parent, 'idivf%d' % index)
                    self.k[-1] /= float(idivf)
                index += 1

            # Check for errors, i.e. 'phase' instead of 'phase1'
            # TODO: Can we raise a more useful error if there is no ``id``?
            if len(self.phase)==0:
               raise Exception("Error: Torsion with id %s has no parseable phase entries." % self.pid)

    _TAGNAME = 'ProperTorsionForce' # SMIRNOFF tag name to process
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
                assert_bonded(topology, atoms[i], atoms[j])

            for (periodicity, phase, k) in zip(torsion.periodicity, torsion.phase, torsion.k):
                force.addTorsion(atom_indices[0], atom_indices[1], atom_indices[2], atom_indices[3], periodicity, phase, k)

        logger.info('{} torsions added'.format(len(torsions)))

        # Check that no topological torsions are missing force parameters
        _check_for_missing_valence_terms('ProperTorsionForce', topology, torsions.keys(), topology.torsions())

parsers[ProperTorsionGenerator._TAGNAME] = ProperTorsionGenerator.parseElement

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
                self.periodicity.append( int(_extractQuantity(node, parent, 'periodicity%d' % index)) )
                self.phase.append( _extractQuantity(node, parent, 'phase%d' % index, unit_name='phase_unit') )
                self.k.append( _extractQuantity(node, parent, 'k%d' % index, unit_name='k_unit') )
                # Optionally handle 'idivf', which divides the periodicity by the specified value
                if ('idivf%d' % index) in node.attrib:
                    idivf = _extractQuantity(node, parent, 'idivf%d' % index)
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
                assert_bonded(topology, atoms[i], atoms[j])

            # Impropers are applied to all six paths around the trefoil
            for (periodicity, phase, k) in zip(improper.periodicity, improper.phase, improper.k):
                # Permute non-central atoms
                others = [ atom_indices[0], atom_indices[2], atom_indices[3] ]
                for p in itertools.permutations( others ):
                    force.addTorsion(p[0], atom_indices[1], p[1], p[2], periodicity, phase, k)

        logger.info('{} impropers added, each applied in a six-fold trefoil' % (len(impropers)))

        # Check that no topological torsions are missing force parameters
        _check_for_missing_valence_terms('ImproperTorsionForce', topology, torsions.keys(), topology.impropers())

parsers[ImproperTorsionGenerator._TAGNAME] = ImproperTorsionGenerator.parseElement

## @private
class NonbondedGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<NonbondedGenerator>`` tags"""

    SCALETOL = 1e-5

    class NonbondedType(ForceType):
        """A SMIRNOFF Lennard-Jones type."""
        def __init__(self, node, parent):
            # NOTE: Currently we support radius definition via 'sigma' or 'rmin_half'.
            super(NonbondedType, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields

            # TODO: Do we want to add an attribute ``type='LennardJones'``?
            # Make sure we don't have BOTH rmin_half AND sigma
            try:
                a = _extractQuantity(node, parent, 'sigma')
                a = _extractQuantity(node, parent, 'rmin_half')
                raise Exception("Error: BOTH sigma and rmin_half cannot be specified simultaneously in the .ffxml file.")
            except:
                pass

            # Handle Lennard-Jones sigma
            try:
                self.sigma = _extractQuantity(node, parent, 'sigma')
            #Handle rmin_half, AMBER-style
            except:
                rmin_half = _extractQuantity(node, parent, 'rmin_half', unit_name='sigma_unit')
                self.sigma = 2.*rmin_half/(2.**(1./6.))
            self.epsilon = _extractQuantity(node, parent, 'epsilon')

    _TAGNAME = 'NonbondedForce' # SMIRNOFF tag name to process
    _INFOTYPE = NonbondedType # info type to store
    _OPENMMTYPE = openmm.NonbondedForce # OpenMM force class to create

    def __init__(self, forcefield, coulomb14scale, lj14scale):
        super(NonbondedForceGenerator, self).__init__(forcefield)
        self.coulomb14scale = coulomb14scale
        self.lj14scale = lj14scale


    # TODO: Handle the case where multiple <NonbondedForce> tags are found
    # if abs(generator.coulomb14scale - float(element.attrib['coulomb14scale'])) > NonbondedGenerator.SCALETOL or \
    #                 abs(generator.lj14scale - float(element.attrib['lj14scale'])) > NonbondedGenerator.SCALETOL:
    #            raise ValueError('Found multiple NonbondedForce tags with different 1-4 scales')
    #    for atom in element.findall('Atom'):
    #        generator.registerAtom(atom, element)

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
        _check_for_missing_valence_terms('NonbondedForce Lennard-Jones parameters', topology, atoms.keys(), topology.atoms())

        # Set the partial charges based on reference molecules.
        for reference_molecule in topology._reference_molecules:
            atom_mappings = topology._reference_to_topology_atom_mappings[reference_molecule]
            # Retrieve charges from reference molecule, stored by atom index
            charge_by_atom = {}
            for atom in reference_molecule.GetAtoms():
                charge_by_atom[atom.GetIdx()] = atom.GetPartialCharge()

            # Loop over mappings and copy NB parameters from reference molecule
            # to other instances of the molecule
            for atom_mapping in atom_mappings:
                for (atom_index, map_atom_index) in atom_mapping.items():
                    # Retrieve NB params for reference atom (charge not set yet)
                    charge, sigma, epsilon = force.getParticleParameters(map_atom_index)
                    # Look up the charge on the atom in the reference molecule
                    charge = charge_by_atom[atom_index]*unit.elementary_charge

                    # Set parameters for equivalent atom in other instance of
                    # this molecule
                    force.setParticleParameters(map_atom_index, charge, sigma, epsilon)
        # TODO: Should we check that there are no missing charges?

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

parsers[NonbondedGenerator._TAGNAME] = NonbondedGenerator.parseElement

## @private
class BondChargeCorrectionGenerator(ForceGenerator):
    """Handle SMIRNOFF ``<BondChargeCorrection>`` tags"""

    class BondChargeCorrectionType(ForceType):
        """A SMIRNOFF bond charge correction type."""
        def __init__(self, node, parent):
            super(BondChargeCorrectionGenerator, self).__init__(node, parent) # base class handles ``smirks`` and ``id`` fields
            self.increment = _extractQuantity(node, parent, 'increment')
            # If no units are specified, assume elementary charge
            if type(self.increment) == float:
                self.increment *= unit.elementary_charge

    _TAGNAME = 'BondChargeCorrection' # SMIRNOFF tag name to process
    _INFOTYPE = BondChargeCorrectionType # info type to store
    _OPENMMTYPE = openmm.NonbondedForce # OpenMM force class to create or utilize

    def __init__(self, forcefield):
        super(NonbondedForceGenerator, self).__init__(forcefield)

    #if element.attrib['method'] != existing[0]._initialChargeMethod:
    #raise Exception("Existing BondChargeCorrectionGenerator uses initial charging method '%s' while new BondChargeCorrectionGenerator requests '%s'" % (existing[0]._initialChargeMethod, element.attrib['method']))

    def createForce(self, system, topology, **args):
        # No forces are created by this generator.
        pass

    # TODO: Move chargeMethod to SMIRNOFF spec
    def postprocessSystem(self, system, topology, chargeMethod=None, **args):
        if chargeMethod != 'BCC':
            # Only apply charge corrections if chargeMethod is 'BCC'
            return

        # Iterate over all defined bond charge corrections, allowing later matches to override earlier ones
        bonds = self.getMatches(topology)

        # Apply bond charge increments
        for force in system.getForces():
            if force.__class__.__name__ in ['NonbondedForce']:
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

parsers[BondChargeCorrectionGenerator._TAGNAME] = BondChargeCorrectionGenerator.parseElement

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
                    value = _extractQuantity(node, parent, name)
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

    # TODO: Generalize this
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

parsers[GBSAForceGenerator._TAGNAME] = GBSAForceGenerator.parseElement
