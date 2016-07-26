#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
forcefield_labeler.py

Class structure mirroring forcefield.py but for simply determining what parameter numbers would be assigned to specified SMIRKS.

AUTHORS

John D. Chodera <john.chodera@choderalab.org>
David L. Mobley <dmobley@mobleylab.org>

Baseed on simtk.openmm.app.forcefield written by Peter Eastman.

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

import lxml.etree as etree

import os
import math
import copy
import re
import numpy
import random

import openeye.oechem
import openeye.oeomega
import openeye.oequacpac

from openeye import oechem

from simtk import openmm, unit

import time

import networkx

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

def _getSMIRKSMatches( oemol, smirks):
        """Find all sets of atoms in the provided oemol that match the provided SMIRKS strings.

            Parameters
        ----------
        smirks : str
            SMIRKS string with tagged atoms.
            If there are N tagged atoms numbered 1..N, the resulting matches will be N-tuples of atoms that match the corresponding tagged atoms.

        Returns
        -------
        matches : list of tuples of atoms numbers
            matches[index] is an N-tuple of atom numbers from the oemol
            Matches are returned in no guaranteed order.
        """

        # Set up query.
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, smirks):
            raise Exception("Error parsing SMIRKS '%s'" % smirks)

        # Perform matching on each unique molecule, unrolling the matches to all matching copies of tha tmolecule in the Topology object.
        matches = list()

        # We require non-unique matches, i.e. all matches
        unique = False
        ss = oechem.OESubSearch(qmol)
        matches = []
        for match in ss.Match( oemol, unique):
            # Compile list of atom indices that match the pattern tags
            atom_indices = dict()
            for ma in match.GetAtoms():
                if ma.pattern.GetMapIdx() != 0:
                    atom_indices[ma.pattern.GetMapIdx()-1] = ma.target.GetIdx()
            # Compress into list
            atom_indices = [ atom_indices[index] for index in range(len(atom_indices)) ]
            # Store
            matches.append( tuple(atom_indices) )

        return matches

#=============================================================================================
# FORCEFIELD LABELER
#=============================================================================================

# Map of functions to parse elements of XML file
parsers = {}

class ForceField_labeler(object):
    """A ForceField_labeler holds parameters and SMIRKS from a SMIRFF XML file and can be used to label OEMols with what force terms would be applied."""


    # Adapted from forcefield.py
    def __init__(self, *files):
        """Load one or more XML parameter definition files and create a SMIRFF ForceField object based on them.

        Parameters
        ----------
        files : list
            A list of XML files defining the SMIRFF force field.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in force fields), or an open file-like object with a read() method from which the forcefield XML data can be loaded.

        """
        self._forces = [] 
        self.loadFile(files)


    def loadFile(self, files):
        """Load a SMIRFF XML file and add the parameter types, SMIRKS, and parameter IDs to this ForceField_labeler.

        Parameters
        ----------
        files : string or file or tuple
            An XML file or tuple of XML files containing SMIRFF force field definitions.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in force fields), or an open file-like object with a read() method from which the forcefield XML data can be loaded.
        """

        # Ensure that we are working with a tuple of files.
        if not isinstance(files, tuple):
            files = (files,)

        # Load in all XML trees.
        trees = list()
        for file in files:
            try:
                # this handles either filenames or open file-like objects
                tree = etree.parse(file)
            except IOError:
                tree = etree.parse(os.path.join(os.path.dirname(__file__), 'data', file))
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


        # Load force definitions
        for tree in trees:
            for child in tree.getroot():
                # Ignore anything we can't parse
                if child.tag in parsers:
                    parsers[child.tag](child, self)


    def getGenerators(self):
        """Get the list of all registered generators."""
        return self._forces

    def registerGenerator(self, generator):
        """Register a new generator."""
        self._forces.append(generator)

    def labelMolecules( self, oemols, verbose = False):
        """Return labels for a set of OEMols corresponding to parameters from this force field.

        Parameters
        ----------
        oemols : list of OEMols
            The OpenEye OEChem OEMol objects as a list; these will be labeled. Should include all atoms with the correct ordering atom atom numbers will be returned along with corresponding labeling.
        verbose : bool
            If True, verbose output will be printed

        Returns
        -------
        molecule_labels : list
            list of labels for molecules. Each entry in the list corresponds to
            one molecule from the provided list of oemols and is a dictionary
            keyed by force type, i.e. molecule_labels[0]['HarmonicBondForce'] 
            gives details for the harmonic bond parameters for the first 
            molecule. Each element is a list of the form [ ( ( atom1, ..., 
            atomN), parameter_id), ... ]
        
        """

        molecule_labels = []

        # Loop over molecules and label
        for idx,mol in enumerate(oemols):
            molecule_labels.append({})
            for force in self._forces:
                # Initialize dictionary storage for this force type
                if isinstance(force, HarmonicBondGenerator):
                    forcelabel = 'HarmonicBondForce'
                # TO DO: Add other force types here
                else:
                    continue
                    if verbose: print("Encountered unimplemented force; skipping.")
   
                # Grab force terms of this type for this molecule and store
                molecule_labels[idx][forcelabel] = force.parseForce( mol, verbose=verbose )
        return molecule_labels

#=============================================================================================
# Utility functions
#=============================================================================================

def _validateSMIRKS(smirks, node=None):
    """Validate the specified SMIRKS string.

    Parameters
    ----------
    smirks : str
       The SMIRKS string to be validated
    node : xml.etree.ElementTree.Element
       Node of etree with 'sourceline' attribute.

    """
    qmol = oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smirks):
        if (node is not None) and ('sourceline' in node.attrib):
            raise Exception("Line %s: Error parsing SMIRKS '%s'" % (node.attrib['sourceline'], node.attrib['smirks']))
        else:
            raise Exception("Error parsing SMIRKS '%s'" % (node.attrib['smirks']))

    return smirks

def _extractQuantity(node, parent, name, unit_name=None):
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

    """
    if name not in node.attrib:
        if 'sourceline' in node.attrib:
            raise Exception("Line %d : Expected XML attribute '%s' not found" % (node.attrib['sourceline'], name))
        else:
            raise Exception("Expected XML attribute '%s' not found" % (name))

    string_names = ['parent_id', 'id']
    # Handle case where this is a normal quantity
    if name not in string_names:
        quantity = float(node.attrib[name])
    # Handle case where it is a label or string
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


#=============================================================================================
# Parsers
#=============================================================================================

# Adapted from forcefield.py

class HarmonicBondGenerator(object):
    """A HarmonicBondGenerator parses a HarmonicBondForce."""

    class BondType(object):
        """A SMIRFF bond type."""
        def __init__(self, node, parent):
            self.smirks = _validateSMIRKS(node.attrib['smirks'], node=node)
            self.pid = _extractQuantity(node, parent, 'id')

    def __init__(self, forcefield):
        self.ff = forcefield
        self._bondtypes = list()

    def registerBond(self, node, parent):
        """Register a SMIRFF bondtype definition."""
        bond = HarmonicBondGenerator.BondType(node, parent)
        self._bondtypes.append(bond)

    @staticmethod
    def parseElement(element, ff):
        # Find existing force generator or create new one.
        existing = [f for f in ff._forces if isinstance(f, HarmonicBondGenerator)]
        if len(existing) == 0:
            generator = HarmonicBondGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]

        # Register all SMIRFF bond definitions.
        for bond in element.findall('Bond'):
            generator.registerBond(bond, element)


    def parseForce(self, oemol, verbose=False, **kwargs):
        """Take a provided OEMol and parse HarmonicBondForce terms for this molecule.

        Parameters
        ----------
            oemol : OEChem OEMol object for molecule to be examined

        Returns
        ---------
            force_terms: list
                Returns a list of tuples, [ ((atom id 1, ... atom id N), parameter id) , (....), ... ] for all forces of this type which would be applied.
        """

        # Iterate over all defined bond SMIRKS, allowing later matches to override earlier ones.
        bonds = ValenceDict()
        for bond in self._bondtypes:
            for atom_indices in _getSMIRKSMatches( oemol, bond.smirks ): 
                bonds[atom_indices] = bond

        if verbose:
            print('')
            print('HarmonicBondGenerator:')
            print('')
            for bond in self._bondtypes:
                print('%64s : %8d matches' % (bond.smirks, len(_getSMIRKSMatches(oemol, bond.smirks))))
            print('')

        # Add all bonds to the output list
        force_terms = []
        for (atom_indices, bond) in bonds.items():
            force_terms.append( ((atom_indices[0], atom_indices[1]), bond.pid) )

        return force_terms


parsers["HarmonicBondForce"] = HarmonicBondGenerator.parseElement
