#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
smarty.py

Example illustrating a scheme to create and destroy atom types automatically using SMARTS.

AUTHORS

John Chodera <jchodera@berkeley.edu>, University of California, Berkeley

The AtomTyper class is based on 'patty' by Pat Walters, Vertex Pharmaceuticals.

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

from optparse import OptionParser # For parsing of command line arguments

import os
import math
import copy
import re
import numpy
import random

import openeye.oechem
import openeye.oeomega
import openeye.oequacpac

from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *

import time

#=============================================================================================
# ATOM TYPER
#=============================================================================================

class AtomTyper(object):
    """
    Atom typer based on SMARTS-defined atom types.

    Based on 'Patty' implementation by Pat Walters.

    """

    class TypingException(Exception):
        """
        Atom typing exception.

        """
        def __init__(self, molecule, atom):
            self.molecule = molecule
            self.atom = atom

        def __str__(self):
            return "Atom not assigned: %6d %8s" % (self.atom.GetIdx(), OEGetAtomicSymbol(self.atom.GetAtomicNum()))

    def __init__(self, typelist, tagname):
        """"
        Create an atom typer instance.

        ARGUMENTS

        typelist (string) - if specified, will read types from list with each element [smarts, typename]
        tagname (string) - tag name

        """

        self.pattyTag = OEGetTag(tagname)

        # Create table of search objects.
        self.smartsList = []
        for [smarts, typename] in typelist:
            pat = OESubSearch()
            pat.Init(smarts)
            pat.SetMaxMatches(0)
            self.smartsList.append([pat,typename,smarts])

        return

    def dump(self):
        for pat,type,smarts in self.smartsList:
            print pat,type,smarts
        return

    def assignTypes(self,mol):
        # Assign null types.
        for atom in mol.GetAtoms():
            atom.SetStringData(self.pattyTag, "")

        # Assign atom types using rules.
        OEAssignAromaticFlags(mol)
        for pat,type,smarts in self.smartsList:
            for matchbase in pat.Match(mol):
                for matchpair in matchbase.GetAtoms():
                    matchpair.target.SetStringData(self.pattyTag,type)

        # Check if any atoms remain unassigned.
        for atom in mol.GetAtoms():
            if atom.GetStringData(self.pattyTag)=="":
                raise AtomTyper.TypingException(mol, atom)
        return

    def debugTypes(self,mol):
        for atom in mol.GetAtoms():
            print "%6d %8s %8s" % (atom.GetIdx(),OEGetAtomicSymbol(atom.GetAtomicNum()),atom.GetStringData(self.pattyTag))
        return

    def getTypeList(self,mol):
        typeList = []
        for atom in mol.GetAtoms():
            typeList.append(atom.GetStringData(self.pattyTag))
        return typeList

    @classmethod
    def read_typelist(cls, filename):
        """
        Read an atomtype or decorator list from a file.

        ARGUMENTS

        filename (string) - the name of the file to be read

        RETURNS

        typelist (list) - typelist[i] is element i of the typelist in format [smarts, typename]

        """
        typelist = list()
        ifs = open(filename)
        lines = ifs.readlines()
        for line in lines:
            # Strip trailing comments
            index = line.find('%')
            if index != -1:
                line = line[0:index]
            # Split into tokens.
            tokens = string.split(line)
            if len(tokens) == 2:
                [smarts,typename] = tokens
                typelist.append([smarts,typename])
        ifs.close()

        return typelist

#=============================================================================================
# ATOMTYPE SAMPLER
#=============================================================================================

class AtomTypeSampler(object):
    """
    Atom type sampler.

    """
    def __init__(self, basetypes_filename, decorators_filename, molecules):
        """
        Initialize an atom type sampler.

        ARGUMENTS

        basetypes_filename - file defining base atom types (which cannot be destroyed)
        decorators_filename - file containing decorators that can be added to existing types to generate subtypes
        molecules - list of molecules for typing

        NOTES

        This is just a proof of concept.  No scoring of molecular properties is performed.

        """

        # Define internal typing tag.
        self.typetag = 'atomtype'

        # Read atomtypes and decorators.
        self.atomtypes = AtomTyper.read_typelist(basetypes_filename)
        self.decorators = AtomTyper.read_typelist(decorators_filename)

        # Store a deep copy of the molecules since they will be annotated
        self.molecules = copy.deepcopy(molecules)

        # Type all molecules with current typelist to ensure that basetypes are sufficient.
        self.type_molecules(self.atomtypes, self.molecules)

        # Compute atomtype statistics on molecules.
        [atom_typecounts, molecule_typecounts] = self.compute_type_statistics(self.atomtypes, self.molecules)
        self.show_type_statistics(self.atomtypes, atom_typecounts, molecule_typecounts)

        return

    def sample_atomtypes(self, verbose=True):
        """
        Perform one step of atom type sampling.

        """
        # Copy current atomtypes for proposal.
        proposed_atomtypes = copy.deepcopy(self.atomtypes)
        proposed_molecules = copy.deepcopy(self.molecules)
        natomtypes = len(proposed_atomtypes)
        ndecorators = len(self.decorators)

        if random.random() < 0.5:
            # Pick an atom type to destroy.
            atomtype_index = random.randint(0, natomtypes-1)
            (atomtype, typename) = proposed_atomtypes[atomtype_index]
            if verbose: print("Attempting to destroy atom type %s : %s..." % (atomtype, typename))
            # Delete the atomtype.
            proposed_atomtypes.remove([atomtype, typename])
            # Try to type all molecules.
            try:
                self.type_molecules(proposed_atomtypes, proposed_molecules)
                # Accept if typing completed.
                self.atomtypes = proposed_atomtypes
                self.molecules = proposed_molecules
                return True
            except Exception as e:
                print e
                # Reject since typing failed.
                if verbose: print("Typing failed; rejecting.")
                return False
        else:
            # Pick an atomtype to subtype.
            atomtype_index = random.randint(0, natomtypes-1)
            # Pick a decorator to add.
            decorator_index = random.randint(0, ndecorators-1)
            # Create new atomtype to insert by appending decorator with 'and' operator.
            (atomtype, atomtype_typename) = self.atomtypes[atomtype_index]
            (decorator, decorator_typename) = self.decorators[decorator_index]
            result = re.match('\[(.+)\]', atomtype)
            proposed_atomtype = '[' + result.groups(1)[0] + '&' + decorator + ']'
            proposed_typename = atomtype_typename + ' ' + decorator_typename
            print("Attempting to create new subtype: '%s' (%s) + '%s' (%s) -> '%s' (%s)" % (atomtype, atomtype_typename, decorator, decorator_typename, proposed_atomtype, proposed_typename))
            # Check if proposed atomtype is already in set.
            existing_atomtypes = set()
            for (a, b) in self.atomtypes:
                existing_atomtypes.add(a)
            if proposed_atomtype in existing_atomtypes:
                if verbose: print("Atom type already exists; rejecting to avoid duplication.")
                return False
            # Insert atomtype immediately after.
            proposed_atomtypes.insert(atomtype_index+1, [proposed_atomtype, proposed_typename])
            print(proposed_atomtypes)
            # Try to type all molecules.
            try:
                # Type molecules.
                self.type_molecules(proposed_atomtypes, proposed_molecules)
                # Compute updated statistics.
                [proposed_atom_typecounts, proposed_molecule_typecounts] = self.compute_type_statistics(proposed_atomtypes, proposed_molecules)
                # Reject if new type is unused.
                if (proposed_atom_typecounts[proposed_typename] == 0):
                    # Reject because new type is unused in dataset.
                    if verbose: print("Atom type '%s' (%s) unused in dataset; rejecting." % (proposed_atomtype, proposed_typename))
                    return False
                # Reject if parent type is now unused.
                if (proposed_atom_typecounts[atomtype_typename] == 0):
                    # Reject because new type is unused in dataset.
                    if verbose: print("Parent type '%s' (%s) now unused in dataset; rejecting." % (atomtype, atomtype_typename))
                    return False
                # Accept.
                self.atomtypes = proposed_atomtypes
                self.molecules = proposed_molecules
                return True
            except Exception as e:
                print("Exception: %s" % str(e))
                # Reject since typing failed.
                if verbose: print("Typing failed for one or more molecules using proposed atomtypes; rejecting.")
                return False

        return

    def type_molecules(self, typelist, molecules):
        """
        Type all molecules with the specified typelist.

        """
        # Create an atom typer.
        atomtyper = AtomTyper(typelist, self.typetag)

        # Type molecules.
        for molecule in molecules:
            atomtyper.assignTypes(molecule)

        return

    def compute_type_statistics(self, typelist, molecules):
        """
        Compute statistics for numnber of molecules assigned each type.

        ARGUMENTS

        typelist
        molecules

        RETURNS

        atom_typecounts (dict) - counts of number of atoms containing each atomtype
        molecule_typecounds (dict) - counts of number of molecules containing each atom type

        """
        # Zero type counts by atom and molecule.
        atom_typecounts = dict()
        molecule_typecounts = dict()
        for [smarts, typename] in typelist:
            atom_typecounts[typename] = 0
            molecule_typecounts[typename] = 0

        # Count number of atoms with each type.
        for molecule in molecules:
            types_in_this_molecule = set()
            for atom in molecule.GetAtoms():
                atomtype = atom.GetStringData(self.typetag)
                types_in_this_molecule.add(atomtype)
                atom_typecounts[atomtype] += 1
            for atomtype in types_in_this_molecule:
                molecule_typecounts[atomtype] += 1

        return (atom_typecounts, molecule_typecounts)

    def show_type_statistics(self, typelist, atom_typecounts, molecule_typecounts):
        """
        Print atom type statistics.

        """
        index = 1
        natoms = 0
        #print "%5s   %10s %10s   %48s %48s" % ('index', 'atoms', 'molecules', 'type name', 'smarts')
        print "%5s   %10s %10s   %48s %48s" % ('INDEX', 'ATOMS', 'MOLECULES', 'TYPE NAME', 'SMARTS')
        for [smarts, typename] in typelist:
            print "%5d : %10d %10d | %48s %48s" % (index, atom_typecounts[typename], molecule_typecounts[typename], typename, smarts)
            natoms += atom_typecounts[typename]
            index += 1

        nmolecules = len(self.molecules)
        print "%5s   %10d %10d" % ('TOTAL', natoms, nmolecules)
        return

    def run(self, niterations, verbose=False):
        """
        Run atomtype sampler for the specified number of iterations.

        Parameters
        ----------
        niterations : int
            The specified number of iterations
        verbose : bool
            If True, print debug info
        """

        for iteration in range(niterations):
            accepted = self.sample_atomtypes(verbose=verbose)
            if verbose:
                print("Iteration %d / %d: %s" % (iteration, niterations, str(accepted)))

                # Compute atomtype statistics on molecules.
                [atom_typecounts, molecule_typecounts] = atomtype_sampler.compute_type_statistics(atomtype_sampler.atomtypes, atomtype_sampler.molecules)
                atomtype_sampler.show_type_statistics(atomtype_sampler.atomtypes, atom_typecounts, molecule_typecounts)
