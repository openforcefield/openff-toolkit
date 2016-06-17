#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
atomtyper.py

Atom type assignment engine using SMARTS strings.

Authors
-------
John Chodera <john.chodera@choderalab.org>, Memorial Sloan Kettering Cancer Center and University of California, Berkeley

The AtomTyper class is based on 'patty' from Pat Walters, Vertex Pharmaceuticals.

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

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
            msg = "Atom not assigned: molecule %s : atom index %6d name %8s element %8s" % (self.molecule.GetTitle(), self.atom.GetIdx(), self.atom.GetName(), OEGetAtomicSymbol(self.atom.GetAtomicNum()))
            msg += '\n'
            for atom in self.molecule.GetAtoms():
                msg += 'atom %8d : name %8s element %8s' % (atom.GetIdx(), atom.GetName(), OEGetAtomicSymbol(self.atom.GetAtomicNum()))
                if atom == self.atom:
                    msg += '  ***'
                msg += '\n'

            return msg

    def __init__(self, typelist, tagname, replacements=None):
        """"
        Create an atom typer instance.

        ARGUMENTS

        typelist : str
            If specified, will read types from list with each element [smarts, typename]
        tagname : str
            Tag name
        replacements : list of [smarts, shortname]
            Substitution/replacement bindings.

        """

        self.pattyTag = OEGetTag(tagname)

        # Create bindings list.
        bindings = list()
        if replacements is not None:
            for [smarts,shortname] in replacements:
                bindings.append( (shortname, smarts) )

        # Create table of search objects.
        self.smartsList = []
        for [smarts, typename] in typelist:
            # Perform binding replacements
            smarts = OESmartsLexReplace(smarts, bindings)
            # Create SMARTS search
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

        Parameters
        ----------
        filename : str
            The name of the file to be read

        Returns
        -------
        typelist : list of tuples
            Typelist[i] is element i of the typelist in format [smarts, typename]

        """
        if filename is None:
            return None

        if not os.path.exists(filename):
            raise Exception("File '%s' not found." % filename)

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
            # Process if we have enough tokens
            if len(tokens) >= 2:
                smarts = tokens[0]
                typename = ' '.join(tokens[1:])
                typelist.append([smarts,typename])
        ifs.close()

        return typelist
