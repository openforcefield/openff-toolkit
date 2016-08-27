#!/usr/bin/env python

#==============================================================================
# MODULE DOCSTRING
#==============================================================================

"""
environment.py

Classes defining a chemical environment for atoms and how they are connected
using networkx graph objects to organize and make changes to the structure.
Output will be in the form of SMARTS and SMIRKS.

AUTHORS

Caitlin Bannan <bannanc@uci.edu>, Mobley Lab, University of California Irvine,
with contributions from John Chodera, Memorial Sloan Kettering Cancer Center
and David Mobley, UC Irvine.

"""
#==============================================================================
# GLOBAL IMPORTS
#==============================================================================

import networkx as nx
import re
import copy

import numpy as np
from numpy import random

import openeye.oechem
from openeye.oechem import *

# -------------
# to do list:
# TODO: should there be print statments in here, should they have a verbose option?
# -------------

class ChemicalEnvironment(object):
    """Chemical environment abstract base class that matches an atom, bond, angle, etc.
    """
    class Atom(object):
        """Atom representation, which may have some ORtypes and ANDtypes properties.

        Properties
        -----------
        ORtypes : list of tuples in the form (base, [list of decorators])
            where bases and decorators are both strings
            The descriptor types that will be combined with logical OR
        ANDtypes : list of string
            The descriptor types  that will be combined with logical AND
        """
        def __init__(self, ORtypes = None, ANDtypes = None, index = None):
            """Initialize an Atom object with optional descriptors.

            Parameters
            -----------
            ORtypes: list of tuples for ORbases and ORdecorators, optional, default = None
                in the form (base, [list of decorators])
            ANDtypes: list of str, optional, default = None
                strings that will be AND'd together in a SMARTS
            index : int, optional, default=None
                If not None, the specified index will be attached as a SMIRKS index (e.g. '[#6:1]')
            """
            # dictionary of ORbases and ORdecorators
            if ORtypes ==  None:
                self.ORtypes = list()
            else:
                self.ORtypes = copy.deepcopy(ORtypes)

            # Set of strings that will be AND'd to the the end
            if ANDtypes == None:
                self.ANDtypes = list()
            else:
                self.ANDtypes = list(ANDtypes)

            self.index = index
            self._atom = True

        def asSMARTS(self):
            """Return the atom representation as SMARTS.

            Returns
            --------
            smarts : str
            The SMARTS string for this atom
            """

            smarts = '['

            # Add the OR'd features
            if self.ORtypes:
                ORList = list()
                for (base, ORdecorators) in self.ORtypes:
                    if base[0] == '$':
                        if ORdecorators:
                            OR = base+'&'+''.join(ORdecorators)
                        else:
                            OR = base
                    else: # base doesn't start with $
                        OR = base+''.join(ORdecorators)
                    ORList.append(OR)
                smarts += ','.join(ORList)
            else:
                smarts += '*'

            if len(self.ANDtypes) > 0:
                smarts += ';' + ';'.join(self.ANDtypes)

            return smarts + ']'

        def asSMIRKS(self):
            """Return the atom representation as SMIRKS.

            Returns
            --------
            smirks : str
            The SMIRKS string for this atom
            """
            smirks = self.asSMARTS()

            # No index specified so SMIRKS = SMARTS
            if self.index == None:
                return smirks

            # Add label to the end of SMARTS
            else:
                return smirks[:-1] + ':' + str(self.index) + smirks[-1]

        def addORtype(self, ORbase, ORdecorators):
            """
            Adds ORtype to the set for this atom.

            Parameters
            --------
            ORbase: string, such as '#6'
            ORdecorators: list of strings, such as ['X4','+0']
            """
            self.ORtypes.append((ORbase, ORdecorators))

        def addANDtype(self, ANDtype):
            """
            Adds ANDtype to the set for this atom.

            Parameters
            --------
            ANDtype: string
                added to the list of ANDtypes for this atom
            """
            self.ANDtypes.append(ANDtype)

        def getORtypes(self):
            """
            returns a copy of the dictionary of ORtypes for this atom
            """
            return copy.deepcopy(self.ORtypes)

        def setORtypes(self, newORtypes):
            """
            sets new ORtypes for this atom

            Parameters
            ----------
            newORtypes: list of tuples in the form (base, [ORdecorators])
                for example: ('#6', ['X4','H0','+0']) --> '#6X4H0+0'
            """
            if newORtypes == None:
                self.ORtypes = list()
            else:
                self.ORtypes = newORtypes

        def getANDtypes(self):
            """
            returns a copy of the list of ANDtypes for this atom
            """
            return list(copy.deepcopy(self.ANDtypes))

        def setANDtypes(self, newANDtypes):
            """
            sets new ANDtypes for this atom

            Parameters
            ----------
            newANDtypes: list of strings
                strings that will be AND'd together in a SMARTS
            """
            if newANDtypes == None:
                self.ANDtypes = list()
            else:
                self.ANDtypes = list(newANDtypes)

    class Bond(Atom):
        """Bond representation, which may have ORtype and ANDtype descriptors.
        Properties
        -----------
        ORtypes : dictionary of bases and ORdecorators in form {base: [list of decorators]}
            The ORtype types that will be combined with logical OR
        ANDtypes : list of string
            The ANDtypes that will be combined with logical AND

        """
        # Implementation identical to atoms apart from what is put in the asSMARTS/asSMIRKS strings

        def __init__(self, ORtypes = None, ANDtypes = None):
            """
            Parameters
            -----------
            ORtypes: list of tuples, optional, default = None
                tuples have form (base, [ORdecorators])
                bond descriptors that will be OR'd together in a SMARTS
            ANDtypes: list of str, optional, default = None
                strings that will be AND'd together in a SMARTS
            index: integer, default = None
                This is for book keeping inside environments and will not be shown in SMARTS or SMIRKS
                example: bond1 in a Bond is the bond between atom1 and atom2
            """
            super(ChemicalEnvironment.Bond,self).__init__(ORtypes, ANDtypes, None)
            self._atom = False
            return

        def asSMARTS(self):
            """Return the atom representation as SMARTS.

            Returns
            --------
            smarts : str
                The SMARTS string for just this atom
            """
            if self.ORtypes:
                smarts = ''
                for (ORbase, ORdecorators) in self.ORtypes:
                    smarts += ORbase+''.join(ORdecorators)
            else:
                smarts = '~'

            if len(self.ANDtypes) > 0:
                smarts += ';' + ';'.join(self.ANDtypes)

            return smarts

        def asSMIRKS(self):
            """
            Returns
            --------
            smarts : str
                The SMIRKS string for just this bond
            """
            #the same as asSMARTS()
            #    for consistency asSMARTS() or asSMIRKS() can be called
            #    for all environment objects
            return self.asSMARTS()

        def getOrder(self):
            """
            Returns a float for the order of this bond
            for multiple ORtypes or ~ it returns the minimum possible order
            """
            orderDict = {'~':1.,
                    '-':1., ':': 1.5, '=':2., '#':3.,
                    '!-':1.5, '!:':1., '!=':1., '!#':1.}
            orderList = [orderDict[base] for (base, decor) in self.ORtypes]
            return min(orderList)

    def __init__(self, smirks = None, label = None):
        """Initialize a chemical environment abstract base class.

        smirks = string, optional
            if smirks is not None, a chemical environment is built
            from the provided SMIRKS string
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        """
        # Create an empty graph which will store Atom objects.
        self._graph = nx.Graph()
        self.label = label

        if smirks is not None:
            # check SMIRKS is parseable
            #mol = OEQMol()
            #if not OEParseSmarts(mol, smirks):
            #    raise Exception("Provides SMIRKS: %s was not parseable" % smirks)

            atoms = dict() # store created atom
            idx = 1 # current atom being created
            store = list() # to store indices while branching
            bondingTo = idx # which atom are we going to bond to

            start = smirks.find('[')
            if start != 0:
                raise Exception("Provided SMIRKS: %s should begin with '[' instead of %s" % (smirks, smirks[0]))
            end = smirks.find(']')

            atom = smirks[start+1:end]
            OR, AND, index = self._getAtomInfo(atom)
            leftover = smirks[end+1:]

            new_atom = self.addAtom(None, newORtypes = OR, newANDtypes = AND, newAtomIndex = index)
            atoms[idx] = new_atom

            while len(leftover) > 0:
                idx += 1

                # Check for branching
                if leftover[0] == ')':
                    bondingTo = store.pop()
                    leftover = leftover[1:]
                if leftover[0] == '(':
                    store.append(bondingTo)
                    leftover = leftover[1:]

                # find beginning and end of atom
                start = leftover.find('[')
                end = leftover.find(']')

                # Get bond and atom info
                bOR, bAND = self._getBondInfo(leftover[:start])
                aOR, aAND, index = self._getAtomInfo(leftover[start+1:end])

                # create new atom
                new_atom = self.addAtom(atoms[bondingTo], bOR, bAND, aOR, aAND, index)

                # update state
                atoms[idx] = new_atom
                bondingTo = idx
                leftover = leftover[end+1:]


    def _getAtomInfo(self, atom):
        """
        given atom string, returns ORtypes, ANDtypes, and index
        """
        # Find atom index
        colon = atom.find(':')
        if colon == -1:
            index = None
        else:
            index = int(atom[colon+1:])
            atom = atom[:colon]

        split = atom.split(';')

        # Get ANDtypes
        ANDtypes = split[1:]
        if len(ANDtypes) == 0:
            ANDtypes = None

        # Get ORtypes
        ORList = split[0].split(',')
        ORtypes = list()
        # Separate ORtypes into bases and decorators
        for OR in ORList:
            ORbase, ORdecors = self._separateORtypes(OR)
            if ORbase != None:
                ORtypes.append( (ORbase, ORdecors) )

        return ORtypes, ANDtypes, index

    def _separateORtypes(self, ORtype):
        """
        Separates ORtype (i.e. "#6X4R+0") into
        a base and decorators (i.e. '#6', ['X4','R','+0'] )
        """
        # special case 1: wild card
        if ORtype == '*':
            return None, []

        # Save regular expression for SMIRKS decorators
        reg = r'(!?[#]\d+|!?[aA]|!?[DHjrVX^]\d+|!?[R+-]\d*|!?[@]\d+|!?[@]@?)'

        # special case 2: replacement string has $ initially
        if ORtype[0] == '$':
            ampersand = ORtype.find('&')
            if ampersand == -1:
                # no & means no OR decorators with the replacement string
                return ORtype, []
            else: # has ampersand
                base = ORtype[:ampersand]
                ORdecorators = re.findall(reg, ORtype[ampersand:])
                return base, ORdecorators

        split = re.findall(reg, ORtype)
        if split:
            return split[0], split[1:]

        return None, []

    def _getBondInfo(self, bond):
        """
        given bond strings returns ORtypes and ANDtypes
        """
        split = bond.split(';')

        ANDtypes = split[1:]
        ORList = split[0].split(',')
        ORtypes = list()
        for OR in ORList:
            if OR[0] != '~':
                ORtypes.append( (OR[0], list(OR[1:]) ) )

        return ORtypes, ANDtypes

    def asSMIRKS(self, smarts = False):
        """
        Returns a SMIRKS representation of the chemical environment

        Parameters
        -----------
        smarts: optional, boolean
            if True, returns a SMARTS instead of SMIRKS without index labels
        """
        return self._asSMIRKS(None, None, smarts)

    def _asSMIRKS(self, initialAtom = None, neighbors = None, smarts = False):
        """Return a SMIRKS representation of the chemical environment.

        Parameters
        -----------
        initalAtom = optional, atom object
            This is randomly selected if not chosen.
        neighbors = optional, list of atom objects
            This is all of the initalAtom neighbors if not specified
            generally this is used only for the recursive calls
            so initial atoms are not reprinted
        smarts = optional, boolean
            if True, returns a SMARTS string instead of SMIRKS
        """
        # If empty chemical environment
        if len(self._graph.nodes()) == 0:
            return ""

        if initialAtom == None:
            initialAtom = self.getAtoms()[0]

        if neighbors == None:
            neighbors = self._graph.neighbors(initialAtom)

        # sort neighbors to guarantee order is constant
        neighbors = sorted(neighbors, key=lambda atom: atom.asSMIRKS())

        # initialize smirks for starting atom
        if smarts:
            smirks = initialAtom.asSMARTS()
        else:
            smirks = initialAtom.asSMIRKS()

        # loop through neighbors
        for idx, neighbor in enumerate(neighbors):
            # get the SMIRKS for the bond between these atoms
            # bonds are the same if smarts or smirks
            bondSMIRKS = self._graph.edge[initialAtom][neighbor]['bond'].asSMIRKS()

            # Get the neighbors for this neighbor
            new_neighbors = self._graph.neighbors(neighbor)
            # Remove initialAtom so it doesn't get reprinted
            new_neighbors.remove(initialAtom)

            # Call asSMIRKS again to get the details for that atom
            atomSMIRKS = self._asSMIRKS(neighbor, new_neighbors, smarts)

            # Use ( ) for branch atoms (all but last)
            if idx < len(neighbors) - 1:
                smirks += '(' + bondSMIRKS + atomSMIRKS + ')'
            # This is for the atoms that are a part of the main chain
            else:
                smirks += bondSMIRKS + atomSMIRKS

        return smirks

    def selectAtom(self, descriptor = None):
        """Select a random atom fitting the descriptor.

        Paramters
        ---------
        descriptor: optional, None
            None - returns any atom with equal probability
            int - will return an atom with that index
            'Indexed' - returns a random indexed atom
            'Unindexed' - returns a random unindexed atom
            'Alpha' - returns a random alpha atom
            'Beta' - returns a random beta atom

        Returns
        --------
        a single Atom object fitting the description
        or None if no such atom exists
        """
        if descriptor == None:
            return random.choice(self._graph.nodes())

        if type(descriptor) is int:
            for atom in self.getIndexedAtoms():
                if atom.index == descriptor:
                    return atom
            return None

        atoms = self.getComponentList(descriptor,'atom')
        if len(atoms) == 0:
            return None

        return random.choice(atoms)

    def getComponentList(component_type, descriptor = None):
        """
        Returns a list of atoms or bonds matching the descriptor

        Parameters
        -----------
        component_type: string: 'atom' or 'bond'
        descriptor: string, optional
            'all', 'Indexed', 'Unindexed', 'Alpha', 'Beta'

        Returns
        -------
        """
        if descriptor != None:
            d = descriptor.lower()

        if component_type.lower() == 'atom':
            if d == 'indexed':
                return self.getIndexedAtoms()
            elif d == 'unindexed':
                return self.getUnindexedAtoms()
            elif d == 'alpha':
                return self.getAlphaAtoms()
            elif d == 'beta':
                return self.getBetaAtoms()
            else:
                return self.getAtoms()
        if d == 'indexed':
            return self.getIndexedBonds()
        elif d == 'unindexed':
            return self.getUnindexedBonds()
        elif d == 'alpha':
            return self.getAlphaBonds()
        elif d == 'beta':
            return self.getBetaBonds()

        return self.getBonds()

    def selectBond(self, descriptor = None):
        """Select a random bond fitting the descriptor.

        Paramters
        ---------
        descriptor: optional, None
            None - returns any bond with equal probability
            int - will return an bond with that index
            'Indexed' - returns a random indexed bond
            'Unindexed' - returns a random unindexed bond
            'Alpha' - returns a random alpha bond
            'Beta' - returns a random beta bond

        Returns
        --------
        a single Bond object fitting the description
        or None if no such atom exists
        """
        if type(descriptor) is int:
            for bond in self.getIndexedBonds():
                if bond.index == descriptor:
                    return bond
            return None

        bonds = self.getComponentList(descriptor,'bond')
        if len(bonds) == 0:
            return None

        return random.choice(bonds)

    def addAtom(self, bondToAtom, bondORtypes= None, bondANDtypes = None,
            newORtypes = None, newANDtypes = None, newAtomIndex = None,
            beyondBeta = False):
        """Add an atom to the specified target atom.

        Parameters
        -----------
        bondToAtom: atom object, required
            atom the new atom will be bound to
        bondORtypes: list of tuples, optional
            strings that will be used for the ORtypes for the new bond
        bondANDtypes: list of strings, optional
            strings that will be used for the ANDtypes for the new bond
        newORtypes: list of strings, optional
            strings that will be used for the ORtypes for the new atom
        newANDtypes: list of strings, optional
            strings that will be used for the ANDtypes for the new atom
        newAtomIndex: int, optional
            integer label that could be used to index the atom in a SMIRKS string
        beyondBeta: boolean, optional
            if True, allows bonding beyond beta position

        Returns
        --------
        newAtom: atom object for the newly created atom
        """
        if bondToAtom == None:
            if len(self._graph.nodes()) > 0:
                return None
            newType = newAtomIndex
            if newType == None:
                newType = 0

            newAtom = self.Atom(newORtypes, newANDtypes, newAtomIndex)
            self._graph.add_node(newAtom, atom_type = newType)
            return newAtom

        # Check if we can get past beta position
        bondToType = self._graph.node[bondToAtom]['atom_type']
        if bondToType < 0 and not beyondBeta:
            return None

        # determine the type integer for the new atom and bond
        if newAtomIndex != None:
            newType = newAtomIndex
            bondType = max(newType, bondToType) - 1
        else:
            if bondToType > 0:
                newType = 0
            else:
                newType = bondToType - 1
            bondType = newType

        # create new bond
        newBond = self.Bond(bondORtypes, bondANDtypes)

        # create new atom
        newAtom = self.Atom(newORtypes, newANDtypes, newAtomIndex)

        # Add node for newAtom
        self._graph.add_node(newAtom, atom_type = newType)

        # Connect original atom and new atom
        self._graph.add_edge(bondToAtom, newAtom, bond = newBond, bond_type = bondType)
        newBond._bond_type = bondType

        return newAtom

    def removeAtom(self, atom, onlyEmpty = False):
        """Remove the specified atom from the chemical environment.
        if the atom is not indexed for the SMIRKS string or
        used to connect two other atoms.

        Parameters
        ----------
        atom: atom object, required
            atom to be removed if it meets the conditions.
        onlyEmpty: boolean, optional
            True only an atom with no ANDtypes and 1 ORtype can be removed

        Returns
        --------
        Boolean True: atom was removed, False: atom was not removed
        """
        # labeled atoms can't be removed
        if atom.index is not None:
            # print("Cannot remove labeled atom %s" % atom.asSMIRKS())
            return False

        # Atom connected to more than one other atom cannot be removed
        if len(self._graph.neighbors(atom)) > 1:
            #print("Cannot remove atom %s because it connects two atoms" % atom.asSMIRKS())
            return False

        # if you can remove "decorated atoms" remove it
        if not onlyEmpty:
            self._graph.remove_node(atom)
            return True

        if len(atom.ANDtypes) > 0:
            return False
        elif len(atom.ORtypes) > 1:
            return False

        self._graph.remove_node(atom)
        return True

    def getAtoms(self):
        """
        Returns
        -------
        list of atoms in the environment
        """
        return self._graph.nodes()

    def getBonds(self):
        """
        Returns
        --------
        a complete list of bonds in the fragment
        """
        return [self._graph.edge[a1][a2]['bond'] for (a1,a2) in self._graph.edges()]

    def getBond(self, atom1, atom2):
        """
        Get bond betwen two atoms

        Parameters
        -----------
        atom1 and atom2: atom objects

        Returns
        --------
        bond object between the atoms or None if no bond there
        """
        if atom2 in self._graph.edge[atom1]:
            return self._graph.edge[atom1][atom2]['bond']
        else:
            return None

    def getIndexedAtoms(self):
        """
        returns the list of Atom objects with an index
        """
        index_atoms = []
        for atom, info in self._graph.nodes_iter(True):
            if info['atom_type'] > 0:
                index_atoms.append(atom)
        return index_atoms

    def getUnindexedAtoms(self):
        """
        returns a list of Atom objects that are not indexed
        """
        unindexed_atoms = []
        for atom, info in self._graph.nodes_iter(True):
            if info['atom_type'] < 1:
                unindexed_atoms.append(atom)
        return unindexed_atoms

    def getAlphaAtoms(self):
        """
        Returns a list of atoms alpha to any indexed atom
            that are not also indexed
        """
        alpha_atoms = []
        for atom, info in self._graph.nodes_iter(True):
            if info['atom_type'] == 0:
                alpha_atoms.append(atom)

        return alpha_atoms

    def getBetaAtoms(self):
        """
        Returns a list of atoms beta to any indexed atom
            that are not alpha or indexed atoms
        """
        beta_atoms = []
        for atom, info in self._graph.nodes_iter(True):
            if info['atom_type'] == -1:
                beta_atoms.append(atom)
        return beta_atoms

    def getIndexedBonds(self):
        """
        Returns a list of Bond objects that connect two indexed atoms
        """
        indexedBonds = []
        for bond in self.getBonds():
            if bond._bond_type > 0:
                indexedBonds.append(bond)
        return indexedBonds

    def getUnindexedBonds(self):
        """
        Returns a list of Bond objects that connect
            an indexed atom to an unindexed atom
            two unindexed atoms
        """
        unindexedBonds = []
        for bond in self.getBonds():
            if bond._bond_type < 1:
                unindexedBonds.append(bond)
        return unindexedBonds

    def getAlphaBonds(self):
        """
        Returns a list of Bond objects that connect
            an indexed atom to alpha atoms
        """
        alphaBonds = []
        for bond in self.getBonds():
            if bond._bond_type == 0:
                alphaBonds.append(bond)
        return alphaBonds

    def getBetaBonds(self):
        """
        Returns a list of Bond objects that connect
            alpha atoms to beta atoms
        """
        betaBonds = []
        for bond in self.getBonds():
            if bond._bond_type == -1:
                betaBonds.append(bond)
        return betaBonds

    def isAlpha(self, component):
        """
        Takes an atom or bond are returns True if it is alpha to an indexed atom
        """
        if component._atom:
            return self._graph.node[component]['atom_type'] == 0
        else:
            return component._bond_type == 0

    def isBeta(self, component):
        """
        Takes an atom or bond are returns True if it is beta to an indexed atom
        """
        if component._atom:
            return self._graph.node[component]['atom_type'] == 1
        else:
            return component._bond_type == 1

    def getType(self):
        """
        Uses number of indexed atoms and bond connectivity
        to determine the type of chemical environment

        Returns
        -------
        chemical environemnt type:
            'VdW', 'Bond', 'Angle', 'Torsion', 'Improper'
            None if number of indexed atoms is 0 or > 4
        """
        index_atoms = self.getIndexedAtoms()
        natoms = len(index_atoms)

        if natoms == 1:
            return "VdW"
        if natoms == 2:
            return "Bond"
        if natoms == 3:
            return "Angle"
        if natoms == 4:
            atom2 = index_atoms[1]
            atom4 = index_atoms[3]
            bond24 = self.getBond(atom2, atom4)
            if bond24 is not None:
                return "Improper"
            return "Torsion"
        else:
            return None

    def getValence(atom):
        """
        Returns the valence (number of neighboring atoms)
        around the given atom
        """
        return len(self._graph.neighbors(atom))

    def getBondOrder(atom):
        """
        Returns minimum bond order around a given atom
        0 if atom has no neighbors
        aromatic bonds count as 1.5
        any bond counts as 1.0
        """
        order = 0.
        for (a1, a2, info) in self._graph.edges_iter(atom, True):
            order += info['bond'].getOrder()
        return order

class AtomChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching one labeled atom.

    """
    def __init__(self, smirks = None, label = None):
        """Initialize a chemical environment corresponding to matching a single atom.

        Parameters
        -----------
        smirks: string, optional
            if not None then the Environment is from this
            otherwise, it is an empty Atom corresponding to "[*:1]"
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything

        For example:
            # create an atom that is carbon, nitrogen, or oxygen with no formal charge
            atom = AtomChemicalEnvironment([['#6', '#7', '#8'], ['+0']])
            print atom.asSMIRKS()
            # prints: "[#6,#7,#8;+0:1]"
        """
        # Initialize base class
        if smirks == None:
            smirks = "[*:1]"

        super(AtomChemicalEnvironment,self).__init__(smirks, label)

        self.atom1 = self.selectAtom(1)
        if self.atom1 == None:
            raise Exception("AtomChemicalEnvironments need an indexed atom None found in %s" % smirks)

    def asSMIRKS(self, smarts = False):
        """
        Returns a SMIRKS representation of the chemical environment

        Parameters
        -----------
        smarts: optional, boolean
            if True, returns a SMARTS instead of SMIRKS without index labels
        """
        return self._asSMIRKS(self.atom1, None, smarts)

    def asAtomtypeSMARTS(self):
        """
        Makes SMARTS string for one atom

        Returns
        --------
        string for the SMARTS string for the first atom (labeled with index :1)

        This is a single atom with neighbors as decorators in the form:
        [atom1$(*~neighbor1)$(*~neighbor2)...]
        """
        # smarts for atom1 without last ']'
        smarts = self.atom1.asSMARTS()[:-1]

        for idx, neighbor in enumerate(self._graph.neighbors(self.atom1)):
            new_neighbors = self._graph.neighbors(neighbor)
            new_neighbors.remove(self.atom1)

            bondSMARTS = self._graph.edge[self.atom1][neighbor]['bond'].asSMARTS()
            neighborSMARTS = self._asSMIRKS(neighbor, new_neighbors, True)

            smarts += '$(*' + bondSMARTS + neighborSMARTS + ')'

        return smarts + ']'

class BondChemicalEnvironment(AtomChemicalEnvironment):
    """Chemical environment matching two labeled atoms (or a bond).
    """
    def __init__(self, smirks = None, label = None):
        """Initialize a chemical environment corresponding to matching two atoms (bond).

        Parameters
        -----------
        smirks: string, optional
            if not None then the Environment is from this
            otherwise, it is an empty Bond corresponding to "[*:1]~[*:2]"
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything

        """
        # Initialize base class
        if smirks == None:
            smirks = "[*:1]~[*:2]"

        super(BondChemicalEnvironment,self).__init__(smirks, label)

        # Add initial atom
        self.atom2 = self.selectAtom(2)
        if self.atom2 == None:
            raise Exception("Error: Bonds need 2 indexed atoms, there were not enough in %s" % smirks)

        self.bond1 = self._graph.edge[self.atom1][self.atom2]['bond']
        self.bond1.index = 1

class AngleChemicalEnvironment(BondChemicalEnvironment):
    """Chemical environment matching three marked atoms (angle).
    """
    def __init__(self, smirks = None, label = None):

        """Initialize a chemical environment corresponding to matching three atoms.

        Parameters
        -----------
        smirks: string, optional
            if not None then the Environment is from this
            otherwise, it is an empty Angle corresponding to "[*:1]~[*:2]~[*:3]"
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        """
        if smirks == None:
            smirks = "[*:1]~[*:2]~[*:3]"

        # Initialize base class
        super(AngleChemicalEnvironment,self).__init__(smirks, label)

        # Add initial atom
        self.atom3 = self.selectAtom(3)
        if self.atom3 == None:
            raise Exception("Error: Angles need 3 indexed atoms, there were not enough in %s" % smirks)
        self.bond2 = self._graph.edge[self.atom2][self.atom3]['bond']
        self.bond2.index = 2

class TorsionChemicalEnvironment(AngleChemicalEnvironment):
    """Chemical environment matching four marked atoms (torsion).
    """
    def __init__(self, smirks = None, label = None):
        """Initialize a chemical environment corresponding to matching four atoms (torsion).

        Parameters
        -----------
        smirks: string, optional
            if not None then the Environment is from this
            otherwise, it is an empty Torsion corresponding to
            "[*:1]~[*:2]~[*:3]~[*:4]"
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        """
        if smirks == None:
            smirks = "[*:1]~[*:2]~[*:3]~[*:4]"
        # Initialize base class
        super(TorsionChemicalEnvironment,self).__init__(smirks, label)

        # Add initial atom
        self.atom4 = self.selectAtom(4)
        if self.atom4 == None:
            raise Exception("Error: Torsion need 4 indexed atoms, there were not enough in %s" % smirks)
        self.bond3 = self._graph.edge[self.atom3][self.atom4]['bond']
        self.bond3.index = 3

class ImproperChemicalEnvironment(AngleChemicalEnvironment):
    """Chemical environment matching four marked atoms (improper).
    """
    def __init__(self, smirks = None, label = None):
        """Initialize a chemical environment corresponding four atoms (improper).

        Parameters
        -----------
        smirks: string, optional
            if not None then the Environment is from this
            otherwise, it is an empty Improper corresponding to
            "[*:1]~[*:2](~[*:3])~[*:4]"
        label = anything, optional
            intended to be used to label this chemical environment
            could be a string, int, or float, or anything
        """
        if smirks == None:
            smirks = "[*:1]~[*:2](~[*:3])~[*:4]"
        # Initialize base class
        super(ImproperChemicalEnvironment,self).__init__(smirks, label)

        # Add initial atom
        self.atom4 = self.selectAtom(4)
        if self.atom4 == None:
            raise Exception("Error: Improper need 4 indexed atoms, there were not enough in %s" % smirks)
        self.bond3 = self._graph.edge[self.atom2][self.atom4]['bond']
        self.bond3.index = 3
