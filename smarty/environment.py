#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

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
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import networkx as nx
import random

class ChemicalEnvironment(object):
    """Chemical environment abstract base class that matches an atom, bond, angle, etc.
    """
    class Atom(object):
        """Atom representation, which may have some base (logical OR) and decorator (logical AND) properties.

        Type will be OR(bases) and AND(decorators).

        Properties
        -----------
        bases : set of string
        The base types that will be combined with logical OR
        decorators : set of string
        The decorators that will be combined with logical AND
        """
        def __init__(self, index = None, bases = None, decorators = None):
            """Initialize an Atom object with empty bases decorators.

            ARGUMENTS:
            index : int, optional, default=None
            If not None, the specified index will be attached as a SMIRKS index (e.g. '[#6:1]')
            bases: set of strings, optional, default = None
            strings that will be OR'd together in a SMARTS
            decorators: set of str, optional, default = None
            strings that will be AND'd together in a SMARTS
            """
            # TODO: Change names for bases and decorators. 

            # Set of strings that will be OR'd together
            if bases is  None:
                self.bases = set()
            else:
                self.bases = set(bases)

            # Set of strings that will be AND'd to the the end 
            # Probably decorators
            if decorators is None:
                self.decorators = set()
            else:
                self.decorators = set(decorators)

            self.index = index

        def asSMARTS(self):
            """Return the atom representation as SMARTS.

            Returns
            --------
            smarts : str
            The SMARTS string for just this atom
            """

            smarts = '['

            # Add the OR'd features
            if self.bases:
                smarts += ','.join(self.bases)
            else:
                smarts += '*'

            if self.decorators:
                smarts += ';' + ';'.join(self.decorators)

            return smarts + ']'

        def asSMIRKS(self):
            """Return the atom representation as SMIRKS.

            Parameters
            ----------

            Returns
            --------
            smirks : str
            The SMIRKS string
            """
            smirks = self.asSMARTS()

            # No index specified so SMIRKS = SMARTS 
            if self.index is None:
                return smirks 

            # Add label to the end of SMARTS
            else: 
                return smirks[:-1] + ':' + str(self.index) + smirks[-1]


    class Bond(object):
        """Bond representation, which may have base (OR) and decorator (AND) types.
        """
        # implementation similar to Atom but for bonds connecting atoms

        def __init__(self, bases = None, decorators = None):
            """
            ARGUMENTS:
            bases: set of strings, optional, default = None
            strings that will be OR'd together in a SMARTS
            decorators: set of str, optional, default = None
            strings that will be AND'd together in a SMARTS 
            """

            # Make set of bases
            if bases is None:
                self.bases = set()
            else:
                self.bases = set(bases)

            # Make set of decorators
            if decorators is None:
                self.decorators = set()
            else:
                self.decorators = set(decorators)

        def asSMARTS(self):
            """Return the atom representation as SMARTS.

            Returns
            --------
            smarts : str
            The SMARTS string for just this atom
            """
            if self.bases:
                smarts = ','.join(self.bases)
            else:
                smarts = '~'

            if self.decorators:
                smarts += ';' + ';'.join(self.decorators)

            return smarts

        def asSMIRKS(self):
            """
            returns the same as asSMARTS
            for consistency asSMARTS() or asSMIRKS() can be called
            """
            return self.asSMARTS()

    def __init__(self):
        """Initialize a chemical environment abstract base class.

        ARGUMENTS

        """
        # Create an empty graph which will store Atom objects.
        self._graph = nx.Graph()

        # Storing sets of information
        # self.atomBaseSet = set(atomBaseSet)
        # self.atomDecorSet = set(atomDecorSet)
        # self.bondBaseSet = set(bondBaseSet)
        # self.bondDecorSet = set(bondDecorSet)

    def asSMIRKS(self, initialAtom = None, neighbors = None):
        """Return a SMIRKS representation of the chemical environment.
        """
        if initialAtom is None:
            initialAtom = self.getAtoms()[0]

        if neighbors is None:
            neighbors = self._graph.neighbors(initialAtom)

        smirks = initialAtom.asSMIRKS()

        for idx, neighbor in enumerate(neighbors):
            bondSMIRKS = self._graph.edge[initialAtom][neighbor]['bond'].asSMIRKS()

            new_neighbors = self._graph.neighbors(neighbor)
            new_neighbors.remove(initialAtom)

            atomSMIRKS = self.asSMIRKS(neighbor, new_neighbors)

            if idx < len(neighbors) - 1:
                smirks += '(' + bondSMIRKS + atomSMIRKS + ')'
            else:
                smirks += bondSMIRKS + atomSMIRKS

        return smirks

    def selectAtom(self):
        """Select an atom with uniform probability.

        Returns atom(node)
        """
        return random.choice(self._graph.nodes())

    def selectBond(self):
        """Select a bond with uniform probability.
        Bonds are found by the two atoms that define them

        RETURNS
        (atom1, atom2) 
        """
        return random.choice(self._graph.edges())

    def addAtom(self, bondToAtom = None, bondBases = None, bondDecorators = None, 
            newAtomIndex = None, newBases = None, newDecorators = None):
        """Add an atom to the specified target atom.

        Returns the new atom
        """
        if bondToAtom is None:
            bondToAtom = self.selectAtom()

        # create new bond
        newBond = self.Bond(bondBases, bondDecorators)

        # create new atom
        newAtom = self.Atom(newAtomIndex, newBases, newDecorators)

        # Add node for newAtom
        self._graph.add_node(newAtom)

        # Connect original atom and new atom
        self._graph.add_edge(bondToAtom, newAtom, bond = newBond)

        return newAtom

    def removeAtom(self, atom):
        """Remove the specified atom from the chemical environment, along with its associated bond.
        """
        pass

    def getAtoms(self):
        return self._graph.nodes()

    def getBonds(self):
        return self._graph.edges()

class AtomChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching one labeled atom.
    """
    def __init__(self, initAtomBases = [None], initAtomDecors = [None]):
        """Initialize a chemical environment corresponding to matching a single atom.
        """
        # Initialize base class
        super(AtomChemicalEnvironment,self).__init__()
        atom1 = self.Atom(1, initAtomBases[0], initAtomDecors[0])
        self._graph.add_node(atom1)

class BondChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching two labeled atoms (or a bond).
    """
    def __init__(self, 
            initAtomBases = [None, None], initAtomDecors = [None, None], 
            initBondBases = [None], initBondDecors = [None]):
        """Initialize a chemical environment corresponding to matching two atoms (bond).
        """
        # Initialize base class
        super(BondChemicalEnvironment,self).__init__()

        # Add initial atom
        atom1 = self.Atom(1, initAtomBases[0], initAtomDecors[0])
        self._graph.add_node(atom1)

        atom2 = self.addAtom(atom1, initBondBases[0], initBondDecors[0], 2, initAtomBases[1], initAtomDecors[1])

class AngleChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching three marked atoms (angle).
    """
    def __init__(self, initAtomBases = [None, None, None], initAtomDecors = [None, None, None], 
            initBondBases = [None, None], initBondDecors = [None, None]):
        """Initialize a chemical environment corresponding to matching three atoms.
        """
        # Initialize base class
        super(AngleChemicalEnvironment,self).__init__()

        # Add initial atom
        atom1 = self.Atom(1, initAtomBases[0], initAtomDecors[0])
        self._graph.add_node(atom1)

        atom2 = self.addAtom(atom1, initBondBases[0], initBondDecors[0], 2, initAtomBases[1], initAtomDecors[1])
        atom3 = self.addAtom(atom2, initBondBases[1], initBondDecors[1], 3, initAtomBases[2], initAtomDecors[2])

class TorsionChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching four marked atoms (torsion).
    """
    def __init__(self, 
            initAtomBases = [None, None, None, None], initAtomDecors = [None, None, None, None], 
            initBondBases = [None, None, None], initBondDecors = [None, None, None]):
        """Initialize a chemical environment corresponding to matching four atoms.
        """
        # Initialize base class
        super(TorsionChemicalEnvironment,self).__init__()

        # Add initial atom
        atom1 = self.Atom(1, initAtomBases[0], initAtomDecors[0])
        self._graph.add_node(atom1)

        atom2 = self.addAtom(atom1, initBondBases[0], initBondDecors[0], 2, initAtomBases[1], initAtomDecors[1])
        atom3 = self.addAtom(atom2, initBondBases[1], initBondDecors[1], 3, initAtomBases[2], initAtomDecors[2])
        atom4 = self.addAtom(atom3, initBondBases[2], initBondDecors[2], 4, initAtomBases[3], initAtomDecors[3])

class ImproperChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching four marked atoms (improper).
    """
    def __init__(self, 
        initAtomBases = [None, None, None, None], initAtomDecors = [None, None, None, None], 
        initBondBases = [None, None, None], initBondDecors = [None, None, None]):
        """Initialize a chemical environment corresponding to matching four atoms
        connected as an improper torsion.
        """
        # Initialize base class
        super(ImproperChemicalEnvironment,self).__init__()

        # Add initial atom
        atom1 = self.Atom(1, initAtomBases[0], initAtomDecors[0])
        self._graph.add_node(atom1)

        atom2 = self.addAtom(atom1, initBondBases[0], initBondDecors[0], 2, initAtomBases[1], initAtomDecors[1])
        atom3 = self.addAtom(atom2, initBondBases[1], initBondDecors[1], 3, initAtomBases[2], initAtomDecors[2])
        atom4 = self.addAtom(atom2, initBondBases[2], initBondDecors[2], 4, initAtomBases[3], initAtomDecors[3])

