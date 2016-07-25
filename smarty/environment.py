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
            # Set of strings that will be OR'd together
            if bases ==  None:
                self.bases = set()
            else:
                self.bases = set(bases)

            # Set of strings that will be AND'd to the the end 
            # Probably decorators
            if decorators == None:
                self.decorators = set()
            else:
                self.decorators = set(decorators)

            self.index = index

        def asSMARTS(self):
            """Return the atom representation as SMARTS.

            Returns
            --------
            smarts : str
            The SMARTS string for this atom
            """

            smarts = '['

            # Add the OR'd features
            if len(self.bases) > 0:
                smarts += ','.join(self.bases)
            else:
                smarts += '*'

            if len(self.decorators) > 0:
                smarts += ';' + ';'.join(self.decorators)

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

        def addBase(self, base):
            """
            Adds base to the set for this atom.

            Parameters
            --------
            base: string
                added to the set of bases that are OR'd together for this atom
            """
            self.bases.add(base)

        def addDecorator(self, decorator):
            """
            Adds decorator to the set for this atom.

            Parameters
            --------
            decorator: string
                added to the set of decorators that are ANd'd together for this atom
            """
            self.decorators.add(decorator)

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
            if bases == None:
                self.bases = set()
            else:
                self.bases = set(bases)

            # Make set of decorators
            if decorators == None:
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
            if len(self.bases) > 0:
                smarts = ','.join(self.bases)
            else:
                smarts = '~'

            if len(self.decorators) > 0:
                smarts += ';' + ';'.join(self.decorators)

            return smarts

        def asSMIRKS(self):
            """
            returns the same as asSMARTS()
                for consistency asSMARTS() or asSMIRKS() can be called
                for all environment objects 
            """
            return self.asSMARTS()

        def addBase(self, base):
            """
            Adds base to the set for this bond.

            Parameters
            --------
            base: string
                added to the set of bases that are OR'd together for this bond
            """
            self.bases.add(base)

        def addDecorator(self, decorator):
            """
            Adds decorator to the set for this bond.

            Parameters
            --------
            decorator: string
                added to the set of decorator that are OR'd together for this bond
            """
            self.decorators.add(decorator)

    def __init__(self):
        """Initialize a chemical environment abstract base class.

        This is an empty chemical environment. 
        Atom, Bond, Angle, Torsion, or Improper Chemical Environment 
        should be used for a filled chemical environment
        """
        # Create an empty graph which will store Atom objects.
        self._graph = nx.Graph()

    def asSMIRKS(self, initialAtom = None, neighbors = None, smarts = False):
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
        if initialAtom == None:
            initialAtom = self.getAtoms()[0]

        if neighbors == None:
            neighbors = self._graph.neighbors(initialAtom)

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
            atomSMIRKS = self.asSMIRKS(neighbor, new_neighbors, smarts)

            # Use ( ) for branch atoms (all but last)
            if idx < len(neighbors) - 1:
                smirks += '(' + bondSMIRKS + atomSMIRKS + ')'
            # This is for the atoms that are a part of the main chain
            else:
                smirks += bondSMIRKS + atomSMIRKS

        return smirks

    def selectAtom(self):
        """Select an atom with uniform probability.

        Returns 
        --------
        a random atom(node)
        """
        return random.choice(self._graph.nodes())

    def selectBond(self):
        """Select a bond with uniform probability.
        Bonds are found by the two atoms that define them

        Returns
        --------
        atom1 and atom2
            Atom objects that are on either end of the bond
        bond
            Bond object connencting atoms 
        """
        # Get a random edge
        (atom1, atom2) = random.choice(self._graph.edges())  
        # Get the bond object for that edge
        bond = self._graph.edge[atom1][atom2]['bond']
        return atom1, atom2, bond 

    def addAtom(self, bondToAtom, bondBases = None, bondDecorators = None, 
            newBases = None, newDecorators = None, newAtomIndex = None):
        """Add an atom to the specified target atom.

        Parameters
        -----------
        bondToAtom: atom object, required
            atom the new atom will be bound to
        bondBases: set of strings, optional
            strings that will be used for the bases (OR set) for the new bond
        bondDecorators: set of strings, optional
            strings that will be used for the decorators (AND set) for the new bond
        newBases: set of strings, optional
            strings that will be used for the bases (OR set) for the new atom
        newDecorators: set of strings, optional
            strings that will be used for the decorates (AND set) fro the new atom
        newAtomIndex: int, optional
            integer label that could be used to index the atom in a SMIRKS string

        Returns 
        --------
        newAtom: atom object for the newly created atom        
        """
        if bondToAtom == None:
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
        """Remove the specified atom from the chemical environment.
        if the atom is not indexed for the SMIRKS string or 
        used to connect two other atoms. 

        Parameters
        ----------
        atom: atom object, required
            atom to be removed if it meets the conditions. 
        """
        if atom.index is not None:
            print("Cannot remove labeled atom %s" % atom.asSMIRKS())

        elif len(self._graph.neighbors(atom)) > 1:
            print("Cannot remove atom %s because it connects two atoms" % atom.asSMIRKS())

        else:
            bondedTo = self._graph.neighbors(atom)[0]
            # Remove atom (removes associated bonds)
            self._graph.remove_node(atom)

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
        list of tuples for each bond in the form
        (atom1, atom2, bond object)
        """
        bondList = []
        for (atom1, atom2) in self._graph.edges():
            bond = self._graph.edge[atom1][atom2]['bond']
            bondList.append((atom1, atom2, bond))

        return bondList

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
        if self._graph.edge[atom1].hasKey(atom2):
            return self._graph.edge[atom1][atom2]['bond']
        else:
            return None
class AtomChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching one labeled atom.
    """
    def __init__(self, initAtomBases = None, initAtomDecors = None):
        """Initialize a chemical environment corresponding to matching a single atom.
        """
        # Initialize base class
        super(AtomChemicalEnvironment,self).__init__()
        self.atom1 = self.Atom(1, initAtomBases, initAtomDecors)
        self._graph.add_node(self.atom1)

    def asSMARTS(self):
        """
        Makes SMARTS string for one atom

        Returns
        --------
        the SMARTS string for the atom for atom ':1'
        """
        # smarts for atom1 without last ']'
        smarts = self.atom1.asSMARTS()[:-1]
        
        for idx, neighbor in enumerate(self._graph.neighbors(self.atom1)):
            new_neighbors = self._graph.neighbors(neighbor)
            new_neighbors.remove(self.atom1)

            bondSMARTS = self._graph.edge[self.atom1][neighbor]['bond'].asSMARTS()
            neighborSMARTS = self.asSMIRKS(neighbor, new_neighbors, True)

            if idx == 0:
                smarts += '$(*' + bondSMARTS + neighborSMARTS + ')'
            else:
                smarts += '(*' + bondSMARTS + neighborSMARTS + ')'

        return smarts + ']'

class BondChemicalEnvironment(AtomChemicalEnvironment):
    """Chemical environment matching two labeled atoms (or a bond).
    """
    def __init__(self, 
            initAtomBases = [None, None], initAtomDecors = [None, None], 
            initBondBases = None, initBondDecors = None):
        """Initialize a chemical environment corresponding to matching two atoms (bond).
        """
        # Initialize base class
        super(BondChemicalEnvironment,self).__init__(initAtomBases[0], initAtomDecors[0])

        # Add initial atom
        self.atom2 = self.addAtom(self.atom1, initBondBases, initBondDecors, initAtomBases[1], initAtomDecors[1], 2)

class AngleChemicalEnvironment(BondChemicalEnvironment):
    """Chemical environment matching three marked atoms (angle).
    """
    def __init__(self, initAtomBases = [None, None, None], initAtomDecors = [None, None, None], 
            initBondBases = [None, None], initBondDecors = [None, None]):
        """Initialize a chemical environment corresponding to matching three atoms.
        """
        # Initialize base class
        super(AngleChemicalEnvironment,self).__init__(initAtomBases[:2], initAtomDecors[:2],
                initBondBases[0], initBondDecors[0])

        # Add initial atom
        self.atom3 = self.addAtom(self.atom2, initBondBases[1], initBondDecors[1], initAtomBases[2], initAtomDecors[2], 3)

class TorsionChemicalEnvironment(AngleChemicalEnvironment):
    """Chemical environment matching four marked atoms (torsion).
    """
    def __init__(self, 
            initAtomBases = [None, None, None, None], initAtomDecors = [None, None, None, None], 
            initBondBases = [None, None, None], initBondDecors = [None, None, None]):
        """Initialize a chemical environment corresponding to matching four atoms.
        """
        # Initialize base class
        super(TorsionChemicalEnvironment,self).__init__(
                initAtomBases[:3], initAtomDecors[:3], 
                initBondBases[:2], initBondDecors[:2])

        # Add initial atom
        self.atom4 = self.addAtom(self.atom3, initBondBases[2], initBondDecors[2], initAtomBases[3], initAtomDecors[3], 4)

class ImproperChemicalEnvironment(AngleChemicalEnvironment):
    """Chemical environment matching four marked atoms (improper).
    """
    def __init__(self, 
        initAtomBases = [None, None, None, None], initAtomDecors = [None, None, None, None], 
        initBondBases = [None, None, None], initBondDecors = [None, None, None]):
        """Initialize a chemical environment corresponding to matching four atoms
        connected as an improper torsion.
        """
        # Initialize base class
        super(ImproperChemicalEnvironment,self).__init__(
                initAtomBases[:3], initAtomDecors[:3],
                initBondBases[:2], initBondDecors[:2])

        # Add initial atom
        self.atom4 = self.addAtom(self.atom2, initBondBases[2], initBondDecors[2], initAtomBases[3], initAtomDecors[3], 4)

