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

            Parameters
            -----------
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
        Properties
        -----------
        bases : set of string
            The base types that will be combined with logical OR
        decorators : set of string
            The decorators that will be combined with logical AND

        """
        # implementation similar to Atom but for bonds connecting atoms

        def __init__(self, bases = None, decorators = None):
            """
            Parameters
            -----------
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
            Returns 
            --------
            the same as asSMARTS()
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

    def selectAtom(self, atomIndex = None):
        """Select an atom with uniform probability.

        Paramters
        ---------
        atomIndex: int, optional
            if None returns a random atom, otherwise returns atom at that index

        Returns 
        --------
        a random atom(node)
        """
        if atomIndex == None:
            return random.choice(self._graph.nodes())
        atoms = self.getAtoms()
        indices = [a.index for a in atoms]
        if atomIndex in indices:
            return atoms.index(atomIndex)
        # TODO: index is not specified raise error instead?
        else:
            return None

    def selectBond(self, atomIndex1 = None, atomIndex2 = None):
        """Select a bond with uniform probability.
        Bonds are found by the two atoms that define them

        Parameters
        ----------
        atomIndex1 and atomIndex2: int, optional
            indices for atoms framing the bond

        Returns
        --------
        atom1 and atom2
            Atom objects that are on either end of the bond
        bond
            Bond object connencting atoms 
        """
        # TODO: Handle exceptions associated with the input atom indices

        # get atom1, if atomIndex1 is None then this is random
        atom1 = self.selectAtom(atomIndex1)

        # if antomIndex2 is None
        if atomIndex2 == None:
            atom2 = random.choice(self._graph.neighbors(atom1))
        else:
            atom2 = self.selectAtom(atomIndex2)
            if not atom2 in self._graph.neighbors(atom1):
                raise Exception("Error Atom Mismatch:
                        Atoms1 (%s) is not bonded to Atom2 (%s)" % (atom1.asSMIRKS(), atom2.asSMIRKS()))

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

        Returns
        --------
        Boolean True: atom was removed, False: atom was not removed
        """
        if atom.index is not None:
            print("Cannot remove labeled atom %s" % atom.asSMIRKS())
            return False

        elif len(self._graph.neighbors(atom)) > 1:
            print("Cannot remove atom %s because it connects two atoms" % atom.asSMIRKS())
            return False

        else:
            # Remove atom (removes associated bonds)
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
    def __init__(self, AtomInfo = [None, None]):
        """Initialize a chemical environment corresponding to matching a single atom.

        Parameters
        -----------
        AtomInfo: list of sets, optional
            Comes in the form [AtomBases, AtomDecors]
            AtomBases: descriptors for the first atom that are connected with logical operation OR
            AtomDecors: descriptors for the first atom that are connected with the logical operation AND

        For example:
            # create an atom that is carbon, nitrogen, or oxygen with no formal charge
            atom = AtomChemicalEnvironment([['#6', '#7', '#8'], ['+0']])
            print atom.asSMIRKS()
            # prints: "[#6,#7,#8;+0:1]"
        """
        # Initialize base class
        super(AtomChemicalEnvironment,self).__init__()
        self.atom1 = self.Atom(1, AtomInfo[0], AtomInfo[1])
        self._graph.add_node(self.atom1)

    def asSMARTS(self):
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
            neighborSMARTS = self.asSMIRKS(neighbor, new_neighbors, True)

            smarts += '$(*' + bondSMARTS + neighborSMARTS + ')'

        return smarts + ']'

class BondChemicalEnvironment(AtomChemicalEnvironment):
    """Chemical environment matching two labeled atoms (or a bond).
    """
    def __init__(self, Atom1Info = [None, None],
            BondInfo = [None, None],
            Atom2Info = [None, None])
        """Initialize a chemical environment corresponding to matching two atoms (bond).

        Parameters
        -----------
        Atom1Info, Atom2Info: list of sets, optional
            Comes in the form [AtomBases, AtomDecors]
            AtomBases: descriptors for the first atom that are connected with logical operation OR
            AtomDecors: descriptors for the first atom that are connected with the logical operation AND
        BondInfo: list of sets, optional
            In the form [BondBases, BondDecors] similar to atom information

        For example:
            # create a tetravalent carbon connected with a single bond to oxygen
            Atom1Info = [['#6'], ['X4']]  
            BondInfo = [['-'], None]  
            Atom2Info = [['#8'], None]  

            bond = BondChemicalEnvironment(Atom1Info, BondInfo, Atom2Info)
            print bond.asSMIRKS()
            # prints: "[#6;X4:1]-[#8:2]" 
        """
        # Initialize base class
        super(BondChemicalEnvironment,self).__init__(Atom1Info)

        # Add initial atom
        self.atom2 = self.addAtom(self.atom1, BondInfo[0], BondInfo[1], AtomInfo[0], AtomInfo[1], 2)

class AngleChemicalEnvironment(BondChemicalEnvironment):
    """Chemical environment matching three marked atoms (angle).
    """
    def __init__(self, Atom1Info = [None, None], Bond1Info = [None, None], 
            Atom2Info = [None, None], Bond2Info = [None, None], Atom3Info = [None, None]):

        """Initialize a chemical environment corresponding to matching three atoms.

        Parameters
        -----------
        Atom1Info, Atom2Info, Atom3Info: list of sets, optional
            Comes in the form [AtomBases, AtomDecors]
            AtomBases: descriptors for the first atom that are connected with logical operation OR
            AtomDecors: descriptors for the first atom that are connected with the logical operation AND
        Bond1Info and Bond2Info: list of sets, optional
            In the form [BondBases, BondDecors] similar to atom information

        For example:
            # create an angle where the center atom is a neutral trivalent carbon 
            Atom2Info = [['#6X3'], ['+0']]
            angle = AngleChemicalEnvironment(Atom2Info = Atom2Info)
            print angle.asSMIRKS()
            # "[*:1]~[#6X3;+0:2]~[*:3]"
        """
        # Initialize base class
        super(AngleChemicalEnvironment,self).__init__(Atom1Info, Bond1Info, Atom2Info)

        # Add initial atom
        self.atom3 = self.addAtom(self.atom2, Bond2Info[0], Bond2Info[1], Atom3Info[0], Atom3Info[1], 3)

class TorsionChemicalEnvironment(AngleChemicalEnvironment):
    """Chemical environment matching four marked atoms (torsion).
    """
    def __init__(self, Atom1Info = [None, None], Bond1Info = [None, None], 
            Atom2Info = [None, None], Bond2Info = [None, None], 
            Atom3Info = [None, None], Bond3Info = [None, None], Atom4Info = [None, None]):
        """Initialize a chemical environment corresponding to matching four atoms (torsion).
        
        Parameters
        -----------
        Atom1Info, Atom2Info, Atom3Info, Atom4Info: list of sets, optional
            Comes in the form [AtomBases, AtomDecors]
            AtomBases: descriptors for the first atom that are connected with logical operation OR
            AtomDecors: descriptors for the first atom that are connected with the logical operation AND
        Bond1Info and Bond2Info, Bond3Info: list of sets, optional
            In the form [BondBases, BondDecors] similar to atom information

        For example:
            # Create a torsion centered around two tetravalent carbons with single ring bonds
            CarbonInfo = [['#6'], ['X4']
            BondInfo = [['-'], ['@']]
            torsion = TorsionChemicalEnvironment(Atom2Info = CarbonInfo, Bond2Info = BondInfo, Atom3Info = CarbonInfo)
            print torsion.asSMIRKS()
            # "[*:1]~[#6X4:2]-;@[#6X4:3]~[*:4]"
        """
        # Initialize base class
        super(TorsionChemicalEnvironment,self).__init__(Atom1Info, Bond1Info, 
                Atom2Info, Bond2Info, Atom3Info)

        # Add initial atom
        self.atom4 = self.addAtom(self.atom3, Bond3Info[0], Bond3Info[1], Atom4Info[0], Atom4Info[1], 4)

class ImproperChemicalEnvironment(AngleChemicalEnvironment):
    """Chemical environment matching four marked atoms (improper).
    """
    def __init__(self, Atom1Info = [None, None], Bond1Info = [None, None], 
            Atom2Info = [None, None], Bond2Info = [None, None], 
            Atom3Info = [None, None], Bond3Info = [None, None], Atom4Info = [None, None]):
        """Initialize a chemical environment corresponding to matching four atoms (improper).
        
        Parameters
        -----------
        Atom1Info, Atom2Info, Atom3Info, Atom4Info: list of sets, optional
            Comes in the form [AtomBases, AtomDecors]
            AtomBases: descriptors for the first atom that are connected with logical operation OR
            AtomDecors: descriptors for the first atom that are connected with the logical operation AND
        Bond1Info and Bond2Info, Bond3Info: list of sets, optional
            In the form [BondBases, BondDecors] similar to atom information

        For example:

        """
        # TODO: add improper example after talking to Christopher about numbering
        # Initialize base class
        super(TorsionChemicalEnvironment,self).__init__(Atom1Info, Bond1Info, 
                Atom2Info, Bond2Info, Atom3Info)

        # Add initial atom
        self.atom4 = self.addAtom(self.atom2, Bond3Info[0], Bond3Info[1], Atom4Info[0], Atom4Info[1], 4)

