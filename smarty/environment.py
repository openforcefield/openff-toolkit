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
import copy

import openeye.oechem
from openeye.oechem import *

class ChemicalEnvironment(object):
    """Chemical environment abstract base class that matches an atom, bond, angle, etc.
    """
    class Atom(object):
        """Atom representation, which may have some ORtypes and ANDtypes properties.

        Properties
        -----------
        ORtypes : list of string
            The descriptor types that will be combined with logical OR
        ANDtypes : list of string
            The descriptor types  that will be combined with logical AND
        """
        def __init__(self, index = None, ORtypes = None, ANDtypes = None):
            """Initialize an Atom object with optional descriptors.

            Parameters
            -----------
            index : int, optional, default=None
                If not None, the specified index will be attached as a SMIRKS index (e.g. '[#6:1]')
            ORtypes: list of strings, optional, default = None
                strings that will be OR'd together in a SMARTS
            ANDtypes: list of str, optional, default = None
                strings that will be AND'd together in a SMARTS
            """
            # Set of strings that will be OR'd together
            if ORtypes ==  None:
                self.ORtypes = list()
            else:
                self.ORtypes = list(ORtypes)

            # Set of strings that will be AND'd to the the end
            if ANDtypes == None:
                self.ANDtypes = list()
            else:
                self.ANDtypes = list(ANDtypes)

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
            if len(self.ORtypes) > 0:
                smarts += ','.join(self.ORtypes)
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

        def addORtype(self, ORtype):
            """
            Adds ORtype to the set for this atom.

            Parameters
            --------
            ORtype: string
                added to the list of ORtypes for this atom
            """
            self.ORtypes.append(ORtype)

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
            returns a copy of the list of ORtypes for this atom
            """
            return list(copy.deepcopy(self.ORtypes))

        def setORtypes(self, newORtypes):
            """
            sets new ORtypes for this atom

            Parameters
            ----------
            newORtypes: list of strings
                strings that will be OR'd together in a SMARTS
            """
            if newORtypes == None:
                self.ORtypes = list()
            else:
                self.ORtypes = list(newORtypes)

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

    class Bond(object):
        """Bond representation, which may have ORtype and ANDtype descriptors.
        Properties
        -----------
        ORtypes : list of string
            The ORtype types that will be combined with logical OR
        ANDtypes : list of string
            The ANDtypes that will be combined with logical AND

        """
        # implementation similar to Atom but for bonds connecting atoms

        def __init__(self, ORtypes = None, ANDtypes = None):
            """
            Parameters
            -----------
            ORtypes: list of strings, optional, default = None
                strings that will be OR'd together in a SMARTS
            ANDtypes: list of str, optional, default = None
                strings that will be AND'd together in a SMARTS
            """

            # Make set of ORtypes
            if ORtypes == None:
                self.ORtypes = list()
            else:
                self.ORtypes = list(ORtypes)

            # Make set of ANDtypes
            if ANDtypes == None:
                self.ANDtypes = list()
            else:
                self.ANDtypes = list(ANDtypes)

        def asSMARTS(self):
            """Return the atom representation as SMARTS.

            Returns
            --------
            smarts : str
                The SMARTS string for just this atom
            """
            if len(self.ORtypes) > 0:
                smarts = ','.join(self.ORtypes)
            else:
                smarts = '~'

            if len(self.ANDtypes) > 0:
                smarts += ';' + ';'.join(self.ANDtypes)

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

        def addORtype(self, ORtype):
            """
            Adds ORtype to the set for this bond.

            Parameters
            --------
            ORtype: string
                added to the list of ORtypes for this bond
            """
            self.ORtypes.append(ORtype)

        def addANDtype(self, ANDtype):
            """
            Adds ANDtype to the set for this bond.

            Parameters
            --------
            ANDtype: string
                added to the list of ANDtype for this bond
            """
            self.ANDtypes.append(ANDtype)

        def getORtypes(self):
            """
            returns a copy of the list of ORtypes for this bond
            """
            return list(copy.deepcopy(self.ORtypes))

        def setORtypes(self, newORtypes):
            """
            sets new ORtypes for this bond

            Parameters
            ----------
            newORtypes: list of strings
                strings that will be OR'd together in a SMARTS
            """
            if newORtypes == None:
                self.ORtypes = list()
            else:
                self.ORtypes = list(newORtypes)

        def getANDtypes(self):
            """
            returns a copy of the list of ANDtypes for this bond
            """
            return list(copy.deepcopy(self.ANDtypes))

        def setANDtypes(self, newANDtypes):
            """
            sets new ANDtypes for this bond

            Parameters
            ----------
            newANDtypes: list of strings
                strings that will be AND'd together in a SMARTS
            """
            if newANDtypes == None:
                self.ANDtypes = list()
            else:
                self.ANDtypes = list(newANDtypes)

    def __init__(self, smirks = None):
        """Initialize a chemical environment abstract base class.
        
        smirks = string, optional
            if smirks is not None, a chemical environment is built 
            from the provided SMIRKS string
        """
        # Create an empty graph which will store Atom objects.
        self._graph = nx.Graph()
        self.label = None

        if smirks is not None:
            # check SMIRKS is parseable
            mol = OEQMol()
            if not OEParseSmarts(mol, smirks):
                raise Exception("Provides SMIRKS: %s was not parseable" % smirks)

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

        ANDtypes = split[1:]
        ORtypes = split[0].split(',')

        if len(ANDtypes) == 0:
            ANDtypes = None

        if len(ORtypes) == 1 and ORtypes[0] == '*':
            ORtypes = None
        
        return ORtypes, ANDtypes, index

    def _getBondInfo(self, bond):
        """
        given bond strings returns ORtypes and ANDtypes
        """
        split = bond.split(';')

        ANDtypes = split[1:]
        ORtypes = split[0].split(',')

        if len(ANDtypes) == 0:
            ANDtypes = None

        if len(ORtypes) == 1 and ORtypes[0] == '~':
            ORtypes = None
        
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
            return atoms[indices.index(atomIndex)]
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
                raise Exception("Error Atom Mismatch: Atoms1 (%s) is not bonded to Atom2 (%s)" % (atom1.asSMIRKS(), atom2.asSMIRKS()))

        # Get the bond object for that edge
        bond = self._graph.edge[atom1][atom2]['bond']
        return atom1, atom2, bond

    def addAtom(self, bondToAtom, bondORtypes= None, bondANDtypes = None,
            newORtypes = None, newANDtypes = None, newAtomIndex = None):
        """Add an atom to the specified target atom.

        Parameters
        -----------
        bondToAtom: atom object, required
            atom the new atom will be bound to
        bondORtypes: list of strings, optional
            strings that will be used for the ORtypes for the new bond
        bondANDtypes: list of strings, optional
            strings that will be used for the ANDtypes for the new bond
        newORtypes: list of strings, optional
            strings that will be used for the ORtypes for the new atom
        newANDtypes: list of strings, optional
            strings that will be used for the ANDtypes for the new atom
        newAtomIndex: int, optional
            integer label that could be used to index the atom in a SMIRKS string

        Returns
        --------
        newAtom: atom object for the newly created atom
        """
        # TODO: determine other requirements to check before adding atom
        if bondToAtom == None:
            if len(self._graph.nodes()) > 0:
                return None
            newAtom = self.Atom(newAtomIndex, newORtypes, newANDtypes)
            self._graph.add_node(newAtom)
            return newAtom

        # create new bond
        newBond = self.Bond(bondORtypes, bondANDtypes)

        # create new atom
        newAtom = self.Atom(newAtomIndex, newORtypes, newANDtypes)

        # Add node for newAtom
        self._graph.add_node(newAtom)

        # Connect original atom and new atom
        self._graph.add_edge(bondToAtom, newAtom, bond = newBond)

        return newAtom

    def removeAtom(self, atom, onlyEmpty = True):
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
            print("Cannot remove labeled atom %s" % atom.asSMIRKS())
            return False

        # Atom connected to more than one other atom cannot be removed
        if len(self._graph.neighbors(atom)) > 1:
            print("Cannot remove atom %s because it connects two atoms" % atom.asSMIRKS())
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
        if atom2 in self._graph.edge[atom1]:
            return self._graph.edge[atom1][atom2]['bond']
        else:
            return None

    def getIndexedAtoms(self):
        """
        returns the list of Atom objects with an index
        """
        index_atoms = []
        for atom in self._graph.nodes():
            if atom.index is not None:
                index_atoms.append([atom.index, atom])
        return [atom for [idx, atom] in sorted(index_atoms)]

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

class AtomChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching one labeled atom.

    """
    def __init__(self, AtomInfo = [None, None]):
        """Initialize a chemical environment corresponding to matching a single atom.

        Parameters
        -----------
        AtomInfo: list of lists, optional
            Comes in the form [AtomORtypes, AtomANDtypes]
            AtomORtypes: descriptors for the first atom that are connected with logical operation OR
            AtomANDtypes: descriptors for the first atom that are connected with the logical operation AND

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
    def __init__(self, Atom1Info = [None, None],
            BondInfo = [None, None],
            Atom2Info = [None, None]):
        """Initialize a chemical environment corresponding to matching two atoms (bond).

        Parameters
        -----------
        Atom1Info, Atom2Info: list of lists, optional
            Comes in the form [AtomORtypes, AtomANDtypes]
            AtomORtypes: descriptors for the first atom that are connected with logical operation OR
            AtomANDtypes: descriptors for the first atom that are connected with the logical operation AND
        BondInfo: list of lists, optional
            In the form [BondORtypes, BondANDtypes] similar to atom information

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
        self.atom2 = self.addAtom(self.atom1, BondInfo[0], BondInfo[1], Atom2Info[0], Atom2Info[1], 2)

class AngleChemicalEnvironment(BondChemicalEnvironment):
    """Chemical environment matching three marked atoms (angle).
    """
    def __init__(self, Atom1Info = [None, None], Bond1Info = [None, None],
            Atom2Info = [None, None], Bond2Info = [None, None], Atom3Info = [None, None]):

        """Initialize a chemical environment corresponding to matching three atoms.

        Parameters
        -----------
        Atom1Info, Atom2Info, Atom3Info: list of lists, optional
            Comes in the form [AtomORtypes, AtomANDtypes]
            AtomORtypes: descriptors for the first atom that are connected with logical operation OR
            AtomANDtypes: descriptors for the first atom that are connected with the logical operation AND
        Bond1Info and Bond2Info: list of lists, optional
            In the form [BondORtypes, BondANDtypes] similar to atom information

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
        Atom1Info, Atom2Info, Atom3Info, Atom4Info: list of lists, optional
            Comes in the form [AtomORtypes, AtomANDtypes]
            AtomORtypes: descriptors for the first atom that are connected with logical operation OR
            AtomANDtypes: descriptors for the first atom that are connected with the logical operation AND
        Bond1Info and Bond2Info, Bond3Info: list of lists, optional
            In the form [BondORtypes, BondANDtypes] similar to atom information

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
        Atom1Info, Atom2Info, Atom3Info, Atom4Info: list of lists, optional
            Comes in the form [AtomORtypes, AtomANDtypes]
            AtomORtypes: descriptors for the first atom that are connected with logical operation OR
            AtomANDtypes: descriptors for the first atom that are connected with the logical operation AND
        Bond1Info and Bond2Info, Bond3Info: list of lists, optional
            In the form [BondORtypes, BondANDtypes] similar to atom information

        For example:

        """
        # TODO: add improper example after talking to Christopher about numbering
        # Initialize base class
        super(ImproperChemicalEnvironment,self).__init__(Atom1Info, Bond1Info,
                Atom2Info, Bond2Info, Atom3Info)

        # Add initial atom
        self.atom4 = self.addAtom(self.atom2, Bond3Info[0], Bond3Info[1], Atom4Info[0], Atom4Info[1], 4)
