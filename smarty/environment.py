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

    def __init__(self, atomBaseSet, atomDecorSet, 
            bondBaseSet = ['~','-',':','=','#'], BondDecorSet = ['@']):
        """Initialize a chemical environment abstract base class.

        ARGUMENTS

        """
        # Create an empty graph which will store Atom objects.
        self._graph = nx.Graph()

        # Storing sets of information
        self.atomBaseSet = set(atomBaseSet)
        self.atomDecorSet = set(atomDecorSet)
        self.bondBaseSet = set(bondBaseSet)
        self.bondDecorSet = set(bondDecorSet)

    def asSMARTS(self):
        """Return a SMARTS/SMIRKS representation of the chemical environment.
        """
        pass

    def selectAtom(self):
        """Select an atom with uniform probability.
        """
        pass

    def selectBond(self):
        """Select a bond with uniform probability.
        """
        pass

    def addAtom(self, atom, bondBases = None, bondDecorators = None, 
            newIndex = None, newBases = None, newDecorators = None):
        """Add an atom to the specified target atom.
        """
        # create new bond
        newBond = self.Bond(bondBases, bondDecorators)
        # create new atom
        newAtom = self.Atom(newAtomIndex, newBases, newDecorators)

        # Add node for newAtom
        self._graph.add_node()
        
        # Connect original atom and new atom
        # HERE
        self._graph.add_edge(atom, newAtom, newBond)

    def removeAtom(self, atom):
        """Remove the specified atom from the chemical environment, along with its associated bond.
        """
        pass

class AtomChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching one labeled atom.
    """
    def __init__(self, atomBaseSet, atomDecorSet, 
            bondBaseSet = ['~','-',':','=','#'], BondDecorSet = ['@']):
        """Initialize a chemical environment corresponding to matching a single atom.
        """
        # Initialize base class
        super(AtomChemicalEnvironment,self).__init__()
        # Add a labeled atom corresponding to :1
        # TODO

        pass

class BondChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching two labeled atoms (or a bond).
    """
    def __init__(self, atomBaseSet, atomDecorSet, 
            bondBaseSet = ['~','-',':','=','#'], BondDecorSet = ['@']):
        """Initialize a chemical environment corresponding to matching two atoms (bond).
        """
        # Initialize base class
        super(AtomChemicalEnvironment,self).__init__()
        # Add a labeled atom corresponding to :1
        # TODO
        pass

class AngleChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching three marked atoms (angle).
    """
    def __init__(self, atomBaseSet, atomDecorSet, 
            bondBaseSet = ['~','-',':','=','#'], BondDecorSet = ['@']):
        """Initialize a chemical environment corresponding to matching three atoms.
        """
        # Initialize base class
        super(AtomChemicalEnvironment,self).__init__()
        # Add a labeled atom corresponding to :1
        # TODO
        pass

class TorsionChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching four marked atoms (torsion).
    """
    def __init__(self, atomBaseSet, atomDecorSet, 
            bondBaseSet = ['~','-',':','=','#'], BondDecorSet = ['@']):
        """Initialize a chemical environment corresponding to matching four atoms.
        """
        # Initialize base class
        super(AtomChemicalEnvironment,self).__init__()
        # Add a labeled atom corresponding to :1
        # TODO
        pass

class ImproperChemicalEnvironment(ChemicalEnvironment):
    """Chemical environment matching four marked atoms (improper).
    """
    def __init__(self, atomBaseSet, atomDecorSet, 
            bondBaseSet = ['~','-',':','=','#'], BondDecorSet = ['@']):
        """Initialize a chemical environment corresponding to matching four atoms
            connected as an improper torsion.
        """
        # Initialize base class
        super(AtomChemicalEnvironment,self).__init__()
        # Add a labeled atom corresponding to :1
        # TODO
        pass

