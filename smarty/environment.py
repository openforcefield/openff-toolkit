#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
environment.py

Classes defining a chemical environment for atoms and how they are connected using SMARTS and SMIRKS. 

AUTHORS

Caitlin Bannan <bannanc@uci.edu>, Mobley Lab, University of California Irvine.

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import networkx

class ChemicalEnvironment(object):
   """Chemical environment abstract base class that matches an atom, bond, angle, etc.
   """
   class Atom(object):
      """Atom representation, which may have some base (logical OR) and decorator (logical AND) properties.

      Type will be OR(basetypes) AND AND(decorators).

      Properties
      -----------
      basetypes : set of string
         The base types that will be combined with logical OR
      decorators : set of string
         The decorators that will be combined with logical AND
      """
      def __init__(self):
         """Initialize an Atom object with empty basetypes and decorators.
         """
         self.basetypes = set()
         self.decorators = set()

      def asSMARTS(self, index=None):
         """Return the atom representation as SMARTS/SMIRKS.

         Parameters
         ----------
         index : int, optional, default=None
            If not None, the specified index will be attached as a SMIRKS index (e.g. '[#6:1]')

         Returns
         --------
         smarts : str
             The SMARTS/SMIRKS string
         """
         pass

   class Bond(object):
      """Bond representation, which may have base (OR) and decorator (AND) types.
      """
      # implementation similar to Atom
      pass

   def __init__(self):
      """Initialize a chemical environment abstract base class.
      """
      # Create an empty graph which will store Atom objects.
      self._graph = nx.Graph()

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

   def addAtom(self, atom, bond, bondedAtom):
      """Add an atom to the specified target atom.
      """
      pass

   def removeAtom(self, atom):
      """Remove the specified atom from the chemical environment, along with its associated bond.
      """
      pass

class AtomChemicalEnvironment(ChemicalEnvironment):
      """Chemical environment matching a single atom.
      """
      def __init__(self):
         """Initialize a chemical environment corresponding to matching a single atom.
         """
         # Initialize base class
         super(AtomChemicalEnvironment,self).__init__()
         # Add a labeled atom corresponding to :1
         # TODO

# and so on for BondChemicalEnvironment, AngleChemicalEnvironment, TorsionChemicalEnvironment
