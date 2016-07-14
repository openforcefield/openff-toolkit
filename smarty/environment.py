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

#=============================================================================================
# AtomEnvironment Object (node of chemical environment)
#=============================================================================================

class AtomEnvironment(object):
    """
    This is an Atom that is a part of a larger chemical environment
    This may be just a dictionary, but I think it is worth having an option to add other methods

    For example, we could eventually have a get valence property that would use the atom, charg

    """
    def __init__(self, base_types = None, properties = None, label = ""):  
        """
        """

        self.base_types = base_types
        self.properties = properties
        self.label = label
        

#=============================================================================================
# Chemical Environment 
#=============================================================================================

class ChemicalEnvironment(object):
    """
    Chemical Environments describe a molecular fragment with
    atoms, their properties, and how they are bonded together

    This includes the properties of those atoms and bonds and the methods 
    associated with changing or adding atoms or bonds to the fragment.
    """

    def __init__(self, connection_type, basetypes, decorators):
        """
        initialize a chemical environment graph with 'empty' atoms for the number of nodes needed
        for the connection type.

        Should we have an option for initial atom/bond identifiers?

        ARGUMENTS

        connection_type: string? number? 
            Determines the type of property being assigned, Lennard-Jones (or atomtype), Bond, Angle, Torsion, or Improper
        basetypes: list of strings             
            List of atom identifiers, this can include element numbers ("#1"), or a list of atoms ("#6,#7,#8")
        decorators: list of strings
            List of properties that can describe an atom, such as "X3" for connectivity 3
        """
        self.connection_type = connection_type
        self.base_types = base_types
        self.decorators = decorators

        self.Graph = self.buildGraph(self.connection_type)

    def buildGraph(self, connection_type):
        """
        Makes a networkx graph that is the minimal size required for the connection_type specified.

        ARGUMENTS
        connection_type: string? number?
            property being determined, atom, bond, angle, torsion, or improper

        Returns a networkx graph object with the assigned shape and atom label numbers assigned. 
        """
