#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Representation of chemical bonds.

.. todo::

   * Make all classes (like Particle, Atom, VirtualSite) hashable and serializable
   * Use class boilerplate suggestion from Kyle?

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import os
import re
import time
import math
import copy
import string
import random
import itertools
import collections

import lxml.etree as etree

import numpy

import networkx as nx

from simtk import openmm, unit
from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology

from openforcefield.utils import generateTopologyFromOEMol, get_data_filename
from openforcefield.typing.chemistry import ChemicalEnvironment, SMIRKSParsingError

#=============================================================================================
# GLOBAL PARAMETERS
#=============================================================================================

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

class Bond(object):
    """
    Chemical bond representation.

    Attributes
    ----------
    atom1, atom2 : Atom
        Atoms involved in the bond
    bondtype : int
        Discrete bond type representation for the Open Forcefield aromaticity model
        TODO: Do we want to pin ourselves to a single standard aromaticity model?
    type : str
        String based bond type
    order : int
        Integral bond order
    fractional_bondorder : float, optional
        Fractional bond order, or None.

    """
    def __init__(self, atom1, atom2, bondtype, fractional_bondorder=None):
        """
        Create a new chemical bond.
        """
        # TODO: Make sure atom1 and atom2 are both Atom types
        self._atom1 = atom1
        self._atom2 = atom2
        self._type = bondtype
        self._fractional_bondorder = fractional_bondorder

    # TODO: add getters for each of these bond properties

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    @property
    def atoms(self):
        return (self._atom1, self._atom2)

    def type(self):
        return self._type

    @property
    def fractional_bondorder(self):
        return self._fractional_bondorder

    @fractional_bondorder.setter
    def fractional_bondorder(self, value):
        self._fractional_bondorder = value
