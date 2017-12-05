#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Representation of molecular topologies.

TODO:
* This will be moved to a more permanent home.

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

import lxml.etree as etree

from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology
from openforcefield.utils import generateTopologyFromOEMol, get_data_filename

import os
import math
import copy
import re
import numpy
import random

from openforcefield.typing.chemistry import ChemicalEnvironment, SMIRKSParsingError

from simtk import openmm, unit

import time

import networkx

import itertools

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================


#=============================================================================================
# TOPOLOGY OBJECTS
#=============================================================================================

class Atom(object):
    """
    Atom representation

    TODO: Should Atom be immutable?

    Parameters
    ----------
    name : str
        An arbitrary label assigned to the atom
    element : str
        Element name

    .. todo::

       Do we want to support the addition of arbitrary additional properties?

    Properties
    ----------
    atomic_number : int
        The atomic number
    mass : float
        The atomic mass

    TODO: Do we want to support the following properties, which would require an Atom to belong to only one Molecule/Residue/Topology?

    bonds : list of Bond
        List of bonds this atom is involved in
    bonded_to : list of Atom
        List of atoms this atom is bonded to
    molecule : Molecule
        The molecule this atom belongs to, if any

    TODO: Should Atoms be mutable or immutable?

    """
    def __init__(self, other=None):
        """
        Create an Atom object.

        Parameters
        ----------
        other : Atom, optional, default=None
            If an Atom is specified, a deep copy will be generated.
        """
        pass

class Bond(object):
    """
    Chemical bond representation

    TODO: Should Bond be immutable?

    Properties
    ----------
    atom1, atom2 : Atom
        Atoms involved in the bond
    bondtype : int
        Discrete bond type representation for the Open Forcefield aromaticity model
        TODO: Do we want to pin ourselves to a single standard aromaticity model?
    order : float
        Fractional bond order

    """
    def __init__(self):
        pass

class ChemicalEntity(object):
    """
    Mixin class for properties shared by chemical entities containing more than one atom.

    Properties
    ----------
    atoms : list of Atom
        Iterate over all Atom objects in the molecule
        TODO:
        * Should we iterate over all atoms in hierarchical order (chains,residues,atoms) or in file order?
        * How can we select different iteration orders?
    bonds :
        Iterate over all Bond objects in the molecule
    angles : tuple of Atom
        Iterate over all angles (Atom tuples) in the molecule
        TODO: Do we need to return an Angle object that collects information about fractional bond orders?
    torsions :
        Iterate over all torsions (propers and impropers) in the molecule
        TODO: Do we need to return a Torsion object that collects information about fractional bond orders?
    propers :
        TODO: Do we need to return a Torsion object that collects information about fractional bond orders?
        Iterate over all proper torsions in the molecule
    impropers :
        Iterate over all improper torsions in the molecule
        TODO: Do we need to return a Torsion object that collects information about fractional bond orders?

    """

class Residue(ChemicalEntity):
    """
    Polymeric residue object

    Properties
    ----------
    residues :
    molecules :

    """
    pass

class Chain(ChemicalEntity):
    """
    Polymeric chain container representation

    (May contain more than one molecule)

    Properties
    ----------
    residues :
    molecules :

    """
    pass

class Molecule(ChemicalEntity):
    """
    Chemical representation of a molecule.

    Properties
    ----------

    These properties are also present in Topology:

    residues : list of Residue
        Iterate over all Residue objects in the molecule
    atoms : list of Atom
        Iterate over all Atom objects in the molecule
        TODO:
        * Should we iterate over all atoms in hierarchical order (chains,residues,atoms) or in file order?
        * How can we select different iteration orders?
    bonds :
        Iterate over all Bond objects in the molecule
    angles : tuple of Atom
        Iterate over all angles (Atom tuples) in the molecule
        TODO: Do we need to return an Angle object that collects information about fractional bond orders?
    torsions :
        Iterate over all torsions (propers and impropers) in the molecule
        TODO: Do we need to return a Torsion object that collects information about fractional bond orders?
    propers :
        TODO: Do we need to return a Torsion object that collects information about fractional bond orders?
        Iterate over all proper torsions in the molecule
    impropers :
        Iterate over all improper torsions in the molecule
        TODO: Do we need to return a Torsion object that collects information about fractional bond orders?

    """
    def __init__(self, other=None):
        """
        Parameters
        ----------
        other : optional, default=None
            If specified, attempt to construct a copy of the Molecule from the specified object.
            This might be a Molecule object, or a file that can be used to construct a Molecule object
            or serialized Molecule object.
        """
        pass

    @staticmethod
    def from_rdkit(rdmol):
        """
        Create a Molecule from an RDKit molecule.

        Requires the RDKit to be installed.

        Parameters
        ----------
        rdmol : rkit.RDMol
            An RDKit molecule

        Returns
        -------
        molecule : openforcefield.Molecule
            An openforcefield molecule

        """
        pass

    def to_rdkit(self):
        """
        Create an RDKit molecule

        Requires the RDKit to be installed.

        Returns
        -------
        rdmol : rkit.RDMol
            An RDKit molecule
        """
        pass

    @staticmethod
    def from_openeye(oemol):
        """
        Create a Molecule from an OpenEye molecule.

        Requires the OpenEye toolkit to be installed.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule

        Returns
        -------
        molecule : openforcefield.Molecule
            An openforcefield molecule

        """
        pass

    def to_openeye(self):
        """
        Create an OpenEye molecule

        Requires the OpenEye toolkit to be installed.

        Returns
        -------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule

        """
        pass

class Topology(ChemicalEntity):
    """
    Chemical representation of a system containing one or more molecules.

    Properties
    ----------

    TODO:
    * Should these properties return deepcopy lists, generators that yield mutable objects, or allow direct mutable access via indexing?
    * Should these be properties or functions?

    chains : list of Chain
        Iterate over all Chain objects in the topology
    molecules : list of Molecule
        Iterate over all Molecule objects in the system in the topology
    unique_molecules : list of Molecule
        Iterate over all unique Molecule objects in the topology

    These properties are also present in Molecule:

    residues : list of Residue
        Iterate over all Residue objects in the topology
    atoms : list of Atom
        Iterate over all Atom objects in the topology
        TODO:
        * Should we iterate over all atoms in hierarchical order (chains,residues,atoms) or in file order?
        * How can we select different iteration orders?
    bonds :
        Iterate over all Bond objects in the topology
    angles : tuple of Atom
        Iterate over all angles (Atom tuples) in the topology
        TODO: Do we need to return an Angle object that collects information about fractional bond orders?
    torsions :
        Iterate over all torsions (propers and impropers) in the topology
        TODO: Do we need to return a Torsion object that collects information about fractional bond orders?
    propers :
        TODO: Do we need to return a Torsion object that collects information about fractional bond orders?
        Iterate over all proper torsions in the topology
    impropers :
        Iterate over all improper torsions in the topology
        TODO: Do we need to return a Torsion object that collects information about fractional bond orders?

    Examples
    --------


    """
    def __init__(self, other=None):
        """
        Parameters
        ----------
        other : optional, default=None
            If specified, attempt to construct a copy of the Topology from the specified object.
            This might be a Topology object, or a file that can be used to construct a Topology object
            or serialized Topology object.

        """
        pass

    @staticmethod
    def from_openmm(openmm_topology, molecules):
        """
        Construct an openforcefield Topology object from an OpenMM Topology object.

        Parameters
        ----------
        openmm_topology : simtk.openmm.app.Topology
            An OpenMM Topology object
        reference_molecules : list of openeye.oechem.OEMol or rdkit.RDMol
            The list of reference molecules in the Topology.

        Returns
        -------
        topology : openforcefield.Topology
            An openforcefield Topology object
        """
        pass

    def to_openmm(self):
        """
        Create an OpenMM Topology object.

        Parameters
        ----------
        openmm_topology : simtk.openmm.app.Topology
            An OpenMM Topology object
        """
        pass

    @staticmethod
    def from_mdtraj(mdtraj_topology, molecules):
        """
        Construct an openforcefield Topology object from an MDTraj Topology object.

        Parameters
        ----------
        mdtraj_topology : mdtraj.Topology
            An MDTraj Topology object
        reference_molecules : list of openeye.oechem.OEMol or rdkit.RDMol
            The list of reference molecules in the Topology.

        Returns
        -------
        topology : openforcefield.Topology
            An openforcefield Topology object
        """
        pass

    def to_mdtraj(self):
        """
        Create an MDTraj Topology object.

        Returns
        ----------
        mdtraj_topology : mdtraj.Topology
            An MDTraj Topology object
        """
        pass

    def is_bonded(self, i, j):
        """Returns True of two atoms are bonded

        Parameters
        ----------
        i, j : int or Atom
            Atoms or atom indices to check

        Returns
        -------
        is_bonded : bool
            True if atoms are bonded, False otherwise.

        """
        pass

    def chemical_environment_matches(self, query):
        """Retrieve all matches for a given chemical environment query.

        TODO:
        * Do we want to generalize this to other kinds of queries too, like mdtraj DSL, pymol selections, atom index slices, etc?
        We could just call it topology.matches(query)

        Parameters
        ----------
        query : str or ChemicalEnvironment
            SMARTS string (with tagged atoms) or ChemicalEnvironment query

        Returns
        -------
        matches : list of Atom tuples
            A list of all matching Atom tuples

        """
        pass
