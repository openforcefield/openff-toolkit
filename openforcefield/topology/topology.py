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
    A chemical atom

    Note that non-chemical virtual sites are represented by the ``VirtualSite`` object

    .. todo::
        * Should ``Atom`` objects be immutable or mutable?
        * Should an ``Atom`` be able to belong to more than one ``Topology`` object?
        * Do we want to support the addition of arbitrary additional properties,
          such as floating point quantities (e.g. ``charge``), integral quantities (such as ``id`` or ``serial`` index in a PDB file),
          or string labels (such as Lennard-Jones types)?
        * Should we be able to create ``Atom`` objects on their own, or only in the context of a ``Topology`` object they belong to?

    """
    def __init__(self, name, element):
        """
        Create an Atom object.

        Parameters
        ----------
        name : str
            A unique name for this atom
        element : str
            The element name

        """
        pass

    @property
    def index(self):
        """
        The unique index of this atom within a Topology

        """
        pass

    @property
    def name(self):
        """
        An arbitrary label assigned to the atom.

        """
        pass

    @property
    def element(self):
        """
        The element name

        """
        pass

    @property
    def atomic_number(self):
        """
        The integer atomic number of the atom.

        """
        pass

    @property
    def mass(self):
        """
        The atomic mass of the atomic site.

        """
        pass

    @property
    def bonds(self):
        """
        The list of ``Bond`` objects this atom is involved in.

        """
        pass

    @property
    def bonded_to(self):
        """
        The list of ``Atom`` objects this atom is involved in

        """
        pass

    @property
    def molecule(self):
        """
        The ``Molecule`` this atom is part of.

        .. todo::
            * Should we have a single unique ``Molecule`` for each molecule type in the system,
            or if we have multiple copies of the same molecule, should we have multiple ``Molecule``s?
        """
        pass

class VirtualSite(object):
    """
    A virtual (non-atom) site.

    .. todo::
        * Should virtual sites be attached to one atom only, or more than one atom?
          OpenMM defines them as belonging to two or more atoms.
        * Should a virtual site be able to belong to more than one Topology?
        * Should virtual sites be immutable or mutable?

    """
    def __init__(self):
        """
        Create a virtual site

        """
        pass

    @property
    def index(self):
        """
        The index of this VirtualSite within a ``Topology``

        """
        pass

class Bond(object):
    """
    Chemical bond representation

    TODO: Should Bond be immutable?

    Attributes
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

    """

    @property
    def atoms(self):
        """
        Iterate over all Atom objects in the molecule

        .. todo::
            * Should we iterate over all atoms in hierarchical order (chains,residues,atoms) or in file order?
            * How can we select different iteration orders?

        """
        pass

    @property
    def bonds(self):
        """
        Iterate over all Bond objects in the molecule

        """
        pass

    def angles(self):
        """
        Iterate over all angles (Atom tuples) in the molecule

        .. todo::
            * Do we need to return an Angle object that collects information about fractional bond orders?

        """
        pass

    @property
    def torsions(self):
        """
        Iterate over all torsions (propers and impropers) in the molecule

        .. todo::
            * Do we need to return a ``Torsion`` object that collects information about fractional bond orders?
            * Should we call this ``dihedrals`` instead of ``torsions``?

        """
        pass

    @property
    def propers(self):
        """
        Iterate over all proper torsions in the molecule

        .. todo::
            * Do we need to return a ``Torsion`` object that collects information about fractional bond orders?
        """
        pass

    @property
    def impropers(self):
        """
        Iterate over all proper torsions in the molecule

        .. todo::
            * Do we need to return a ``Torsion`` object that collects information about fractional bond orders?
        """
        pass

class Residue(ChemicalEntity):
    """
    Polymeric residue object

    Attributes
    ----------
    atoms : list of Atom
        The atoms that belong to this residue
    molecule : Molecule
        The Molecule that this residue belongs to

    """
    pass

class Chain(ChemicalEntity):
    """
    Polymeric chain container representation

    (May contain more than one molecule, as in a PDB chain containing waters)

    .. todo::
        * It seems like the hierarchical view (chains, residues, atoms) is arbitrary and should be simply
          left as annotation added to atoms to allow iterating over them in this hierarchical way if desired,
          rather than as first-class objects

    Attributes
    ----------
    residues : list of Residue
        The residues within this chain
    molecule : list of Molecules
        The molecules associated with this chain

    """
    pass

class Molecule(ChemicalEntity):
    """
    Chemical representation of a molecule.

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

    .. todo::
        * Should these properties return deepcopy lists, generators that yield mutable objects, or allow direct mutable access via indexing?
        * Should these be properties or functions?


    Attributes
    ----------
    chains : list of Chain
        Iterate over all Chain objects in the topology
    molecules : list of Molecule
        Iterate over all Molecule objects in the system in the topology
    unique_molecules : list of Molecule
        Iterate over all unique Molecule objects in the topology


    Examples
    --------
    Create a Topology copy
    >>> topology_copy = Topology(topology)

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
