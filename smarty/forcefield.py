#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
forcefield.py

OpenMM ForceField replacement using SMIRKS-based matching.

AUTHORS

John D. Chodera <john.chodera@choderalab.org>

Baseed on simtk.openmm.app.forcefield written by Peter Eastman.

TODO
* Add general constraint syntax using SMIRKS-based matching:
```XML
<Constraints>
  <!-- bonds involving hydrogen -->
  <Constraint smirks="[#1:1]-[*:2]"/> <!-- add constraint between atoms 1 and 2 using equilibrium bond distance -->
  <!-- angles involving hydrogen -->
  <Constraint smirks="[#1:1]-[*:2]-[#1:3]"/> <!-- add constraint between atoms 1 and 3 using auto-calculated distance from equilibrium bond and angles -->
</Constraints>
```
* Move utility functions like 'generateTopologyFromOEMol()' elsewhere?
"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

#import xml.etree.ElementTree as etree
import lxml.etree as etree

from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology

import os
import math
import copy
import re
import numpy
import random

import openeye.oechem
import openeye.oeomega
import openeye.oequacpac

from openeye import oechem

from simtk import openmm, unit

import time

import networkx

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

def _convertParameterToNumber(param):
    if unit.is_quantity(param):
        if param.unit.is_compatible(unit.bar):
            return param / unit.bar
        return param.value_in_unit_system(unit.md_unit_system)
    return float(param)

#=============================================================================================
# Augmented Topology
#=============================================================================================

def generateTopologyFromOEMol(molecule):
    """
    Generate an OpenMM Topology object from an OEMol molecule.

    Parameters
    ----------
    molecule : openeye.oechem.OEMol
        The molecule from which a Topology object is to be generated.

    Returns
    -------
    topology : simtk.openmm.app.Topology
        The Topology object generated from `molecule`.

    """
    # Create a Topology object with one Chain and one Residue.
    from simtk.openmm.app import Topology
    topology = Topology()
    chain = topology.addChain()
    resname = molecule.GetTitle()
    residue = topology.addResidue(resname, chain)

    # Create atoms in the residue.
    for atom in molecule.GetAtoms():
        name = atom.GetName()
        element = elem.Element.getByAtomicNumber(atom.GetAtomicNum())
        atom = topology.addAtom(name, element, residue)

    # Create bonds.
    atoms = { atom.name : atom for atom in topology.atoms() }
    for bond in molecule.GetBonds():
        topology.addBond(atoms[bond.GetBgn().GetName()], atoms[bond.GetEnd().GetName()])

    return topology

def generateGraphFromTopology(topology):
    """Geneate a NetworkX graph from a Topology object.

    Parameters
    ----------
    topology : simtk.openmm.app.Topology
        The source topology.

    Returns
    -------
    graph : networkx.Graph
        The resulting graph, with nodes labeled with atom indices and elements

    """
    import networkx as nx
    # Create graph of atoms connected by bonds.
    G = nx.Graph()
    for atom in topology.atoms():
        G.add_node(atom.index, element=atom.element)
    for (atom1, atom2) in topology.bonds():
        G.add_edge(atom1.index, atom2.index)

    return G

class _Topology(Topology):
    """Augmented Topology object which adds:

    self._reference_molecules is a list of OEMol for the reference molecules
    self._reference_to_topology_atom_mappings[reference_molecule] is a list of atom indices mapping a reference molecule atom index to the topology atom index
    """
    def __init__(self, topology, reference_molecules):
        """
        Parameters
        ----------
        topology : simtk.openmm.app.Topology
            The Topology object to initialize this one from.
        reference_molecules : list of openeye.oechem.OEMol
            The list of reference molecules in the Topology.

        """
        # Initialize.
        super(_Topology, self).__init__()

        # TODO: Find a way to avoid having this be fragile based on internal representation of Topology.
        # TODO: Should this also use a deepcopy of 'topology' first?
        self._chains = topology._chains
        self._numResidues = topology._numResidues
        self._numAtoms = topology._numAtoms
        self._bonds = topology._bonds
        self._periodicBoxVectors = topology._periodicBoxVectors

        # Store reference molecules.
        # TODO: Deep copy?
        self._reference_molecules = reference_molecules

        # Identify all molecules and atom mappings.
        self._identifyMolecules()

    def _identifyMolecules(self):
        """Identify all unique reference molecules and atom mappings to all instances in the Topology.
        """
        import networkx as nx
        from networkx.algorithms import isomorphism

        # Generate list of topology atoms.
        atoms = [ atom for atom in self.atoms() ]

        # Generate graphs for reference molecules.
        self._reference_molecule_graphs = list()
        for reference_molecule in self._reference_molecules:
            # Generate Topology
            reference_molecule_topology = generateTopologyFromOEMol(reference_molecule)
            # Generate Graph
            reference_molecule_graph = generateGraphFromTopology(reference_molecule_topology)
            self._reference_molecule_graphs.append(reference_molecule_graph)

        # Generate a graph for the current topology.
        G = generateGraphFromTopology(self)

        # Extract molecules (as connected component subgraphs).
        self._reference_to_topology_atom_mappings = { reference_molecule : list() for reference_molecule in self._reference_molecules }
        for molecule_graph in nx.connected_component_subgraphs(G):
            # Check if we have already stored a reference molecule for this molecule.
            reference_molecule_exists = False
            for (reference_molecule_graph, reference_molecule) in zip(self._reference_molecule_graphs, self._reference_molecules):
                GM = isomorphism.GraphMatcher(molecule_graph, reference_molecule_graph)
                if GM.is_isomorphic():
                    # This molecule is present in the list of unique reference molecules.
                    reference_molecule_exists = True
                    # Add the reference atom mappings.
                    reference_to_topology_atom_mapping = dict()
                    for (topology_atom, reference_atom) in GM.mapping.iteritems():
                        reference_to_topology_atom_mapping[reference_atom] = topology_atom
                    self._reference_to_topology_atom_mappings[reference_molecule].append(reference_to_topology_atom_mapping)
                    # Break out of the search loop.
                    break

            # If the reference molecule could not be found, throw an exception.
            if not reference_molecule_exists:
                msg = 'No provided molecule matches topology molecule:\n'
                for index in sorted(list(molecule_graph)):
                    msg += 'Atom %8d %5s %5d %3s\n' % (atoms[index].index, atoms[index].name, atoms[index].residue.index, atoms[index].residue.name)
                raise Exception(msg)

    def getSMIRKSMatches(self, smirks):
        """Find all sets of atoms in the topology that match the provided SMIRKS strings.

        Parameters
        ----------
        smirks : str
            SMIRKS string with tagged atoms.
            If there are N tagged atoms numbered 1..N, the resulting matches will be N-tuples of atoms that match the corresponding tagged atoms.

        Returns
        -------
        matches : list of tuples of Atom
            matches[index] is an N-tuple of Atom entries from the topology
            Matches are returned in no guaranteed order.

        """
        # Set up query.
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, smirks):
            raise Exception("Error parsing SMIRKS '%s'" % smirks)

        # Perform matching on each unique molecule, unrolling the matches to all matching copies of tha tmolecule in the Topology object.
        matches = list()
        unique = True # give unique matches
        for reference_molecule in self._reference_molecules:
            # Find all atomsets that match this definition in the reference molecule
            ss = oechem.OESubSearch(qmol)
            for match in ss.Match(reference_molecule, unique):
                # Compile list of reference atom indices that match the pattern tags.
                reference_atom_indices = dict()
                for ma in match.GetAtoms():
                    if ma.pattern.GetMapIdx() != 0:
                        reference_atom_indices[ma.pattern.GetMapIdx()-1] = ma.target.GetIdx()
                # Compress into list.
                reference_atom_indices = [ reference_atom_indices[index] for index in range(len(reference_atom_indices)) ]
                # Unroll all instances of this molecule.
                for reference_to_topology_atom_mapping in self._reference_to_topology_atom_mappings[reference_molecule]:
                    # Create match.
                    atom_indices = tuple([ reference_to_topology_atom_mapping[atom_index] for atom_index in reference_atom_indices ])
                    matches.append(atom_indices)

        return matches

#=============================================================================================
# FORCEFIELD
#=============================================================================================

# Enumerated values for nonbonded method

class NoCutoff(object):
    def __repr__(self):
        return 'NoCutoff'
NoCutoff = NoCutoff()

class CutoffNonPeriodic(object):
    def __repr__(self):
        return 'CutoffNonPeriodic'
CutoffNonPeriodic = CutoffNonPeriodic()

class CutoffPeriodic(object):
    def __repr__(self):
        return 'CutoffPeriodic'
CutoffPeriodic = CutoffPeriodic()

class Ewald(object):
    def __repr__(self):
        return 'Ewald'
Ewald = Ewald()

class PME(object):
    def __repr__(self):
        return 'PME'
PME = PME()

# Enumerated values for constraint type

class HBonds(object):
    def __repr__(self):
        return 'HBonds'
HBonds = HBonds()

class AllBonds(object):
    def __repr__(self):
        return 'AllBonds'
AllBonds = AllBonds()

class HAngles(object):
    def __repr__(self):
        return 'HAngles'
HAngles = HAngles()

# A map of functions to parse elements of the XML file.

parsers = {}

class ForceField(object):
    """A ForceField constructs OpenMM System objects based on a Topology.
    """

    def __init__(self, *files):
        """Load one or more XML parameter definition files and create a SMIRFF ForceField object based on them.

        Parameters
        ----------
        files : list
            A list of XML files defining the SMIRFF force field.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in force fields), or an open file-like object with a read() method from which the forcefield XML data can be loaded.

        """
        self._forces = []
        self.loadFile(files)

    def loadFile(self, files):
        """Load a SMIRFF XML file and add the definitions from it to this ForceField.

        Parameters
        ----------
        files : string or file or tuple
            An XML file or tuple of XML files containing SMIRFF force field definitions.
            Each entry may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in force fields), or an open file-like object with a read() method from which the forcefield XML data can be loaded.
        """

        # Ensure that we are working with a tuple of files.
        if not isinstance(files, tuple):
            files = (files,)

        # Load in all XML trees.
        trees = list()
        for file in files:
            try:
                # this handles either filenames or open file-like objects
                tree = etree.parse(file)
            except IOError:
                tree = etree.parse(os.path.join(os.path.dirname(__file__), 'data', file))
            except Exception as e:
                # Fail with an error message about which file could not be read.
                # TODO: Also handle case where fallback to 'data' directory encounters problems,
                # but this is much less worrisome because we control those files.
                msg  = str(e) + '\n'
                if hasattr(file, 'name'):
                    filename = file.name
                else:
                    filename = str(file)
                msg += "ForceField.loadFile() encountered an error reading file '%s'\n" % filename
                raise Exception(msg)

            trees.append(tree)

        # Load the atom masses.
        for tree in trees:
            if tree.getroot().find('AtomTypes') is not None:
                for type in tree.getroot().find('AtomTypes').findall('Type'):
                    self.registerAtomType(type.attrib)

        # Load force definitions
        for tree in trees:
            for child in tree.getroot():
                if child.tag in parsers:
                    parsers[child.tag](child, self)

    def getGenerators(self):
        """Get the list of all registered generators."""
        return self._forces

    def registerGenerator(self, generator):
        """Register a new generator."""
        self._forces.append(generator)

    def createSystem(self, topology, molecules, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*unit.nanometer,
                     constraints=None, rigidWater=True, removeCMMotion=True, hydrogenMass=None, residueTemplates=dict(), verbose=False, **kwargs):
        """Construct an OpenMM System representing a Topology with this force field.

        Parameters
        ----------
        topology : Topology
            The Topology for which to create a System
        molecules : list of openeye.oechem.OEMol
            List of molecules appearing in the topology
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, or PME.
        nonbondedCutoff : distance=1*nanometer
            The cutoff distance to use for nonbonded interactions
        constraints : object=None
            Specifies which bonds and angles should be implemented with constraints.
            Allowed values are None, HBonds, AllBonds, or HAngles.
        rigidWater : boolean=True
            If true, water molecules will be fully rigid regardless of the value
            passed for the constraints argument
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep
            their total mass the same.
        residueTemplates : dict=dict()
           Key: Topology Residue object
           Value: string, name of _TemplateData residue template object to use for
                  (Key) residue
           This allows user to specify which template to apply to particular Residues
           in the event that multiple matching templates are available (e.g Fe2+ and Fe3+
           templates in the ForceField for a monoatomic iron ion in the topology).
        verbose : bool
           If True, verbose output will be printed.
        kwargs
             Arbitrary additional keyword arguments may also be specified.
             This allows extra parameters to be specified that are specific to
             particular force fields.

        Returns
        -------
        system
            the newly created System
        """
        # Work with a modified form of the topology that provides additional accessors.
        topology = _Topology(topology, molecules)

        # Create the System and add atoms
        system = openmm.System()
        for atom in topology.atoms():
            # Add the particle to the OpenMM system.
            system.addParticle(atom.element.mass) # TODO: Add option to use a different mass than the integral elemental mass?

        # Adjust hydrogen masses if requested.
        if hydrogenMass is not None:
            if not unit.is_quantity(hydrogenMass):
                hydrogenMass *= unit.dalton
            for atom1, atom2 in topology.bonds():
                if atom1.element == elem.hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == elem.hydrogen and atom1.element not in (elem.hydrogen, None):
                    transferMass = hydrogenMass-system.getParticleMass(atom2.index)
                    system.setParticleMass(atom2.index, hydrogenMass)
                    system.setParticleMass(atom1.index, system.getParticleMass(atom1.index)-transferMass)

        # Set periodic boundary conditions.
        boxVectors = topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            system.setDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
        elif nonbondedMethod not in [NoCutoff, CutoffNonPeriodic]:
            raise ValueError('Requested periodic boundary conditions for a Topology that does not specify periodic box dimensions')

        # TODO: Convert requested bonds and angles to use constraints

        # Set nonbonded method.
        kwargs['nonbondedMethod'] = nonbondedMethod
        kwargs['nonbondedCutoff'] = nonbondedCutoff

        # Add forces to the System
        for force in self._forces:
            force.createForce(system, topology, verbose=verbose, **kwargs)
        if removeCMMotion:
            system.addForce(openmm.CMMotionRemover())

        # Let force generators do postprocessing
        for force in self._forces:
            if 'postprocessSystem' in dir(force):
                force.postprocessSystem(system, data, **kwargs)

        return sys

#=============================================================================================
# The following classes are generators that know how to create Force subclasses and add them to a System that is being
# created.  Each generator class must define two methods: 1) a static method that takes an etree Element and a ForceField,
# and returns the corresponding generator object; 2) a createForce() method that constructs the Force object and adds it
# to the System.  The static method should be added to the parsers map.
#=============================================================================================

def _validateSMIRKS(smirks, node=None):
    """Validate the specified SMIRKS string.

    Parameters
    ----------
    smirks : str
       The SMIRKS string to be validated
    node : xml.etree.ElementTree.Element
       Node of etree with 'sourceline' attribute.

    """
    qmol = oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smirks):
        if (node is not None) and ('sourceline' in note.attrib):
            raise Exception("Line %s: Error parsing SMIRKS '%s'" % (node.attrib['sourceline'], node.attrib['smirks']))
        else:
            raise Exception("Error parsing SMIRKS '%s'" % (node.attrib['smirks']))

    return smirks

import collections
class TransformedDict(collections.MutableMapping):
    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key

class ValenceDict(TransformedDict):
    """Enforce uniqueness in atom indices"""
    def __keytransform__(self, key):
        """Reverse tuple if first element is larger than last element."""
        # Ensure key is a tuple.
        key = tuple(key)
        # Reverse the key if the first element is bigger than the last.
        if key[0] > key[-1]:
            key = tuple(reversed(key))
        return key

#=============================================================================================
# Force generators
#=============================================================================================

## @private
class HarmonicBondGenerator(object):
    """A HarmonicBondGenerator constructs a HarmonicBondForce."""

    class BondType(object):
        """A SMIRFF bond type."""
        def __init__(self, node):
            self.smirks = _validateSMIRKS(node.attrib['smirks'], node=node)
            self.length = _convertParameterToNumber(node.attrib['length'])
            self.k = _convertParameterToNumber(node.attrib['k'])

    def __init__(self, forcefield):
        self.ff = forcefield
        self._bondtypes = list()

    def registerBond(self, node):
        """Register a SMIRFF bondtype definition."""
        bond = HarmonicBondGenerator.BondType(node)
        self._bondtypes.append(bond)

    @staticmethod
    def parseElement(element, ff):
        # Find existing force generator or create new one.
        existing = [f for f in ff._forces if isinstance(f, HarmonicBondGenerator)]
        if len(existing) == 0:
            generator = HarmonicBondGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]

        # Register all SMIRFF bond definitions.
        for bond in element.findall('Bond'):
            generator.registerBond(bond)

    def createForce(self, system, topology, verbose=False, **kwargs):
        # Find existing force or create new one.
        existing = [system.getForce(i) for i in range(system.getNumForces())]
        existing = [f for f in existing if type(f) == openmm.HarmonicBondForce]
        if len(existing) == 0:
            force = openmm.HarmonicBondForce()
            system.addForce(force)
        else:
            force = existing[0]

        # Iterate over all defined bond types, allowing later matches to override earlier ones.
        bonds = ValenceDict()
        for bond in self._bondtypes:
            for atom_indices in topology.getSMIRKSMatches(bond.smirks):
                bonds[atom_indices] = bond

        if verbose:
            print('')
            print('HarmonicBondForceGenerator:')
            print('')
            for bond in self._bondtypes:
                print('%32s : %8d matches' % (bond.smirks, len(topology.getSMIRKSMatches(bond.smirks))))
            print('')

        # Add all bonds to the system.
        for (atom_indices, bond) in bonds.items():
            force.addBond(atom_indices[0], atom_indices[1], bond.length, bond.k)

        if verbose: print('%d bonds added' % (len(bonds)))

        # DEBUG: Check that no topology bonds aren't missing force bonds.
        atoms = [ atom for atom in topology.atoms() ]
        topology_bonds = ValenceDict()
        for (atom1, atom2) in topology.bonds():
            topology_bonds[(atom1.index,atom2.index)] = True
        if set(bonds.keys()) != set(topology_bonds.keys()):
            msg = 'Mismatch between bonds added and topological bonds.\n'
            created_bondset = set(bonds.keys())
            topology_bondset = set(topology_bonds.keys())
            msg += 'Bonds created that are not present in Topology:\n'
            msg += str(created_bondset.difference(topology_bondset)) + '\n'
            msg += 'Topology bonds not assigned parameters:\n'
            for (a1, a2) in topology_bondset.difference(created_bondset):
                atom1 = atoms[a1]
                atom2 = atoms[a2]
                msg += '(%8d,%8d) : %5s %3s %3s - %5s %3s %3s' % (a1, a2, atom1.residue.index, atom1.residue.name, atom1.name, atom2.residue.index, atom2.residue.name, atom2.name)
                msg += '\n'
            raise Exception(msg)

parsers["HarmonicBondForce"] = HarmonicBondGenerator.parseElement

#=============================================================================================
