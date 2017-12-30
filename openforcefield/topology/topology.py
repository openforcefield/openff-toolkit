#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Representation of molecular topologies.

.. todo::

   * Make all classes (like Particle, Atom, VirtualSite) hashable
     Use class boilerplate suggestion from Kyle?

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

import lxml.etree as etree

import numpy

import networkx

from simtk import openmm, unit
from simtk.openmm.app import element as elem
from simtk.openmm.app import Topology

from openforcefield.utils import generateTopologyFromOEMol, get_data_filename
from openforcefield.typing.chemistry import ChemicalEnvironment, SMIRKSParsingError

#=============================================================================================
# GLOBAL PARAMETERS
#=============================================================================================

DEFAULT_AROMATICITY_MODEL = 'MDL' # TODO: Is there a more specific spec for this?
DEFAULT_FRACTIONAL_BONDORDER_MODEL = 'Wiberg' # TODO: Is there a more specific spec for this?
DEFAULT_CHARGE_MODEL = 'AM1-BCC' # TODO: Should this be `AM1` so that BCCs can appear in the SMIRNOFF forcefield?

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

def _getSMIRKSMatches_OEMol(oemol, smirks, aromaticity_model=None):
    """Find all sets of atoms in the provided oemol that match the provided SMIRKS strings.

    Parameters
    ----------
    oemol : OpenEye oemol
        oemol to process with the SMIRKS in order to find matches
    smirks : str
        SMIRKS string with tagged atoms.
        If there are N tagged atoms numbered 1..N, the resulting matches will be N-tuples of atoms that match the corresponding tagged atoms.
    aromaticity_model : str (optional)
        OpenEye aromaticity model designation as a string, such as "OEAroModel_MDL". Default: None. If none is provided, molecule is processed exactly as provided; otherwise it is prepared with this aromaticity model prior to querying.

    Returns
    -------
    matches : list of tuples of atoms numbers
        matches[index] is an N-tuple of atom numbers from the oemol
        Matches are returned in no guaranteed order.
    """
    from openeye import oechem

    # Make a copy of molecule so we don't influence original (probably safer than deepcopy per C Bayly)
    mol = oechem.OEMol(oemol)

    # Set up query.
    qmol = oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smirks):
        raise Exception("Error parsing SMIRKS '%s'" % smirks)

    # Determine aromaticity model
    if aromaticity_model:
        if type(aromaticity_model) == str:
            # Check if the user has provided a manually-specified aromaticity_model
            if hasattr(oechem, aromaticity_model):
                oearomodel = getattr(oechem, 'OEAroModel_' + aromaticity_model)
            else:
                raise ValueError("Error: provided aromaticity model not recognized by oechem.")
        else:
            raise ValueError("Error: provided aromaticity model must be a string.")

        # If aromaticity model was provided, prepare molecule
        oechem.OEClearAromaticFlags( mol)
        oechem.OEAssignAromaticFlags( mol, oearomodel)
        # avoid running OEPrepareSearch or we lose desired aromaticity, so instead:
        oechem.OEAssignHybridization( mol)

    # Perform matching on each mol
    matches = list()

    # We require non-unique matches, i.e. all matches
    unique = False
    ss = oechem.OESubSearch(qmol)
    matches = []
    for match in ss.Match( mol, unique):
        # Compile list of atom indices that match the pattern tags
        atom_indices = dict()
        for ma in match.GetAtoms():
            if ma.pattern.GetMapIdx() != 0:
                atom_indices[ma.pattern.GetMapIdx()-1] = ma.target.GetIdx()
        # Compress into list
        atom_indices = [ atom_indices[index] for index in range(len(atom_indices)) ]
        # Store
        matches.append( tuple(atom_indices) )

    return matches

#=============================================================================================
# Augmented Topology
#=============================================================================================

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
    self._reference_to_topology_atom_mappings[reference_molecule] is a list of dicts, where each dict maps the atom indices of atoms in the reference molecule onto an equivalent atom index for a topology atom.
    self._bondorders is a list of floating point bond orders for the bonds in the Topology.
    self._bondorders_by_atomindices is a dict of floating point bond orders for the bonds in the Topology, keyed by indices of the atoms involved.

    Assumes class is immutable.

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

        # Get/initialize bond orders
        self._updateBondOrders()

        # Track constraints
        self._constrained_atom_pairs = dict()

    def angles(self):
        """
        Get an iterator over all i-j-k angles.
        """
        if not hasattr(self, '_angles'):
            self._construct_bonded_atoms_list()
            self._angles = set()
            for atom1 in self._atoms:
                for atom2 in self._bondedAtoms[atom1]:
                    for atom3 in self._bondedAtoms[atom2]:
                        if atom1 == atom3:
                            continue
                        if atom1.index < atom3.index:
                            self._angles.add( (atom1, atom2, atom3) )
                        else:
                            self._angles.add( (atom3, atom2, atom1) )

        return iter(self._angles)

    def torsions(self):
        """
        Get an iterator over all i-j-k-l torsions.
        Note that i-j-k-i torsions are excluded.
        """
        if not hasattr(self, '_torsions'):
            self._construct_bonded_atoms_list()

            self._torsions = set()
            for atom1 in self._atoms:
                for atom2 in self._bondedAtoms[atom1]:
                    for atom3 in self._bondedAtoms[atom2]:
                        if atom1 == atom3:
                            continue
                        for atom4 in self._bondedAtoms[atom3]:
                            if atom4 == atom2:
                                continue
                            # Exclude i-j-k-i
                            if atom1 == atom4:
                                continue
                            if atom1.index < atom4.index:
                                self._torsions.add( (atom1, atom2, atom3, atom4) )
                            else:
                                self._torsions.add( (atom4, atom3, atom2, atom1) )

        return iter(self._torsions)

    def _construct_bonded_atoms_list(self):
        """
        Construct list of all atoms each atom is bonded to.
        """
        if not hasattr(self, '_bondedAtoms'):
            self._atoms = [ atom for atom in self.atoms() ]
            self._bondedAtoms = dict()
            for atom in self._atoms:
                self._bondedAtoms[atom] = set()
            for bond in self._bonds:
                self._bondedAtoms[bond[0]].add(bond[1])
                self._bondedAtoms[bond[1]].add(bond[0])

    def _isBonded(self, atom_index_1, atom_index_2):
        """Return True if atoms are bonded, False if not.

        Parameters
        ----------
        atom_index_1 : int
        atom_index_2 : int
            Atom indices

        Returns
        -------
        is_bonded : bool
            True if atoms are bonded, False otherwise

        TODO
        ----
        This assumes _Topology is immutable.
        """
        self._construct_bonded_atoms_list()
        atom1 = self._atoms[atom_index_1]
        atom2 = self._atoms[atom_index_2]
        return atom2 in self._bondedAtoms[atom1]

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
                    for (topology_atom, reference_atom) in GM.mapping.items():
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

    def _updateBondOrders(self, Wiberg = False):
        """Update and store list of bond orders for the molecules in this Topology. Can be used for initialization of bondorders list, or for updating bond orders in the list.

        Parameters:
        ----------
        Wiberg : bool (optional)
            Default False. If False, uses bond orders OEChem assigns to bonds on the molecule. If True, instead uses Wiberg bond orders stored on bonds in the molecule. These must already be present, i.e. from assignPartialCharges with an AM1 method.

        """
        # Initialize
        self._bondorders=list()
        self._bondorders_by_atomindices = {}
        # Loop over reference molecules and pull bond orders

        for mol in self._reference_molecules:
            # Pull mappings for this molecule
            mappings = self._reference_to_topology_atom_mappings[mol]
            # Loop over bonds
            for idx,bond in enumerate(mol.GetBonds()):
                # Get atom indices involved in bond
                at1 = bond.GetBgn().GetIdx()
                at2 = bond.GetEnd().GetIdx()
                # Get bond order
                if not Wiberg:
                    order = bond.GetOrder()
                else:
                    order = bond.GetData('WibergBondOrder')
                # Convert atom numbers to topology atom numbers; there may be multiple matches
                for mapping in mappings:
                    topat1 = None
                    topat2 = None
                    for mapatom in mapping:
                        if mapatom==at1:
                            topat1 = mapping[mapatom]
                        elif mapatom==at2:
                            topat2 = mapping[mapatom]
                    if topat1==None or topat2==None:
                        raise ValueError("No mapping found for these topology atoms (indices %s-%s)." % (at1, at2))
                    # Store bond order to re-use below and elsewhere; store in both directions
                    if not topat1 in self._bondorders_by_atomindices:
                        self._bondorders_by_atomindices[topat1] = {}
                    if not topat2 in self._bondorders_by_atomindices:
                        self._bondorders_by_atomindices[topat2] = {}
                    self._bondorders_by_atomindices[topat2][topat1] = order
                    self._bondorders_by_atomindices[topat1][topat2] = order

        # Loop over bonds in topology and store orders in the same order
        for bond in self._bonds:
            # See if we have in the 0-1 order and store
            topat1 = bond[0].index
            topat2 = bond[1].index
            order = self._bondorders_by_atomindices[topat1][topat2]
            self._bondorders.append(order)

    def _unrollSMIRKSMatches(self, smirks, aromaticity_model=None):
        """Find all sets of atoms in the topology that match the provided SMIRKS strings.

        Parameters
        ----------
        smirks : str
            SMIRKS string with tagged atoms.
            If there are N tagged atoms numbered 1..N, the resulting matches will be N-tuples of atoms that match the corresponding tagged atoms.
        aromaticity_model : str (optional)
            Default None. Aromaticity model used in SMIRKS matching, as per getSMIRKSMatches_OEMol docs. If provided, pre-processes molecule with this model prior to matching. Otherwise, uses provided oemol.

        Returns
        -------
        matches : list of tuples of Atom
            matches[index] is an N-tuple of Atom entries from the topology
            Matches are returned in no guaranteed order.

        """

        # Perform matching on each unique molecule, unrolling the matches to all matching copies of that molecule in the Topology object.
        matches = list()
        for reference_molecule in self._reference_molecules:
            # Find all atomsets that match this definition in the reference molecule
            refmol_matches = getSMIRKSMatches_OEMol( reference_molecule, smirks, aromaticity_model = aromaticity_model)

            # Loop over matches
            for reference_atom_indices in refmol_matches:
                # Unroll corresponding atom indices over all instances of this molecule
                for reference_to_topology_atom_mapping in self._reference_to_topology_atom_mappings[reference_molecule]:
                    # Create match.
                    atom_indices = tuple([ reference_to_topology_atom_mapping[atom_index] for atom_index in reference_atom_indices ])
                    matches.append(atom_indices)

        return matches

    def _assignPartialCharges(self, molecule, oechargemethod, modifycharges = True):
        """Assign partial charges to the specified molecule using best practices and the OpenEye toolkit.

        Parameters
        ----------
        molecule : OEMol
            The molecule to be charged.
            NOTE: The molecule will be modified when charges are added.
        oechargemethod : str
            The name of the charge method from oequacpac to use (e.g. 'OECharges_AM1BCCSym')
        modifycharges : bool (optional)
            If False, don't actually assign partial charges; use the charge calculation solely to update the Wiberg bond orders.

        .. notes::
            * As per Christopher Bayly and the `Canonical AM1-BCC documentation <http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html>`_, ``OEAssignPartialCharges`` needs multiple conformations to ensure well-behaved charges. This implements that recipe for conformer generation.
            This conformer generation may or may not be necessary if the calculation is only to obtain bond orders; this will have to be investigated separately so it is retained for now.

        """
        # TODO: Cache charged molecules here to save time in future calls to createSystem
        from openeye import oechem, oeomega, oequacpac

        # Expand conformers
        if not oechem.OEChemIsLicensed(): raise(ImportError("Need License for OEChem!"))
        if not oeomega.OEOmegaIsLicensed(): raise(ImportError("Need License for OEOmega!"))
        omega = oeomega.OEOmega()
        omega.SetMaxConfs(800)
        omega.SetCanonOrder(False)
        omega.SetSampleHydrogens(True)
        omega.SetEnergyWindow(15.0)
        omega.SetRMSThreshold(1.0)
        omega.SetStrictStereo(True) #Don't generate random stereoisomer if not specified
        charged_copy = oechem.OEMol(molecule)
        status = omega(charged_copy)
        if not status:
            raise(RuntimeError("Omega returned error code %s" % status))

        # Assign charges
        status = oequacpac.OEAssignPartialCharges(charged_copy, getattr(oequacpac, oechargemethod), False, False)
        if not status:
            raise(RuntimeError("OEAssignPartialCharges returned error code %s" % status))
        # Our copy has the charges we want but not the right conformation. Copy charges over. Also copy over Wiberg bond orders if present
        partial_charges = []
        partial_bondorders = []
        if modifycharges:
            for atom in charged_copy.GetAtoms():
                partial_charges.append( atom.GetPartialCharge() )
            for (idx,atom) in enumerate(molecule.GetAtoms()):
                atom.SetPartialCharge( partial_charges[idx] )
        for bond in charged_copy.GetBonds():
            partial_bondorders.append( bond.GetData("WibergBondOrder"))
        for (idx, bond) in enumerate(molecule.GetBonds()):
            bond.SetData("WibergBondOrder", partial_bondorders[idx])

        # If the charge method was not an OpenEye AM1 method and we need Wiberg bond orders, obtain Wiberg bond orders
        if not (type(chargeMethod) == str and 'AM1' in chargeMethod) and self._use_fractional_bondorder:
            if verbose: print("Doing an AM1 calculation to get Wiberg bond orders.")
            for molecule in molecules:
                # Do AM1 calculation just to get bond orders on moleules (discarding charges)
                self._assignPartialCharges(molecule, "OECharges_AM1", modifycharges = False)


    def _charge_molecules(self, chargeMethod='AM1-BCC'):
        # Charge molecules, if needed
        if chargeMethod == None:
            # Don't charge molecules
            if verbose: print('Charges specified in provided molecules will be used.')
            oechargemethod = None
        elif chargeMethod == 'BCC':
            # Check if we have a BondChargeCorrectionGenerator populated
            force_generators = { force.__class__.__name__ : force for force in self._forces }
            if ('BondChargeCorrectionGenerator' in force_generators):
                oechargemethod = force_generators['BondChargeCorrectionGenerator']._oechargemethod
                if verbose: print('Applying oechem.oequacpac.OEAssignPartialCharges with initial charge method "%s" followed by bond charge corrections.' % oechargemethod)
            else:
                # Don't charge molecules if no bond charge corrections were found
                oechargemethod = None
        elif type(chargeMethod) == str:
            from openeye import oequacpac
            # Check if the user has provided a manually-specified charge method
            if hasattr(oequacpac, chargeMethod):
                oechargemethod = chargeMethod
                if verbose: print('Applying oechem.oequacpac.OEAssignPartialCharges with specified charge method "%s".' % oechargemethod)
            else:
                raise Exception("Unknown chargeMethod '%s'"% chargeMethod)
        else:
            raise Exception("Unknown chargeMethod ''%s'"% str(chargeMethod))

#=============================================================================================
# TOPOLOGY OBJECTS
#=============================================================================================

class Particle(object):
    """
    A particle in a system

    This could be an ``Atom`` or a ``VirtualSite``

    """
    def __init__(self):
        """
        Create a particle.
        """
        pass

    @property
    def particle_index(self):
        """
        Index of this particle within the ``Topology``

        .. todo::

           Should this just be called ``index``, or does that risk confusion within
           the index within ``topology.atoms``?

        """
        pass

    def __repr__(self):
        pass

    def __str__(self):
        pass

class Atom(Particle):
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

    @property
    def atom_index(self):
        """
        The index of this Atom within the the list of atoms in ``Topology``.
        Note that this can be different from ``particle_index``.

        .. todo::

           Do we need this?

        """
        pass

class VirtualSite(Particle):
    """
    A virtual (non-atom) site.

    .. todo::
        * Should virtual sites be attached to one atom only, or more than one atom?
          OpenMM defines them as belonging to two or more atoms.
        * Should a virtual site be able to belong to more than one Topology?
        * Should virtual sites be immutable or mutable?

    """
    def __init__(self, weights, particles):
        """
        Create a virtual site, defined by a linear combination of multiple particles.

        Parameters
        ----------
        weights : list of floats of shape [N]
            weights[index] is the weight of particles[index] contributing to the position of the virtual site.
        particles : list of Particle of shape [N]
            particles[index] is the corresponding particle for weights[index]
            The ``VirtualSite`` is bound to the ``Particle``s in the list specified here.

        """
        self._weights = np.array(weights) # make a copy and convert to array internally
        self._particle = [ particle for particle in particles ] # create a list of Particles

    @property
    def virtual_site_index(self):
        """
        The index of this VirtualSite within the list of virtual sites within ``Topology``
        Note that this can be different from ``particle_index``.

        .. todo::

           Do we need this?

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

    def chemical_environment_matches(self, query):
        """Retrieve all matches for a given chemical environment query.

        TODO:
        * Do we want to generalize this to other kinds of queries too, like mdtraj DSL, pymol selections, atom index slices, etc?
          We could just call it topology.matches(query)

        Parameters
        ----------
        query : str or ChemicalEnvironment
            SMARTS string (with one or more tagged atoms) or ``ChemicalEnvironment`` query
            Query will internally be resolved to SMIRKS using ``query.asSMIRKS()`` if it has an ``.asSMIRKS`` method.

        Returns
        -------
        matches : list of Atom tuples
            A list of all matching Atom tuples

        """
        # Resolve to SMIRKS if needed
        if hasattr(query, 'asSMIRKS'):
            smirks = query.asSMIRKS()
        else:
            smirks = query

        # TODO: Enumerate matches using the currently selected toolkit.
        if self._toolkit == 'oechem':

            oemol = molecule.as_oemol()
            matches = _getSMIRKSMatches_OEMol(oemol, smirks, aromaticity_model=self._aromaticity_model)

        # Perform matching on each unique molecule, unrolling the matches to all matching copies of that molecule in the Topology object.
        matches = list()
        for molecule in self.unique_molecules:
            # Find all atomsets that match this definition in the reference molecule
            refmol_matches = molecule.chemical_environment_matches(smirks)

            # Loop over matches
            for reference_atom_indices in refmol_matches:
                # Unroll corresponding atom indices over all instances of this molecule
                for reference_to_topology_atom_mapping in self._reference_to_topology_atom_mappings[reference_molecule]:
                    # Create match.
                    atom_indices = tuple([ reference_to_topology_atom_mapping[atom_index] for atom_index in reference_atom_indices ])
                    matches.append(atom_indices)

        return matches

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

    def set_aromaticity_model(self, aromaticity_model):
        """
        Set the aromaticity model to use for representing this molecule.

        Parameters
        ----------
        aromaticity_model : str
            Aromaticity model to use for this molecule. One of ['MDL']

        """
        # TODO: Validate aromaticity model against allowed models
        self._aromaticity_model = aromaticity_model

    def get_aromaticity_model(self):
        """
        Retrieve aromaticity model for this molecule.

        Returns
        -------
        aromaticity_model : str
            The aromaticity model in use for this molecule.

        """
        return self._aromaticity_model

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

    def assign_partial_charges(self, method='AM1-BCC', toolkit=None, **kwargs):
        """Assign partial atomic charges.

        .. todo::
            * Is it OK that the ``Molecule`` object does not store geometry, but will create it using ``openeye.omega`` or ``rdkit``?
            * Should this method assign charges to the ``Atom``s in the molecule, a separate ``charges`` molecule property,
              or just return the charge array? Is it OK that the Topology is modified?
            * How do we add enough flexibility to specify the toolkit and optional parameters, such as:
              ``oequacpac.OEAssignPartialCharges(charged_copy, getattr(oequacpac, 'OECharges_AM1BCCSym'), False, False)``
            * What will we do about virtual sites, since this only assigns partial atomic charges?

        Parameters
        ----------
        method : str, optional, default='am1-bcc'
            The name of the charge method to use.
            Options are:
            * `AM1` : symmetrized AM1 charges without BCCs
            * 'AM1-BCC' : symmetrized ELF AM1-BCC charges using best practices
        toolkit : str, optional, default=None
            If specified, the provided toolkit module will be used; otherwise, all toolkits will be tried in undefined order.
            Currently supported options:
            * 'openeye' : generate conformations with ``openeye.omega`` and assign charges with ``openeye.oequacpac``
            * 'rdkit' : generate conformations with ``rdkit`` and assign charges with ``antechamber``
            ``kwargs`` will be passed to the toolkit.
        """
        pass

    def assign_fractional_bond_orders(self, method='Wiberg', toolkit=None, **kwargs):
        """Assign fractional bond orders.

        .. todo::
            * Is it OK that the ``Molecule`` object does not store geometry, but will create it using ``openeye.omega`` or ``rdkit``?
            * Should this method assign fractional bond orders to the ``Bond``s in the molecule, a separate ``bond_orders`` molecule property,
              or just return the array of bond orders?
            * How do we add enough flexibility to specify the toolkit and optional parameters, such as:
              ``oequacpac.OEAssignPartialCharges(charged_copy, getattr(oequacpac, 'OECharges_AM1BCCSym'), False, False)``

        method : str, optional, default='Wiberg'
            The name of the charge method to use.
            Options are:
            * 'Wiberg' : Wiberg bond order
        toolkit : str, optional, default=None
            If specified, the provided toolkit module will be used; otherwise, all toolkits will be tried in undefined order.
            Currently supported options:
            * 'openeye' : generate conformations with ``openeye.omega`` and assign Wiberg bond order with ``openeye.oequacpac`` using OECharges_AM1BCCSym
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
        self._aromaticity_model = DEFAULT_AROMATICITY_MODEL
        self._fractional_bondorder_model = DEFAULT_FRACTIONAL_BONDORDER_MODEL
        self._charge_model = DEFAULT_CHARGE_MODEL
        self._constrained_atom_pairs = dict()
        pass

    # QUESTION: Should aromaticity model instead be specified only when getting SMARTS matches?
    # QUESTION: Should we allow representations with mutliple aromaticity models to be cached?
    def set_aromaticity_model(self, aromaticity_model):
        """
        Set the aromaticity model applied to all molecules in the topology.

        Parameters
        ----------
        aromaticity_model : str
            Aromaticity model to use. One of: ['MDL']

        """
        self._aromaticity_model = aromaticity_model
        for molecule in self.molecules:
            molecule.set_aromaticity_model(aromaticity_model)

    # QUESTION: Should fractional bond order model instead be specified only when retrieving fractional bond orders?
    # QUESTION: Should we allow fractional bond orders with multiple bondorder models to be cached?
    def set_fractional_bondorder_model(self, fractional_bondorder_model):
        """
        Set the fractional bond order model applied to all molecules in the topology.

        Parameters
        ----------
        fractional_bondorder_model : str
            Fractional bond order model to use. One of: ['Wiberg']

        """
        self._fractional_bondorder_model = fractional_bondorder_model
        for molecule in self.molecules:
            molecule.set_fractional_bondorder_model(fractional_bondorder_model)

    # QUESTION: Should charge model instead be specified only when retrieving partial charges?
    # QUESTION: Should we allow partial charge sets with multiple charge models to be cached?
    def set_charge_model(self, charge_model):
        """
        Set the fractional bond order model applied to all molecules in the topology.

        Parameters
        ----------
        charge_model : str
            Charge model to use. One of: ['AM1', 'AM1-BCC', 'Mulliken']

        """
        self._charge_model = charge_model
        for molecule in self.molecules:
            molecule.set_charge_model(charge_model)

    def chemical_environment_matches(self, query):
        """Retrieve all matches for a given chemical environment query.

        TODO:
        * Do we want to generalize this to other kinds of queries too, like mdtraj DSL, pymol selections, atom index slices, etc?
          We could just call it topology.matches(query)

        Parameters
        ----------
        query : str or ChemicalEnvironment
            SMARTS string (with one or more tagged atoms) or ``ChemicalEnvironment`` query
            Query will internally be resolved to SMIRKS using ``query.asSMIRKS()`` if it has an ``.asSMIRKS`` method.

        Returns
        -------
        matches : list of Atom tuples
            A list of all matching Atom tuples
        """
        # Perform matching on each unique molecule, unrolling the matches to all matching copies of that molecule in the Topology object.
        matches = list()
        for molecule in self.unique_molecules:
            # Find all atomsets that match this definition in the reference molecule
            refmol_matches = molecule.chemical_environment_matches(smirks)

            # Loop over matches
            for reference_match in refmol_matches:
                # Unroll corresponding atom indices over all instances of this molecule
                for reference_to_topology_atom_mapping in self._reference_to_topology_atom_mappings[reference_molecule]:
                    # Create match.
                    match = tuple([ reference_to_topology_atom_mapping[atom] for atom in reference_match ])
                    matches.append(match)

        return matches

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

    def add_particle(self, particle):
        """Add a Particle to the Topology.

        Parameters
        ----------
        particle : Particle
            The Particle to be added.
            The Topology will take ownership of the Particle.

        """
        pass

    def add_molecule(self, molecule):
        """Add a Molecule to the Topology.

        Parameters
        ----------
        molecule : Molecule
            The Molecule to be added.
            The Topology will take ownership of the Molecule.
        """
        molecule.set_aromaticity_model(self._aromaticity_model)
        pass

    def add_constraint(self, iatom, jatom, distance=True):
        """
        Mark a pair of atoms as constrained.

        Constraints between atoms that are not bonded (e.g., rigid waters) are permissible.

        Parameters
        ----------
        iatom, jatom : Atom
            Atoms to mark as constrained
            These atoms may be bonded or not in the Topology
        distance : simtk.unit.Quantity, optional, default=True
            Constraint distance
            ``True`` if distance has yet to be determined
            ``False`` if constraint is to be removed

        """
        # Check that constraint hasn't already been specified.
        if (iatom,jatom) in self._constrained_atom_pairs:
            existing_distance = self._constrained_atom_pairs[(iatom,jatom)]
            if unit.is_quantity(existing_distance) and (distance is True):
                raise Exception('Atoms (%d,%d) already constrained with distance %s but attempting to override with unspecified distance' % (iatom, jatom, existing_distance))
            if (existing_distance is True) and (distance is True):
                raise Exception('Atoms (%d,%d) already constrained with unspecified distance but attempting to override with unspecified distance' % (iatom, jatom))
            if distance is False:
                del self._constrained_atom_pairs[(iatom,jatom)]
                del self._constrained_atom_pairs[(jatom,iatom)]
                return

        self._constrained_atom_pairs[(iatom,jatom)] = distance
        self._constrained_atom_pairs[(jatom,iatom)] = distance

    def is_constrained(self, iatom, jatom):
        """
        Check if a pair of atoms are marked as constrained.

        Parameters
        ----------
        iatom, jatom : int
            Indices of atoms to mark as constrained.

        Returns
        -------
        distance : simtk.unit.Quantity or bool
            True if constrained but constraints have not yet been applied
            Distance if constraint has already been added to System

        """
        if (iatom,jatom) in self._constrained_atom_pairs:
            return self._constrained_atom_pairs[(iatom,jatom)]
        else:
            return False

    def get_fractional_bond_order(self, iatom, jatom):
        """
        Retrieve the fractional bond order for a bond.

        An Exception is raised if it cannot be determined.

        Parameters
        ----------
        iatom, jatom : Atom
            Atoms for which a fractional bond order is to be retrieved.

        Returns
        -------
        order : float
            Fractional bond order between the two specified atoms.

        """
        # TODO: Look up fractional bond order in corresponding list of unique molecules,
        # computing it lazily if needed.

        pass
