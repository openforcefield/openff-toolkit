#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Representation of molecular chemical entities.

.. todo::

   * Make all classes (like Particle, Atom, VirtualSite) hashable and serializable.
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

# TODO: Can we have the `ALLOWED_*_MODELS` list automatically appear in the docstrings below?
# TODO: Should `ALLOWED_*_MODELS` be objects instead of strings?
# TODO: Should these be imported from `openforcefield.cheminformatics.aromaticity_models` and `.bondorder_models`?

DEFAULT_AROMATICITY_MODEL = 'MDL' # TODO: Is there a more specific name and reference for the aromaticity model?
ALLOWED_AROMATICITY_MODELS = ['MDL']

DEFAULT_FRACTIONAL_BONDORDER_MODEL = 'Wiberg' # TODO: Is there a more specific name and reference for the aromatciity model?
ALLOWED_FRACTIONAL_BONDORDER_MODELS = ['Wiberg']

DEFAULT_CHARGE_MODEL = 'AM1-CM2' # TODO: Should this be `AM1-BCC`, or should we encode BCCs explicitly?
ALLOWED_CHARGE_MODELS = ['AM1-CM2', 'AM1-BCC', 'Mulliken'] # TODO: Which models do we want to support?

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

from openforcefield.utils import requires_openeye_licenses

@requires_openeye_licenses('oechem')
def _getSMARTSMatches_OEMol(oemol, smarts, aromaticity_model=None):
    """Find all sets of atoms in the provided oemol that match the provided SMARTS string.

    Parameters
    ----------
    oemol : openeye.oechem.OEMol or similar
        oemol to process with the SMIRKS in order to find matches
    smarts : str
        SMARTS string with any number of sequentially tagged atoms.
        If there are N tagged atoms numbered 1..N, the resulting matches will be N-tuples of atoms that match the corresponding tagged atoms.
    aromaticity_model : str, optional, default=None
        OpenEye aromaticity model designation as a string, such as ``OEAroModel_MDL``.
        If ``None``, molecule is processed exactly as provided; otherwise it is prepared with this aromaticity model prior to querying.

    Returns
    -------
    matches : list of tuples of atoms indices within the ``oemol``
        matches[index] is an N-tuple of atom numbers from the ``oemol``
        Matches are returned in no guaranteed order.
        # TODO: What is returned if no matches are found? An empty list, or None?
        # TODO: Ensure that SMARTS numbers 1, 2, 3... are rendered into order of returnd matches indexed by 0, 1, 2...


    .. notes ::

       * Raises ``LicenseError`` if valid OpenEye tools license is not found, rather than causing program to terminate
       * Raises ``ValueError`` if ``smarts`` query is malformed

    """
    # We have wrapped the function with @requires_openeye_licenses so that if valid license is not found,
    # a LicenseError will be raised instead of the program abruptly terminating.
    from openeye import oechem

    # Make a copy of molecule so we don't influence original (probably safer than deepcopy per C Bayly)
    mol = oechem.OEMol(oemol)

    # Set up query
    qmol = oechem.OEQMol()
    if not oechem.OEParseSmarts(qmol, smarts):
        raise ValueError("Error parsing SMARTS '%s'" % smarts)

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
        oechem.OEClearAromaticFlags(mol)
        oechem.OEAssignAromaticFlags(mol, oearomodel)
        # Avoid running OEPrepareSearch or we lose desired aromaticity, so instead:
        oechem.OEAssignHybridization(mol)

    # Build list of matches
    # TODO: Update the logic here to preserve ordering of template molecule for equivalent atoms
    #       and speed matching for larger molecules.
    unique = False # We require all matches, not just one of each kind
    substructure_search = oechem.OESubSearch(qmol)
    matches = list()
    for match in substructure_search.Match(mol, unique):
        # Compile list of atom indices that match the pattern tags
        atom_indices = dict()
        for matched_atom in match.GetAtoms():
            if matched_atom.pattern.GetMapIdx() != 0:
                atom_indices[matched_atom.pattern.GetMapIdx()-1] = matched_atom.target.GetIdx()
        # Compress into list
        atom_indices = [ atom_indices[index] for index in range(len(atom_indices)) ]
        # Convert to tuple
        matches.append( tuple(atom_indices) )

    return matches

#=============================================================================================
# TOPOLOGY OBJECFTS
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

# TODO: Should this be a mixin?
class ChemicalEntity(object):
    """
    Mixin class for properties shared by chemical entities containing more than one atom.

    A ``ChemicalEntity`` can be queried for SMARTS matches.
    # TODO: Should only molecules be queryable for SMARTS matches?

    """
    def __init__(self, other=None):
        """
        Create a new ChemicalEntity.
        """
        self._particles = list()
        self._bonds = None

    def _invalidate_cached_properties(self):
        """
        Indicate that the chemical entity has been altered.
        """
        # List of all possible cached property names
        CACHED_PROPERTY_NAMES = ['_angles', '_propers', '_impropers', '_charges']
        # Delete any cached proprties
        for property_name in CACHED_PROPERTY_NAMES:
            if hasattr(self, property_name):
                delattr(self, property_name)

    def add_atom(self, atom):
        """
        Add an Atom.

        Parameters
        ----------
        atom : Atom
            The Atom to add.
        """
        self._particles.append(atom)
        self._invalidate_cached_properties()

    def add_virtual_site(self, virtual_site):
        """
        Add a Virtual Site.

        Parameters
        ----------
        virtual_site : VirtualSite
            The VirtualSite to add.
        """
        # Make sure that all Atoms referenced in the virtual site are already present in the entity.
        for atom in virtual_site.atoms:
            if atom not in self._particles:
                raise Exception("{} depends on {}, which is not present in the chemical entity".format(virtual_site, atom))
        self._particles.append(virtual_site)

    def n_particles(self):
        """
        The number of Particle objects, which corresponds to how many positions must be used.
        """
        return sum([1 for particle in self.particles])

    @property
    def n_atoms(self):
        """
        The number of Atom objects.
        """
        return sum([1 for atom in self.atoms])

    @property
    def n_virtual_sites(self):
        """
        The number of VirtualSite objects.
        """
        return sum([1 for virtual_site in self.virtual_sites])

    @property
    def n_bonds(self):
        """
        The number of Bond objects.
        """
        return sum([1 for bond in self.bonds])

    @property
    def particles(self):
        """
        Iterate over all Particle objects.
        """
        for particle in self._particles:
            yield particle

    @property
    def atoms(self):
        """
        Iterate over all Atom objects.

        """
        for particle in self._particles:
            if isinstance(particle, Atom):
                yield particle

    @property
    def virtual_sites(self):
        """
        Iterate over all VirtualSite objects.
        """
        for particle in self._particles:
            if isinstance(particle, VirtualSite):
                yield particle

    @property
    def bonds(self):
        """
        Iterate over all Bond objects.

        """
        for bond in self._bonds:
            yield bond

    @property
    def angles(self):
        """
        Iterate over all angles (Atom tuples) in the molecule

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
            This might be a :class:`Molecule` object, a file that can be used to construct a :class:`Molecule` object,
            a serialized :class:`Molecule` object, an ``openeye.oechem.OEMol``, or an ``rdkit.Chem.rdchem.Mol``.
        """
        # Initialize base class
        super(self, Molecule).__init__(other=other)

        if other is not None:
            # TODO: Can we check interface compliance (in a try..except) instead of checking instances?
            if isinstance(other, openforcefield.topology.Molecule):
                self._copy_initializer(other)
            elif issubclass(other, openeye.oechem.OEMolBase):
                mol = Molecule.from_openeye(other)
                self._copy_initializer(mol)
            elif isinstance(other, rdkit.Chem.rdchem.Mol):
                mol = Molecule.from_rdkit(other)
                self._copy_initializer(mol)
            elif isinstance(other, str):
                self.__setstate__(other)
            else:
                msg = 'Cannot construct openforcefield.topology.Molecule from {}\n'.format(other)
                msg += 'other must be '
                raise Exception(msg)

    def __setstate__(self, state):
        # TODO: Implement deserialization
        pass

    def __getstate__(self):
        # TODO: Implement serialization
        pass

    @staticmethod
    def from_file(filename):
        """
        Create one or more molecules from a file

        Parameters
        ----------
        filename : str
            The name of the file to stream one or more molecules from.

        Returns
        -------
        molecules : Molecule or list of Molecules
            If there is a single molecule in the file, a Molecule is returned;
            otherwise, a list of Molecule objects is returned.
            
        """
        raise NotImplementedError()

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
        # TODO: Implement this.
        raise NotImplementedError("RDKit functionality not yet implemented")

    def to_rdkit(self, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an RDKit molecule

        Requires the RDKit to be installed.

        Returns
        -------
        rdmol : rkit.RDMol
            An RDKit molecule
        """
        # TODO: Implement this.
        raise NotImplementedError("RDKit functionality not yet implemented")

    # TODO: Should this method be called from_oemol instead?
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
        # OE Hierarchical molecule view
        hv = oechem.OEHierView(mol, oechem.OEAssumption_BondedResidue +
                               oechem.OEAssumption_ResPerceived +
                               oechem.OEAssumption_PDBOrder)

        # Create empty OpenMM Topology
        topology = app.Topology()
        # Dictionary used to map oe atoms to openmm atoms
        oe_atom_to_openmm_at = {}

        for chain in hv.GetChains():
            # TODO: Fail if hv contains more than one molecule.

            # Create empty OpenMM Chain
            openmm_chain = topology.addChain(chain.GetChainID())

            for frag in chain.GetFragments():

                for hres in frag.GetResidues():

                    # Get OE residue
                    oe_res = hres.GetOEResidue()
                    # Create OpenMM residue
                    openmm_res = topology.addResidue(oe_res.GetName(), openmm_chain)

                    for oe_at in hres.GetAtoms():
                        # Select atom element based on the atomic number
                        element = app.element.Element.getByAtomicNumber(oe_at.GetAtomicNum())
                        # Add atom OpenMM atom to the topology
                        openmm_at = topology.addAtom(oe_at.GetName(), element, openmm_res)
                        openmm_at.index = oe_at.GetIdx()
                        # Add atom to the mapping dictionary
                        oe_atom_to_openmm_at[oe_at] = openmm_at

        if topology.getNumAtoms() != mol.NumAtoms():
            oechem.OEThrow.Error("OpenMM topology and OEMol number of atoms mismatching: "
                                 "OpenMM = {} vs OEMol  = {}".format(topology.getNumAtoms(), mol.NumAtoms()))

        # Count the number of bonds in the openmm topology
        omm_bond_count = 0

        def IsAmideBond(oe_bond):
            # TODO: Can this be replaced by a SMARTS query?

            # This supporting function checks if the passed bond is an amide bond or not.
            # Our definition of amide bond C-N between a Carbon and a Nitrogen atom is:
            #          O
            #          ║
            #  CA or O-C-N-
            #            |

            # The amide bond C-N is a single bond
            if oe_bond.GetOrder() != 1:
                return False

            atomB = oe_bond.GetBgn()
            atomE = oe_bond.GetEnd()

            # The amide bond is made by Carbon and Nitrogen atoms
            if not (atomB.IsCarbon() and atomE.IsNitrogen() or
                    (atomB.IsNitrogen() and atomE.IsCarbon())):
                return False

            # Select Carbon and Nitrogen atoms
            if atomB.IsCarbon():
                C_atom = atomB
                N_atom = atomE
            else:
                C_atom = atomE
                N_atom = atomB

            # Carbon and Nitrogen atoms must have 3 neighbour atoms
            if not (C_atom.GetDegree() == 3 and N_atom.GetDegree() == 3):
                return False

            double_bonds = 0
            single_bonds = 0

            for bond in C_atom.GetBonds():
                # The C-O bond can be single or double.
                if (bond.GetBgn() == C_atom and bond.GetEnd().IsOxygen()) or \
                        (bond.GetBgn().IsOxygen() and bond.GetEnd() == C_atom):
                    if bond.GetOrder() == 2:
                        double_bonds += 1
                    if bond.GetOrder() == 1:
                        single_bonds += 1
                # The CA-C bond is single
                if (bond.GetBgn() == C_atom and bond.GetEnd().IsCarbon()) or \
                        (bond.GetBgn().IsCarbon() and bond.GetEnd() == C_atom):
                    if bond.GetOrder() == 1:
                        single_bonds += 1
            # Just one double and one single bonds are connected to C
            # In this case the bond is an amide bond
            if double_bonds == 1 and single_bonds == 1:
                return True
            else:
                return False

        # Creating bonds
        for oe_bond in mol.GetBonds():
            # Set the bond type
            if oe_bond.GetType() is not "":
                if oe_bond.GetType() in ['Single', 'Double', 'Triple', 'Aromatic', 'Amide']:
                    off_bondtype = oe_bond.GetType()
                else:
                    off_bondtype = None
            else:
                if oe_bond.IsAromatic():
                    oe_bond.SetType("Aromatic")
                    off_bondtype = "Aromatic"
                elif oe_bond.GetOrder() == 2:
                    oe_bond.SetType("Double")
                    off_bondtype = "Double"
                elif oe_bond.GetOrder() == 3:
                    oe_bond.SetType("Triple")
                    off_bond_type = "Triple"
                elif IsAmideBond(oe_bond):
                    oe_bond.SetType("Amide")
                    off_bond_type = "Amide"
                elif oe_bond.GetOrder() == 1:
                    oe_bond.SetType("Single")
                    off_bond_type = "Single"
                else:
                    off_bond_type = None

            molecule.add_bond(oe_atom_to_openmm_at[oe_bond.GetBgn()], oe_atom_to_openmm_at[oe_bond.GetEnd()],
                              type=off_bondtype, order=oe_bond.GetOrder())

        if molecule.n_bondsphe != mol.NumBonds():
            oechem.OEThrow.Error("OpenMM topology and OEMol number of bonds mismatching: "
                                 "OpenMM = {} vs OEMol  = {}".format(omm_bond_count, mol.NumBonds()))

        dic = mol.GetCoords()
        positions = [Vec3(v[0], v[1], v[2]) for k, v in dic.items()] * unit.angstrom

        return topology, positions

    # TODO: Should this method be called to_oemol instead?
    def to_openeye(self, positions=None, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an OpenEye molecule

        Requires the OpenEye toolkit to be installed.

        Returns
        -------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule
        positions : simtk.unit.Quantity with shape [nparticles,3], optional, default=None
            Positions to use in constructing OEMol.
            If virtual sites are present in the Topology, these indices will be skipped.

        NOTE: This comes from https://github.com/oess/oeommtools/blob/master/oeommtools/utils.py

        """
        oe_mol = oechem.OEMol()
        molecule_atom_to_oe_atom = {} # Mapping dictionary between Molecule atoms and oe atoms

        # Python set used to identify atoms that are not in protein residues
        keep = set(proteinResidues).union(dnaResidues).union(rnaResidues)

        for chain in topology.chains():
            for res in chain.residues():
                # Create an OEResidue
                oe_res = oechem.OEResidue()
                # Set OEResidue name
                oe_res.SetName(res.name)
                # If the atom is not a protein atom then set its heteroatom
                # flag to True
                if res.name not in keep:
                    oe_res.SetFragmentNumber(chain.index + 1)
                    oe_res.SetHetAtom(True)
                # Set OEResidue Chain ID
                oe_res.SetChainID(chain.id)
                # res_idx = int(res.id) - chain.index * len(chain._residues)
                # Set OEResidue number
                oe_res.SetResidueNumber(int(res.id))

                for openmm_at in res.atoms():
                    # Create an OEAtom  based on the atomic number
                    oe_atom = oe_mol.NewAtom(openmm_at.element._atomic_number)
                    # Set atom name
                    oe_atom.SetName(openmm_at.name)
                    # Set Symbol
                    oe_atom.SetType(openmm_at.element.symbol)
                    # Set Atom index
                    oe_res.SetSerialNumber(openmm_at.index + 1)
                    # Commit the changes
                    oechem.OEAtomSetResidue(oe_atom, oe_res)
                    # Update the dictionary OpenMM to OE
                    openmm_atom_to_oe_atom[openmm_at] = oe_atom

        if self.n_atoms != oe_mol.NumAtoms():
            raise Exception("OEMol has an unexpected number of atoms: "
                            "Molecule has {} atoms, while OEMol has {} atoms".format(topology.n_atom, oe_mol.NumAtoms()))

        # Create bonds
        for off_bond in self.bonds():
            oe_mol.NewBond(oe_atoms[bond.atom1], oe_atoms[bond.atom2], bond.bond_order)
            if off_bond.type:
                if off_bond.type == 'Aromatic':
                    oe_atom0.SetAromatic(True)
                    oe_atom1.SetAromatic(True)
                    oe_bond.SetAromatic(True)
                    oe_bond.SetType("Aromatic")
                elif off_bond.type in ["Single", "Double", "Triple", "Amide"]:
                    oe_bond.SetType(omm_bond.type)
                else:
                    oe_bond.SetType("")

        if self.n_bonds != oe_mol.NumBonds():
            oechem.OEThrow.Erorr("OEMol has an unexpected number of bonds:: "
                                 "Molecule has {} bonds, while OEMol has {} bonds".format(self.n_bond, oe_mol.NumBonds()))

        if positions is not None:
            # Set the OEMol positions
            particle_indices = [ atom.particle_index for atom in self.atoms ] # get particle indices
            pos = positions[particle_indices].value_in_units_of(unit.angstrom)
            pos = list(itertools.chain.from_iterable(pos))
            oe_mol.SetCoords(pos)
            oechem.OESetDimensionFromCoords(oe_mol)

        return oe_mol

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

    # TODO: Compute terms for each unique molecule, then use mapping to molecules to enumerate all terms
    def angles(self):
        """
        Get an iterator over all i-j-k angles.
        """
        # TODO: This assumes molecules are immutable. If they are mutable, we have to delete ``_angles`` when the atom/bond table is modified.
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

    # TODO: Compute terms for each unique molecule, then use mapping to molecules to enumerate all terms
    # TODO: This assumes molecules are immutable. If they are mutable, we have to delete ``_torsions`` when the atom/bond table is modified.
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
