#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================
"""
Class definitions to represent a molecular system and its chemical components

.. todo::

   * Create MoleculeImage, ParticleImage, AtomImage, VirtualSiteImage here. (Or ``MoleculeInstance``?)
   * Create ``MoleculeGraph`` to represent fozen set of atom elements and bonds that can used as a key for compression
   * Add hierarchical way of traversing Topology (chains, residues)
   * Make all classes hashable and serializable.
   * JSON/BSON representations of objects?
   * Use `attrs <http://www.attrs.org/>`_ for object setter boilerplate?

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import itertools

from collections import MutableMapping
from collections import OrderedDict

from simtk import unit
from simtk.openmm import app

from openforcefield.typing.chemistry import ChemicalEnvironment
from openforcefield.utils.toolkits import (
    DEFAULT_AROMATICITY_MODEL, ALLOWED_AROMATICITY_MODELS, ALLOWED_FRACTIONAL_BOND_ORDER_MODELS,
    GLOBAL_TOOLKIT_REGISTRY, ALLOWED_CHARGE_MODELS
)
from openforcefield.utils.serialization import Serializable
from openforcefield.utils import MessageException

#=============================================================================================
# Exceptions
#=============================================================================================


class DuplicateUniqueMoleculeError(MessageException):
    """
    Exception for when the user provides indistinguishable unique molecules when trying to identify atoms from a PDB
    """
    pass

class NotBondedError(MessageException):
    """
    Exception for when a function requires a bond between two atoms, but none is present
    """
    pass


#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

class _TransformedDict(MutableMapping):
    """A dictionary that transform and sort keys.

    The function __keytransform__ can be inherited to apply an arbitrary
    key-altering function before accessing the keys.

    The function __sortfunc__ can be inherited to specify a particular
    order over which to iterate over the dictionary.

    """

    def __init__(self, *args, **kwargs):
        self.store = OrderedDict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(sorted(self.store, key=self.__sortfunc__))

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key

    @staticmethod
    def __sortfunc__(key):
        return key

# TODO: Encapsulate this atom ordering logic directly into Atom/Bond/Angle/Torsion classes?
class ValenceDict(_TransformedDict):
    """Enforce uniqueness in atom indices."""

    def __keytransform__(self, key):
        """Reverse tuple if first element is larger than last element."""
        # Ensure key is a tuple.
        key = tuple(key)
        # Reverse the key if the first element is bigger than the last.
        if key[0] > key[-1]:
            key = tuple(reversed(key))
        return key


class ImproperDict(_TransformedDict):
    """Symmetrize improper torsions."""

    def __keytransform__(self, key):
        """Reorder tuple in numerical order except for element[1] which is the central atom; it retains its position."""
        # Ensure key is a tuple
        key = tuple(key)
        # Retrieve connected atoms
        connectedatoms = [key[0], key[2], key[3]]
        # Sort connected atoms
        connectedatoms.sort()
        # Re-store connected atoms
        key = tuple(
            [connectedatoms[0], key[1], connectedatoms[1], connectedatoms[2]])
        return (key)


#=============================================================================================
# TOPOLOGY OBJECTS
#=============================================================================================

# =============================================================================================
# TopologyAtom
# =============================================================================================


class TopologyAtom(Serializable):
    """
    A TopologyAtom is a lightweight data structure that represents a single openforcefield.topology.molecule.Atom in
    a Topology. A TopologyAtom consists of two references -- One to its fully detailed "atom", an
    openforcefield.topology.molecule.Atom, and another to its parent "topology_molecule", which occupies a spot in
    the parent Topology's TopologyMolecule list.

    As some systems can be very large, there is no always-existing representation of a TopologyAtom. They are created on
    demand as the user requests them.

    .. warning :: This API is experimental and subject to change.

    """

    def __init__(self, atom, topology_molecule):
        """
        Create a new TopologyAtom.

        Parameters
        ----------
        atom : An openforcefield.topology.molecule.Atom
            The reference atom
        topology_molecule : An openforcefield.topology.TopologyMolecule
            The TopologyMolecule that this TopologyAtom belongs to
        """
        # TODO: Type checks
        self._atom = atom
        self._topology_molecule = topology_molecule

    @property
    def atom(self):
        """
        Get the reference Atom for this TopologyAtom.

        Returns
        -------
        an openforcefield.topology.molecule.Atom
        """
        return self._atom

    @property
    def atomic_number(self):
        """
        Get the atomic number of this atom

        Returns
        -------
        int
        """
        return self._atom.atomic_number

    @property
    def topology_molecule(self):
        """
        Get the TopologyMolecule that this TopologyAtom belongs to.

        Returns
        -------
        openforcefield.topology.TopologyMolecule
        """
        return self._topology_molecule

    @property
    def molecule(self):
        """
        Get the reference Molecule that this TopologyAtom belongs to.

        Returns
        -------
        openforcefield.topology.molecule.Molecule
        """
        return self._topology_molecule.molecule

    @property
    def topology_atom_index(self):
        """
        Get the index of this atom in its parent Topology.

        Returns
        -------
        int
            The index of this atom in its parent topology.
        """
        mapped_molecule_atom_index = self._topology_molecule._ref_to_top_index[
            self._atom.molecule_atom_index]
        return self._topology_molecule.atom_start_topology_index + mapped_molecule_atom_index

    @property
    def topology_particle_index(self):
        """
        Get the index of this particle in its parent Topology.

        Returns
        -------
        int
            The index of this atom in its parent topology.
        """
        # This assumes that particles in a molecule are ordered with all the Atoms first and VirtualSites last.
        mapped_molecule_particle_index = self._topology_molecule._ref_to_top_index[
            self._atom.molecule_atom_index]
        return self._topology_molecule.particle_start_topology_index + mapped_molecule_particle_index

    @property
    def topology_bonds(self):
        """
        Get the TopologyBonds connected to this TopologyAtom.

        Returns
        -------
        iterator of openforcefield.topology.TopologyBonds
        """

        for bond in self._atom.bonds:
            reference_mol_bond_index = bond.molecule_bond_index
            yield self._topology_molecule.bond(reference_mol_bond_index)

    def __eq__(self, other):
        return ((self._atom == other._atom)
                and (self._topology_molecule == other._topology_molecule))

    def __repr__(self):
        return "TopologyAtom {} with reference atom {} and parent TopologyMolecule {}".format(
            self.topology_atom_index, self._atom, self._topology_molecule)

    def to_dict(self):
        """Convert to dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO

    @classmethod
    def from_dict(cls, d):
        """Static constructor from dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO

    #@property
    #def bonds(self):
    #    """
    #    Get the Bonds connected to this TopologyAtom.
    #
    #    Returns
    #    -------
    #    iterator of openforcefield.topology.molecule.Bonds
    #    """
    #    for bond in self._atom.bonds:
    #        yield bond

    # TODO: Add all atom properties here? Or just make people use TopologyAtom.atom for that?


#=============================================================================================
# TopologyBond
#=============================================================================================


class TopologyBond(Serializable):
    """
    A TopologyBond is a lightweight data structure that represents a single openforcefield.topology.molecule.Bond in
    a Topology. A TopologyBond consists of two references -- One to its fully detailed "bond", an
    openforcefield.topology.molecule.Bond, and another to its parent "topology_molecule", which occupies a spot in
    the parent Topology's TopologyMolecule list.

    As some systems can be very large, there is no always-existing representation of a TopologyBond. They are created on
    demand as the user requests them.

    .. warning :: This API is experimental and subject to change.

    """

    def __init__(self, bond, topology_molecule):
        """

        Parameters
        ----------
        bond : An openforcefield.topology.molecule.Bond
            The reference bond.
        topology_molecule : An openforcefield.topology.TopologyMolecule
            The TopologyMolecule that this TopologyBond belongs to.
        """
        # TODO: Type checks
        self._bond = bond
        self._topology_molecule = topology_molecule

    @property
    def bond(self):
        """
        Get the reference Bond for this TopologyBond.

        Returns
        -------
        an openforcefield.topology.molecule.Bond
        """
        return self._bond

    @property
    def topology_molecule(self):
        """
        Get the TopologyMolecule that this TopologyBond belongs to.

        Returns
        -------
        openforcefield.topology.TopologyMolecule
        """
        return self._topology_molecule

    @property
    def topology_bond_index(self):
        """
        Get the index of this bond in its parent Topology.

        Returns
        -------
        int
            The index of this bond in its parent topology.
        """
        return self._topology_molecule.bond_start_topology_index + self._bond.molecule_bond_index

    @property
    def molecule(self):
        """
        Get the reference Molecule that this TopologyBond belongs to.

        Returns
        -------
        openforcefield.topology.molecule.Molecule
        """
        return self._topology_molecule.molecule

    @property
    def bond_order(self):
        """
        Get the order of this TopologyBond.

        Returns
        -------
        int : bond order
        """
        return self._bond.bond_order

    @property
    def atoms(self):
        """
        Get the TopologyAtoms connected to this TopologyBond.

        Returns
        -------
        iterator of openforcefield.topology.TopologyAtom
        """
        for ref_atom in self._bond.atoms:
            yield TopologyAtom(ref_atom, self._topology_molecule)

    def to_dict(self):
        """Convert to dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO

    @classmethod
    def from_dict(cls, d):
        """Static constructor from dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO


#=============================================================================================
# TopologyVirtualSite
#=============================================================================================


class TopologyVirtualSite(Serializable):
    """
    A TopologyVirtualSite is a lightweight data structure that represents a single
    openforcefield.topology.molecule.VirtualSite in a Topology. A TopologyVirtualSite consists of two references --
    One to its fully detailed "VirtualSite", an openforcefield.topology.molecule.VirtualSite, and another to its parent
    "topology_molecule", which occupies a spot in the parent Topology's TopologyMolecule list.

    As some systems can be very large, there is no always-existing representation of a TopologyVirtualSite. They are
    created on demand as the user requests them.

    .. warning :: This API is experimental and subject to change.

    """

    def __init__(self, virtual_site, topology_molecule):
        """

        Parameters
        ----------
        virtual_site : An openforcefield.topology.molecule.VirtualSite
            The reference virtual site
        topology_molecule : An openforcefield.topology.TopologyMolecule
            The TopologyMolecule that this TopologyVirtualSite belongs to
        """
        # TODO: Type checks
        self._virtual_site = virtual_site
        self._topology_molecule = topology_molecule

    def atom(self, index):
        """
        Get the atom at a specific index in this TopologyVirtualSite

        Parameters
        ----------
        index : int
            The index of the atom in the reference VirtualSite to retrieve

        Returns
        -------
        TopologyAtom

        """
        return TopologyAtom(self._virtual_site.atoms[index],
                            self.topology_molecule)
        #for atom in self._virtual_site.atoms:
        #    reference_mol_atom_index = atom.molecule_atom_index
        #    yield self._topology_molecule.atom(reference_mol_atom_index)

    @property
    def atoms(self):
        """
        Get the TopologyAtoms involved in this TopologyVirtualSite.

        Returns
        -------
        iterator of openforcefield.topology.TopologyAtom
        """
        for atom in self._virtual_site.atoms:
            reference_mol_atom_index = atom.molecule_atom_index
            yield self._topology_molecule.atom(reference_mol_atom_index)

    @property
    def virtual_site(self):
        """
        Get the reference VirtualSite for this TopologyVirtualSite.

        Returns
        -------
        an openforcefield.topology.molecule.VirtualSite
        """
        return self._virtual_site

    @property
    def topology_molecule(self):
        """
        Get the TopologyMolecule that this TopologyVirtualSite belongs to.

        Returns
        -------
        openforcefield.topology.TopologyMolecule
        """
        return self._topology_molecule

    @property
    def topology_virtual_site_index(self):
        """
        Get the index of this virtual site in its parent Topology.

        Returns
        -------
        int
            The index of this virtual site in its parent topology.
        """
        return self._topology_molecule.virtual_site_start_topology_index + self._virtual_site.molecule_virtual_site_index

    @property
    def topology_particle_index(self):
        """
        Get the index of this particle in its parent Topology.

        Returns
        -------
        int
            The index of this particle in its parent topology.
        """
        return self._topology_molecule.particle_start_topology_index + self._virtual_site.molecule_particle_index

    @property
    def molecule(self):
        """
        Get the reference Molecule that this TopologyVirtualSite belongs to.

        Returns
        -------
        openforcefield.topology.molecule.Molecule
        """
        return self._topology_molecule.molecule

    @property
    def type(self):
        """
        Get the type of this virtual site

        Returns
        -------
        str : The class name of this virtual site
        """
        return self._virtual_site.type

    def __eq__(self, other):
        return ((self._virtual_site == other._virtual_site)
                and (self._topology_molecule == other._topology_molecule))

    def to_dict(self):
        """Convert to dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO

    @classmethod
    def from_dict(cls, d):
        """Static constructor from dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO


# =============================================================================================
# TopologyMolecule
# =============================================================================================


class TopologyMolecule:
    """
    TopologyMolecules are built to be an efficient way to store large numbers of copies of the same molecule for
    parameterization and system preparation.

    .. warning :: This API is experimental and subject to change.

    """

    def __init__(self,
                 reference_molecule,
                 topology,
                 local_topology_to_reference_index=None):
        """
        Create a new TopologyMolecule.

        Parameters
        ----------
        reference_molecule : an openforcefield.topology.molecule.Molecule
            The reference molecule, with details like formal charges, partial charges, bond orders, partial bond orders,
            and atomic symbols.
        topology : an openforcefield.topology.Topology
            The topology that this TopologyMolecule belongs to
        local_topology_to_reference_index : dict, optional, default=None
            Dictionary of {TopologyMolecule_atom_index : Molecule_atom_index} for the TopologyMolecule that will be built
        """
        # TODO: Type checks
        self._reference_molecule = reference_molecule
        self._topology = topology
        if local_topology_to_reference_index is None:
            local_topology_to_reference_index = dict(
                [(i, i) for i in range(reference_molecule.n_atoms)])
        self._top_to_ref_index = local_topology_to_reference_index
        self._ref_to_top_index = dict(
            (k, j) for j, k in local_topology_to_reference_index.items())

    @property
    def topology(self):
        """
        Get the topology that this TopologyMolecule belongs to

        Returns
        -------
        an openforcefield.topology.Topology
        """
        return self._topology

    @property
    def reference_molecule(self):
        """
        Get the reference molecule for this TopologyMolecule

        Returns
        -------
        an openforcefield.topology.molecule.Molecule
        """
        return self._reference_molecule

    @property
    def n_atoms(self):
        """
        The number of atoms in this topology.

        Returns
        -------
        int
        """
        return self._reference_molecule.n_atoms

    def atom(self, index):
        """
        Get the TopologyAtom with a given topology molecule index in this TopologyMolecule.

        Parameters
        ----------
        index : int
            Index of the TopologyAtom within this TopologyMolecule to retrieve

        Returns
        -------
        an openforcefield.topology.TopologyAtom
        """
        ref_mol_atom_index = self._top_to_ref_index[index]
        return TopologyAtom(self._reference_molecule.atoms[ref_mol_atom_index],
                            self)

    @property
    def atoms(self):
        """
        Return an iterator of all the TopologyAtoms in this TopologyMolecule

        Returns
        -------
        an iterator of openforcefield.topology.TopologyAtoms
        """
        iterate_order = list(self._top_to_ref_index.items())
        # Sort by topology index
        iterate_order.sort(key=lambda x: x[0])
        for top_index, ref_index in iterate_order:
            #self._reference_molecule.atoms:
            yield TopologyAtom(self._reference_molecule.atoms[ref_index], self)

    @property
    def atom_start_topology_index(self):
        """
        Get the topology index of the first atom in this TopologyMolecule

        """
        atom_start_topology_index = 0
        for topology_molecule in self._topology.topology_molecules:
            if self == topology_molecule:
                return atom_start_topology_index
            atom_start_topology_index += topology_molecule.n_atoms

    def bond(self, index):
        """
        Get the TopologyBond with a given reference molecule index in this TopologyMolecule

        Parameters
        ----------
        index : int
            Index of the TopologyBond within this TopologyMolecule to retrieve

        Returns
        -------
        an openforcefield.topology.TopologyBond
        """
        return TopologyBond(self.reference_molecule.bonds[index], self)

    @property
    def bonds(self):
        """
        Return an iterator of all the TopologyBonds in this TopologyMolecule

        Returns
        -------
        an iterator of openforcefield.topology.TopologyBonds
        """
        for bond in self._reference_molecule.bonds:
            yield TopologyBond(bond, self)

    @property
    def n_bonds(self):
        """Get the number of bonds in this TopologyMolecule

        Returns
        -------
        int : number of bonds
        """
        return self._reference_molecule.n_bonds

    @property
    def bond_start_topology_index(self):
        """Get the topology index of the first bond in this TopologyMolecule

        """
        bond_start_topology_index = 0
        for topology_molecule in self._topology.topology_molecules:
            if self == topology_molecule:
                return bond_start_topology_index
            bond_start_topology_index += topology_molecule.n_bonds

    def particle(self, index):
        """
        Get the TopologyParticle with a given reference molecule index in this TopologyMolecule

        Parameters
        ----------
        index : int
            Index of the TopologyParticle within this TopologyMolecule to retrieve

        Returns
        -------
        an openforcefield.topology.TopologyParticle
        """
        return TopologyParticle(self.reference_molecule.particles[index], self)

    @property
    def particles(self):
        """
        Return an iterator of all the TopologyParticles in this TopologyMolecules

        Returns
        -------
        an iterator of openforcefield.topology.TopologyParticle
        """
        # Note: This assumes that particles are all Atoms (in topology map order), and then virtualsites
        yield_order = list(self._top_to_ref_index.items())
        # Sort by topology atom index
        yield_order.sort(key=lambda x: x[0])
        for top_atom_index, ref_mol_atom_index in yield_order:
            ref_atom = self._reference_molecule.atoms[ref_mol_atom_index]
            yield TopologyAtom(ref_atom, self)

        # TODO: Add ordering scheme here
        for vsite in self.reference_molecule.virtual_sites:
            yield TopologyVirtualSite(vsite, self)

        #for particle in self._reference_molecule.particles:
        #    if isinstance(particle, Atom):
        #        yield TopologyAtom(particle, self)
        #    elif isinstance(particle, VirtualSite):
        #        yield TopologyVirtualSite(particle, self)

    @property
    def n_particles(self):
        """Get the number of particles in this TopologyMolecule

        Returns
        -------
        int : The number of particles
        """
        return self._reference_molecule.n_particles

    @property
    def particle_start_topology_index(self):
        """Get the topology index of the first particle in this TopologyMolecule.
        """
        particle_start_topology_index = 0
        for topology_molecule in self._topology.topology_molecules:
            if self == topology_molecule:
                return particle_start_topology_index
            particle_start_topology_index += topology_molecule.n_particles

    def virtual_site(self, index):
        """
        Get the TopologyVirtualSite with a given reference molecule index in this TopologyMolecule

        Parameters
        ----------
        index : int
            Index of the TopologyVirtualSite within this TopologyMolecule to retrieve

        Returns
        -------
        an openforcefield.topology.TopologyVirtualSite
        """
        return TopologyVirtualSite(
            self.reference_molecule.virtual_sites[index], self)

    @property
    def virtual_sites(self):
        """
        Return an iterator of all the TopologyVirtualSites in this TopologyMolecules

        Returns
        -------
        an iterator of openforcefield.topology.TopologyVirtualSite
        """
        for vs in self._reference_molecule.virtual_sites:
            yield TopologyVirtualSite(vs, self)

    @property
    def n_virtual_sites(self):
        """Get the number of virtual sites in this TopologyMolecule

        Returns
        -------
        int
        """
        return self._reference_molecule.n_virtual_sites

    @property
    def angles(self):
        """Iterable of Tuple[TopologyAtom]: iterator over the angles in this Topology."""
        return self._convert_to_topology_atom_tuples(self._reference_molecule.angles)

    @property
    def n_angles(self):
        """int: number of angles in this Topology."""
        return self._reference_molecule.n_angles

    @property
    def propers(self):
        """Iterable of Tuple[TopologyAtom]: iterator over the proper torsions in this Topology."""
        return self._convert_to_topology_atom_tuples(self._reference_molecule.propers)

    @property
    def n_propers(self):
        """int: number of proper torsions in this Topology."""
        return self._reference_molecule.n_propers

    @property
    def impropers(self):
        """Iterable of Tuple[TopologyAtom]: iterator over the improper torsions in this Topology."""
        return self._convert_to_topology_atom_tuples(self._reference_molecule.impropers)

    @property
    def n_impropers(self):
        """int: number of proper torsions in this Topology."""
        return self._reference_molecule.n_impropers

    @property
    def virtual_site_start_topology_index(self):
        """Get the topology index of the first virtual site in this TopologyMolecule
        """
        virtual_site_start_topology_index = 0
        for topology_molecule in self._topology.topology_molecules:
            if self == topology_molecule:
                return virtual_site_start_topology_index
            virtual_site_start_topology_index += topology_molecule.n_virtual_sites

    def to_dict(self):
        """Convert to dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO

    @classmethod
    def from_dict(cls, d):
        """Static constructor from dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO

    def _convert_to_topology_atom_tuples(self, molecule_atom_tuples):
        for atom_tuple in molecule_atom_tuples:
            mol_atom_indices = (a.molecule_atom_index for a in atom_tuple)
            top_mol_atom_indices = (self._ref_to_top_index[mol_idx] for mol_idx in mol_atom_indices)
            yield tuple(self.atom(i) for i in top_mol_atom_indices)


# TODO: pick back up figuring out how we want TopologyMolecules to know their starting TopologyParticle indices

# =============================================================================================
# Topology
# =============================================================================================

# TODO: Revise documentation and remove chains


class Topology(Serializable):
    """
    A Topology is a chemical representation of a system containing one or more molecules appearing in a specified order.

    .. warning :: This API is experimental and subject to change.

    Examples
    --------

    Import some utilities

    >>> from simtk.openmm import app
    >>> from openforcefield.tests.utils import get_data_file_path, get_packmol_pdb_file_path
    >>> pdb_filepath = get_packmol_pdb_file_path('cyclohexane_ethanol_0.4_0.6')
    >>> monomer_names = ('cyclohexane', 'ethanol')

    Create a Topology object from a PDB file and sdf files defining the molecular contents

    >>> from openforcefield.topology import Molecule, Topology
    >>> pdbfile = app.PDBFile(pdb_filepath)
    >>> sdf_filepaths = [get_data_file_path(f'systems/monomers/{name}.sdf') for name in monomer_names]
    >>> unique_molecules = [Molecule.from_file(sdf_filepath) for sdf_filepath in sdf_filepaths]
    >>> topology = Topology.from_openmm(pdbfile.topology, unique_molecules=unique_molecules)

    Create a Topology object from a PDB file and IUPAC names of the molecular contents

    >>> pdbfile = app.PDBFile(pdb_filepath)
    >>> unique_molecules = [Molecule.from_iupac(name) for name in monomer_names]
    >>> topology = Topology.from_openmm(pdbfile.topology, unique_molecules=unique_molecules)

    Create an empty Topology object and add a few copies of a single benzene molecule

    >>> topology = Topology()
    >>> molecule = Molecule.from_iupac('benzene')
    >>> molecule_topology_indices = [topology.add_molecule(molecule) for index in range(10)]

    """

    def __init__(self, other=None):
        """
        Create a new Topology.

        Parameters
        ----------
        other : optional, default=None
            If specified, attempt to construct a copy of the Topology from the specified object.
            This might be a Topology object, or a file that can be used to construct a Topology object
            or serialized Topology object.

        """
        from openforcefield.topology.molecule import FrozenMolecule

        # Assign cheminformatics models
        model = DEFAULT_AROMATICITY_MODEL
        self._aromaticity_model = model
        #self._fractional_bond_order_model = DEFAULT_FRACTIONAL_BOND_ORDER_MODEL
        #self._charge_model = DEFAULT_CHARGE_MODEL

        # Initialize storage
        self._initialize()

        # TODO: Try to construct Topology copy from `other` if specified
        if isinstance(other, Topology):
            self.copy_initializer(other)
        elif isinstance(other, FrozenMolecule):
            self.from_molecules([other])
        elif isinstance(other, OrderedDict):
            self._initialize_from_dict(other)

    def _initialize(self):
        """
        Initializes a blank Topology.
        """
        self._aromaticity_model = DEFAULT_AROMATICITY_MODEL
        self._constrained_atom_pairs = dict()
        self._box_vectors = None
        #self._is_periodic = False
        #self._reference_molecule_dicts = set()
        # TODO: Look into weakref and what it does. Having multiple topologies might cause a memory leak.
        self._reference_molecule_to_topology_molecules = OrderedDict()
        self._topology_molecules = list()

    @property
    def reference_molecules(self):
        """
        Get an iterator of reference molecules in this Topology.

        Returns
        -------
        iterable of openforcefield.topology.Molecule
        """
        for ref_mol in self._reference_molecule_to_topology_molecules.keys():
            yield ref_mol

    @classmethod
    def from_molecules(cls, molecules):
        """
        Create a new Topology object containing one copy of each of the specified molecule(s).

        Parameters
        ----------
        molecules : Molecule or iterable of Molecules
            One or more molecules to be added to the Topology

        Returns
        -------
        topology : Topology
            The Topology created from the specified molecule(s)
        """
        # Ensure that we are working with an iterable
        try:
            iter(molecules)
        except TypeError as te:
            # Make iterable object
            molecules = [molecules]

        # Create Topology and populate it with specified molecules
        topology = cls()
        for molecule in molecules:
            topology.add_molecule(molecule)

        return topology

    def assert_bonded(self, atom1, atom2):
        """
        Raise an exception if the specified atoms are not bonded in the topology.

        Parameters
        ----------
        atom1, atom2 : openforcefield.topology.Atom or int
            The atoms or atom topology indices to check to ensure they are bonded

        """
        if (type(atom1) is int) and (type(atom2) is int):
            atom1 = self.atom(atom1)
            atom2 = self.atom(atom2)

        #else:
        if not (self.is_bonded(atom1, atom2)):
            # TODO: Raise more specific exception.
            raise Exception(
                'Atoms {} and {} are not bonded in topology'.format(
                    atom1, atom2))

    @property
    def aromaticity_model(self):
        """
        Get the aromaticity model applied to all molecules in the topology.

        Returns
        -------
        aromaticity_model : str
            Aromaticity model in use.
        """
        return self._aromaticity_model

    @aromaticity_model.setter
    def aromaticity_model(self, aromaticity_model):
        """
        Set the aromaticity model applied to all molecules in the topology.

        Parameters
        ----------
        aromaticity_model : str
            Aromaticity model to use. One of: ['OEAroModel_MDL']
        """

        if not aromaticity_model in ALLOWED_AROMATICITY_MODELS:
            msg = "Aromaticity model must be one of {}; specified '{}'".format(
                ALLOWED_AROMATICITY_MODELS, aromaticity_model)
            raise ValueError(msg)
        self._aromaticity_model = aromaticity_model

    @property
    def box_vectors(self):
        """Return the box vectors of the topology, if specified
        Returns
        -------
        box_vectors : simtk.unit.Quantity wrapped numpy array
            The unit-wrapped box vectors of this topology
        """
        return self._box_vectors

    @box_vectors.setter
    def box_vectors(self, box_vectors):
        """
        Sets the box vectors to be used for this topology.

        Parameters
        ----------
        box_vectors : simtk.unit.Quantity wrapped numpy array
            The unit-wrapped box vectors

        """
        if box_vectors is None:
            self._box_vectors = None
            return
        if not hasattr(box_vectors, 'unit'):
            raise ValueError("Given unitless box vectors")
        if not (unit.angstrom.is_compatible(box_vectors.unit)):
            raise ValueError(
                "Attempting to set box vectors in units that are incompatible with simtk.unit.Angstrom"
            )

        if hasattr(box_vectors, 'shape'):
            assert box_vectors.shape == (3, )
        else:
            assert len(box_vectors) == 3
        self._box_vectors = box_vectors

    @property
    def charge_model(self):
        """
        Get the partial charge model applied to all molecules in the topology.

        Returns
        -------
        charge_model : str
            Charge model used for all molecules in the Topology.

        """
        return self._charge_model

    @charge_model.setter
    def charge_model(self, charge_model):
        """
        Set the partial charge model used for all molecules in the topology.

        Parameters
        ----------
        charge_model : str
            Charge model to use for all molecules in the Topology.
            Allowed values: ['AM1-BCC']
            * ``AM1-BCC``: Canonical AM1-BCC scheme
        """
        if not charge_model in ALLOWED_CHARGE_MODELS:
            raise ValueError(
                "Charge model must be one of {}; specified '{}'".format(
                    ALLOWED_CHARGE_MODELS, charge_model))
        self._charge_model = charge_model

    @property
    def constrained_atom_pairs(self):
        """Returns the constrained atom pairs of the Topology

        Returns
        -------
        constrained_atom_pairs : dict
             dictionary of the form d[(atom1_topology_index, atom2_topology_index)] = distance (float)
        """
        return self._constrained_atom_pairs

    @property
    def fractional_bond_order_model(self):
        """
        Get the fractional bond order model for the Topology.

        Returns
        -------
        fractional_bond_order_model : str
            Fractional bond order model in use.

        """
        return self._fractional_bond_order_model

    @fractional_bond_order_model.setter
    def fractional_bond_order_model(self, fractional_bond_order_model):
        """
        Set the fractional bond order model applied to all molecules in the topology.

        Parameters
        ----------
        fractional_bond_order_model : str
            Fractional bond order model to use. One of: ['Wiberg']

        """
        if not fractional_bond_order_model in ALLOWED_FRACTIONAL_BOND_ORDER_MODELS:
            raise ValueError(
                "Fractional bond order model must be one of {}; specified '{}'"
                .format(ALLOWED_FRACTIONAL_BOND_ORDER_MODELS,
                        fractional_bond_order_model))
        self._fractional_bond_order_model = fractional_bond_order_model

    @property
    def n_reference_molecules(self):
        """
        Returns the number of reference (unique) molecules in in this Topology.

        Returns
        -------
        n_reference_molecules : int
        """
        count = 0
        for i in self.reference_molecules:
            count += 1
        return count

    @property
    def n_topology_molecules(self):
        """
        Returns the number of topology molecules in in this Topology.

        Returns
        -------
        n_topology_molecules : int
        """
        return len(self._topology_molecules)

    @property
    def topology_molecules(self):
        """Returns an iterator over all the TopologyMolecules in this Topology

        Returns
        -------
        topology_molecules : Iterable of TopologyMolecule
        """
        return self._topology_molecules

    @property
    def n_topology_atoms(self):
        """
        Returns the number of topology atoms in in this Topology.

        Returns
        -------
        n_topology_atoms : int
        """
        n_atoms = 0
        for reference_molecule in self.reference_molecules:
            n_atoms_per_topology_molecule = reference_molecule.n_atoms
            n_instances_of_topology_molecule = len(
                self.
                _reference_molecule_to_topology_molecules[reference_molecule])
            n_atoms += n_atoms_per_topology_molecule * n_instances_of_topology_molecule
        return n_atoms

    @property
    def topology_atoms(self):
        """Returns an iterator over the atoms in this Topology. These will be in ascending order of topology index (Note
        that this is not necessarily the same as the reference molecule index)

        Returns
        -------
        topology_atoms : Iterable of TopologyAtom
        """
        for topology_molecule in self._topology_molecules:
            for atom in topology_molecule.atoms:
                yield atom

    @property
    def n_topology_bonds(self):
        """
        Returns the number of TopologyBonds in in this Topology.

        Returns
        -------
        n_bonds : int
        """
        n_bonds = 0
        for reference_molecule in self.reference_molecules:
            n_bonds_per_topology_molecule = reference_molecule.n_bonds
            n_instances_of_topology_molecule = len(
                self.
                _reference_molecule_to_topology_molecules[reference_molecule])
            n_bonds += n_bonds_per_topology_molecule * n_instances_of_topology_molecule
        return n_bonds

    @property
    def topology_bonds(self):
        """Returns an iterator over the bonds in this Topology

        Returns
        -------
        topology_bonds : Iterable of TopologyBond
        """
        for topology_molecule in self._topology_molecules:
            for bond in topology_molecule.bonds:
                yield bond

    @property
    def n_topology_particles(self):
        """
        Returns the number of topology particles (TopologyAtoms and TopologyVirtualSites) in in this Topology.

        Returns
        -------
        n_topology_particles : int
        """
        n_particles = 0
        for reference_molecule in self.reference_molecules:
            n_particles_per_topology_molecule = reference_molecule.n_particles
            n_instances_of_topology_molecule = len(
                self.
                _reference_molecule_to_topology_molecules[reference_molecule])
            n_particles += n_particles_per_topology_molecule * n_instances_of_topology_molecule
        return n_particles

    @property
    def topology_particles(self):
        """Returns an iterator over the particles (TopologyAtoms and TopologyVirtualSites) in this Topology. The
        TopologyAtoms will be in order of ascending Topology index (Note that this may differ from the
        order of atoms in the reference molecule index).

        Returns
        --------
        topology_particles : Iterable of TopologyAtom and TopologyVirtualSite
        """
        for topology_molecule in self._topology_molecules:
            for particle in topology_molecule.particles:
                yield particle

    @property
    def n_topology_virtual_sites(self):
        """
        Returns the number of TopologyVirtualSites in in this Topology.

        Returns
        -------
        n_virtual_sites : iterable of TopologyVirtualSites
        """
        n_virtual_sites = 0
        for reference_molecule in self.reference_molecules:
            n_virtual_sites_per_topology_molecule = reference_molecule.n_virtual_sites
            n_instances_of_topology_molecule = len(
                self.
                _reference_molecule_to_topology_molecules[reference_molecule])
            n_virtual_sites += n_virtual_sites_per_topology_molecule * n_instances_of_topology_molecule
        return n_virtual_sites

    @property
    def topology_virtual_sites(self):
        """Get an iterator over the virtual sites in this Topology

        Returns
        -------
        topology_virtual_sites : Iterable of TopologyVirtualSite
        """
        for topology_molecule in self._topology_molecules:
            for virtual_site in topology_molecule.virtual_sites:
                yield virtual_site

    @property
    def n_angles(self):
        """int: number of angles in this Topology."""
        return sum(mol.n_angles for mol in self._topology_molecules)

    @property
    def angles(self):
        """Iterable of Tuple[TopologyAtom]: iterator over the angles in this Topology."""
        for topology_molecule in self._topology_molecules:
            for angle in topology_molecule.angles:
                yield angle

    @property
    def n_propers(self):
        """int: number of proper torsions in this Topology."""
        return sum(mol.n_propers for mol in self._topology_molecules)

    @property
    def propers(self):
        """Iterable of Tuple[TopologyAtom]: iterator over the proper torsions in this Topology."""
        for topology_molecule in self._topology_molecules:
            for proper in topology_molecule.propers:
                yield proper

    @property
    def n_impropers(self):
        """int: number of improper torsions in this Topology."""
        return sum(mol.n_impropers for mol in self._topology_molecules)

    @property
    def impropers(self):
        """Iterable of Tuple[TopologyAtom]: iterator over the improper torsions in this Topology."""
        for topology_molecule in self._topology_molecules:
            for improper in topology_molecule.impropers:
                yield improper

    def chemical_environment_matches(self,
                                     query,
                                     aromaticity_model='MDL',
                                     toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Retrieve all matches for a given chemical environment query.

        TODO:
        * Do we want to generalize this to other kinds of queries too, like mdtraj DSL, pymol selections, atom index slices, etc?
          We could just call it topology.matches(query)

        Parameters
        ----------
        query : str or ChemicalEnvironment
            SMARTS string (with one or more tagged atoms) or ``ChemicalEnvironment`` query
            Query will internally be resolved to SMARTS using ``query.as_smarts()`` if it has an ``.as_smarts`` method.
        aromaticity_model : str
            Override the default aromaticity model for this topology and use the specified aromaticity model instead.
            Allowed values: ['MDL']

        Returns
        -------
        matches : list of TopologyAtom tuples
            A list of all matching Atom tuples

        """
        # Render the query to a SMARTS string
        if type(query) is str:
            smarts = query
        elif type(query) is ChemicalEnvironment:
            smarts = query.as_smarts()
        else:
            raise ValueError(
                "Don't know how to convert query '%s' into SMARTS string" %
                query)

        # Perform matching on each unique molecule, unrolling the matches to all matching copies of that molecule in the Topology object.
        matches = list()
        for ref_mol in self.reference_molecules:
            # Find all atomsets that match this definition in the reference molecule
            # This will automatically attempt to match chemically identical atoms in a canonical order within the Topology
            refmol_matches = ref_mol.chemical_environment_matches(
                smarts, toolkit_registry=toolkit_registry)

            # Loop over matches
            for reference_match in refmol_matches:
                #mol_dict = molecule.to_dict
                # Unroll corresponding atom indices over all instances of this molecule
                for topology_molecule in self._reference_molecule_to_topology_molecules[
                        ref_mol]:
                    match = list()
                    # Create match TopologyAtoms.
                    for reference_molecule_atom_index in reference_match:
                        reference_atom = topology_molecule._reference_molecule.atoms[
                            reference_molecule_atom_index]
                        topology_atom = TopologyAtom(reference_atom,
                                                     topology_molecule)
                        match.append(topology_atom)
                        #topology_molecule_atom_index = topology_molecule._ref_to_top_index[reference_molecule_atom_index]
                        #atom_topology_index = topology_molecule.atom_start_topology_index+topology_molecule_atom_index
                        #match.append(self.atom(atom_topology_index))
                    match = tuple(match)
                    #match = tuple([topology_molecule.atom_start_topology_index+ref_mol_atom_index for ref_mol_atom_index in reference_match])
                    matches.append(match)

        return matches

    def to_dict(self):
        """Convert to dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO

    @classmethod
    def from_dict(cls, d):
        """Static constructor from dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO

    # TODO: Merge this into Molecule.from_networkx if/when we implement that.
    @staticmethod
    def _networkx_to_hill_formula(mol_graph):
        """
        Convert a networkX representation of a molecule to a molecule formula. Used in printing out
        informative error messages when a molecule from an openmm topology can't be matched.

        Parameters
        ----------
        mol_graph : a networkX graph
            The graph representation of a molecule

        Returns
        -------
        formula : str
            The molecular formula of the graph molecule
        """
        from simtk.openmm.app import Element

        # Make a flat list of all atomic numbers in the molecule
        atom_nums = []
        for idx in mol_graph.nodes:
            atom_nums.append(mol_graph.node[idx]['atomic_number'])

        # Count the number of instances of each atomic number
        at_num_to_counts = dict([(unq, atom_nums.count(unq)) for unq in atom_nums])

        symbol_to_counts = {}
        # Check for C and H first, to make a correct hill formula (remember dicts in python 3.6+ are ordered)
        if 6 in at_num_to_counts:
            symbol_to_counts['C'] = at_num_to_counts[6]
            del at_num_to_counts[6]

        if 1 in at_num_to_counts:
            symbol_to_counts['H'] = at_num_to_counts[1]
            del at_num_to_counts[1]

        # Now count instances of all elements other than C and H, in order of ascending atomic number
        sorted_atom_nums = sorted(at_num_to_counts.keys())
        for atom_num in sorted_atom_nums:
            symbol_to_counts[Element.getByAtomicNumber(atom_num).symbol] = at_num_to_counts[atom_num]

        # Finally format the formula as string
        formula = ''
        for ele, count in symbol_to_counts.items():
            formula += f'{ele}{count}'
        return(formula)


    @classmethod
    def from_openmm(cls, openmm_topology, unique_molecules=None):
        """
        Construct an openforcefield Topology object from an OpenMM Topology object.

        Parameters
        ----------
        openmm_topology : simtk.openmm.app.Topology
            An OpenMM Topology object
        unique_molecules : iterable of objects that can be used to construct unique Molecule objects
            All unique molecules must be provided, in any order, though multiple copies of each molecule are allowed.
            The atomic elements and bond connectivity will be used to match the reference molecules
            to molecule graphs appearing in the OpenMM ``Topology``. If bond orders are present in the
            OpenMM ``Topology``, these will be used in matching as well.
            If all bonds have bond orders assigned in ``mdtraj_topology``, these bond orders will be used to attempt to construct
            the list of unique Molecules if the ``unique_molecules`` argument is omitted.

        Returns
        -------
        topology : openforcefield.topology.Topology
            An openforcefield Topology object
        """
        import networkx as nx
        from networkx.algorithms.isomorphism import GraphMatcher

        # Check to see if the openMM system has defined bond orders, by looping over all Bonds in the Topology.
        omm_has_bond_orders = True
        for omm_bond in openmm_topology.bonds():
            if omm_bond.order is None:
                omm_has_bond_orders = False

        # Set functions for determining equality between nodes and edges
        node_match_func = lambda x, y: x['atomic_number'] == y['atomic_number']
        if omm_has_bond_orders:
            edge_match_func = lambda x, y: x['bond_order'] == y['bond_order']
        else:
            edge_match_func = None

        # Convert all unique mols to graphs
        topology = cls()
        graph_to_unq_mol = {}
        for unq_mol in unique_molecules:
            unq_mol_graph = unq_mol.to_networkx()
            for existing_graph in graph_to_unq_mol.keys():
                if nx.is_isomorphic(
                        existing_graph,
                        unq_mol_graph,
                        node_match=node_match_func,
                        edge_match=edge_match_func):
                    msg = "Error: Two unique molecules have indistinguishable " \
                          "graphs: {} and {}".format(unq_mol, graph_to_unq_mol[existing_graph])
                    raise DuplicateUniqueMoleculeError(msg)
            graph_to_unq_mol[unq_mol_graph] = unq_mol

        # Convert all openMM mols to graphs
        omm_topology_G = nx.Graph()
        for atom in openmm_topology.atoms():
            omm_topology_G.add_node(
                atom.index, atomic_number=atom.element.atomic_number)
        for bond in openmm_topology.bonds():
            omm_topology_G.add_edge(
                bond.atom1.index, bond.atom2.index, bond_order=bond.order)

        # For each connected subgraph (molecule) in the topology, find its match in unique_molecules
        topology_molecules_to_add = list()
        for omm_mol_G in (omm_topology_G.subgraph(c).copy()
                          for c in nx.connected_components(omm_topology_G)):
            match_found = False
            for unq_mol_G in graph_to_unq_mol.keys():
                if nx.is_isomorphic(
                        unq_mol_G,
                        omm_mol_G,
                        node_match=node_match_func,
                        edge_match=edge_match_func):
                    # Take the first valid atom indexing map
                    GM = GraphMatcher(
                        omm_mol_G,
                        unq_mol_G,
                        node_match=node_match_func,
                        edge_match=edge_match_func)
                    for mapping in GM.isomorphisms_iter():
                        topology_atom_map = mapping
                        break
                    first_topology_atom_index = min(topology_atom_map.keys())
                    topology_molecules_to_add.append(
                        (first_topology_atom_index, unq_mol_G,
                         topology_atom_map.items()))
                    match_found = True
                    break
            if not (match_found):
                hill_formula = Topology._networkx_to_hill_formula(omm_mol_G)
                raise ValueError(f'No match found for molecule {hill_formula}')

        # The connected_component_subgraph function above may have scrambled the molecule order, so sort molecules
        # by their first atom's topology index
        topology_molecules_to_add.sort(key=lambda x: x[0])
        for first_index, unq_mol_G, top_to_ref_index in topology_molecules_to_add:
            local_top_to_ref_index = dict(
                [(top_index - first_index, ref_index)
                 for top_index, ref_index in top_to_ref_index])
            topology.add_molecule(
                graph_to_unq_mol[unq_mol_G],
                local_topology_to_reference_index=local_top_to_ref_index)

        topology.box_vectors = openmm_topology.getPeriodicBoxVectors()
        # TODO: How can we preserve metadata from the openMM topology when creating the OFF topology?
        return topology

    def to_openmm(self):
        """
        Create an OpenMM Topology object.

        The OpenMM ``Topology`` object will have one residue per topology
        molecule. Currently, the number of chains depends on how many copies
        of the same molecule are in the ``Topology``. Molecules with more
        than 5 copies are all assigned to a single chain, otherwise one
        chain is created for each molecule. This behavior may change in
        the future.

        Parameters
        ----------
        openmm_topology : simtk.openmm.app.Topology
            An OpenMM Topology object
        """
        from simtk.openmm.app import Topology as OMMTopology
        from simtk.openmm.app import Single, Double, Triple, Aromatic
        from simtk.openmm.app.element import Element as OMMElement

        omm_topology = OMMTopology()

        # Keep track of which chains and residues have been added.
        mol_to_chains = {}
        mol_to_residues = {}

        # Go through atoms in OpenFF to preserve the order.
        omm_atoms = []
        # We need to iterate over the topology molecules if we want to
        # keep track of chains/residues as Atom.topology_molecule is
        # instantiated every time and can't be used as a key.
        for topology_molecule in self.topology_molecules:
            for atom in topology_molecule.atoms:
                reference_molecule = topology_molecule.reference_molecule
                n_molecules = len(self._reference_molecule_to_topology_molecules[reference_molecule])

                # Add 1 chain per molecule unless there are more than 5 copies,
                # in which case we add a single chain for all of them.
                if n_molecules <= 5:
                    # We associate a chain to each molecule.
                    key_molecule = topology_molecule
                else:
                    # We associate a chain to all the topology molecule.
                    key_molecule = reference_molecule

                # Create a new chain if it doesn't exit.
                try:
                    chain = mol_to_chains[key_molecule]
                except KeyError:
                    chain = omm_topology.addChain()
                    mol_to_chains[key_molecule] = chain

                # Add one molecule for each topology molecule.
                try:
                    residue = mol_to_residues[topology_molecule]
                except KeyError:
                    residue = omm_topology.addResidue(reference_molecule.name, chain)
                    mol_to_residues[topology_molecule] = residue

                # Add atom.
                element = OMMElement.getByAtomicNumber(atom.atomic_number)
                omm_atom = omm_topology.addAtom(atom.atom.name, element, residue)

                # Make sure that OpenFF and OpenMM Topology atoms have the same indices.
                assert atom.topology_atom_index == int(omm_atom.id)-1
                omm_atoms.append(omm_atom)

        # Add all bonds.
        bond_types = {
            1: Single,
            2: Double,
            3: Triple
        }
        for bond in self.topology_bonds:
            atom1, atom2 = bond.atoms
            atom1_idx, atom2_idx = atom1.topology_atom_index, atom2.topology_atom_index
            bond_type = Aromatic if bond.bond.is_aromatic else bond_types[bond.bond_order]
            omm_topology.addBond(omm_atoms[atom1_idx], omm_atoms[atom2_idx],
                                 type=bond_type, order=bond.bond_order)
        return omm_topology

    @staticmethod
    def from_mdtraj(mdtraj_topology, unique_molecules=None):
        """
        Construct an openforcefield Topology object from an MDTraj Topology object.

        Parameters
        ----------
        mdtraj_topology : mdtraj.Topology
            An MDTraj Topology object
        unique_molecules : iterable of objects that can be used to construct unique Molecule objects
            All unique molecules mult be provided, in any order, though multiple copies of each molecule are allowed.
            The atomic elements and bond connectivity will be used to match the reference molecules
            to molecule graphs appearing in the MDTraj ``Topology``. If bond orders are present in the
            MDTraj ``Topology``, these will be used in matching as well.
            If all bonds have bond orders assigned in ``mdtraj_topology``, these bond orders will be used to attempt to construct
            the list of unique Molecules if the ``unique_molecules`` argument is omitted.

        Returns
        -------
        topology : openforcefield.Topology
            An openforcefield Topology object
        """
        return Topology.from_openmm(
            mdtraj_topology.to_openmm(), unique_molecules=unique_molecules)

    # TODO: Jeff prepended an underscore on this before 0.2.0 release to remove it from the API.
    #       Before exposing this, we should look carefully at the information that is preserved/lost during this
    #       conversion, and make it clear what would happen to this information in a round trip. For example,
    #       we should know what would happen to formal and partial bond orders and charges, stereochemistry, and
    #       multi-conformer information. It will be important to document these risks to users, as all of these
    #       factors could lead to unintended behavior during system parameterization.
    def _to_mdtraj(self):
        """
        Create an MDTraj Topology object.

        Returns
        ----------
        mdtraj_topology : mdtraj.Topology
            An MDTraj Topology object
        """
        import mdtraj as md
        return md.Topology.from_openmm(self.to_openmm())

    @staticmethod
    def from_parmed(parmed_structure, unique_molecules=None):
        """
        .. warning:: This functionality will be implemented in a future toolkit release.

        Construct an openforcefield Topology object from a ParmEd Structure object.

        Parameters
        ----------
        parmed_structure : parmed.Structure
            A ParmEd structure object
        unique_molecules : iterable of objects that can be used to construct unique Molecule objects
            All unique molecules mult be provided, in any order, though multiple copies of each molecule are allowed.
            The atomic elements and bond connectivity will be used to match the reference molecules
            to molecule graphs appearing in the structure's ``topology`` object. If bond orders are present in the
            structure's ``topology`` object, these will be used in matching as well.
            If all bonds have bond orders assigned in the structure's ``topology`` object,
            these bond orders will be used to attempt to construct
            the list of unique Molecules if the ``unique_molecules`` argument is omitted.

        Returns
        -------
        topology : openforcefield.Topology
            An openforcefield Topology object
        """
        import parmed
        # TODO: Implement functionality
        raise NotImplementedError

    def to_parmed(self):
        """

        .. warning:: This functionality will be implemented in a future toolkit release.

        Create a ParmEd Structure object.

        Returns
        ----------
        parmed_structure : parmed.Structure
            A ParmEd Structure objecft
        """
        import parmed
        # TODO: Implement functionality
        raise NotImplementedError

    # TODO: Jeff prepended an underscore on this before 0.2.0 release to remove it from the API.
    #       This function is deprecated and expects the OpenEye toolkit. We need to discuss what
    #       to do with this functionality in light of our move to the ToolkitWrapper architecture.
    #       Also, as written, this function implies several things about our toolkit's ability to
    #       handle biopolymers. We need to discuss how much of this functionality we will expose
    #       and how we can make our toolkit's current scope clear to users..
    @staticmethod
    def _from_openeye(oemol):
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
        # TODO: Convert this to cls.from_molecules(Molecule.from_openeye())?
        # OE Hierarchical molecule view
        hv = oechem.OEHierView(
            oemol, oechem.OEAssumption_BondedResidue +
            oechem.OEAssumption_ResPerceived + oechem.OEAssumption_PDBOrder)

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
                    openmm_res = topology.addResidue(oe_res.GetName(),
                                                     openmm_chain)

                    for oe_at in hres.GetAtoms():
                        # Select atom element based on the atomic number
                        element = app.element.Element.getByAtomicNumber(
                            oe_at.GetAtomicNum())
                        # Add atom OpenMM atom to the topology
                        openmm_at = topology.addAtom(oe_at.GetName(), element,
                                                     openmm_res)
                        openmm_at.index = oe_at.GetIdx()
                        # Add atom to the mapping dictionary
                        oe_atom_to_openmm_at[oe_at] = openmm_at

        if topology.getNumAtoms() != mol.NumAtoms():
            oechem.OEThrow.Error(
                "OpenMM topology and OEMol number of atoms mismatching: "
                "OpenMM = {} vs OEMol  = {}".format(topology.getNumAtoms(),
                                                    mol.NumAtoms()))

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
                if oe_bond.GetType() in [
                        'Single', 'Double', 'Triple', 'Aromatic', 'Amide'
                ]:
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
                    off_bondtype = "Triple"
                elif IsAmideBond(oe_bond):
                    oe_bond.SetType("Amide")
                    off_bondtype = "Amide"
                elif oe_bond.GetOrder() == 1:
                    oe_bond.SetType("Single")
                    off_bondtype = "Single"
                else:
                    off_bondtype = None

            molecule.add_bond(
                oe_atom_to_openmm_at[oe_bond.GetBgn()],
                oe_atom_to_openmm_at[oe_bond.GetEnd()],
                type=off_bondtype,
                order=oe_bond.GetOrder())

        if molecule.n_bondsphe != mol.NumBonds():
            oechem.OEThrow.Error(
                "OpenMM topology and OEMol number of bonds mismatching: "
                "OpenMM = {} vs OEMol  = {}".format(omm_bond_count,
                                                    mol.NumBonds()))

        dic = mol.GetCoords()
        positions = [Vec3(v[0], v[1], v[2])
                     for k, v in dic.items()] * unit.angstrom

        return topology, positions

    # TODO: Jeff prepended an underscore on this before 0.2.0 release to remove it from the API.
    #       This function is deprecated and expects the OpenEye toolkit. We need to discuss what
    #       to do with this functionality in light of our move to the ToolkitWrapper architecture.
    #       It also expects Topology to be organized by chain, which is not currently the case.
    #       Bringing this function back would require non-trivial engineering and testing, and we
    #       would want to discuss what "guarantee" of correctness it could offer.
    def _to_openeye(self,
                   positions=None,
                   aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an OpenEye OEMol from the topology

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
            raise Exception(
                "OEMol has an unexpected number of atoms: "
                "Molecule has {} atoms, while OEMol has {} atoms".format(
                    topology.n_atom, oe_mol.NumAtoms()))

        # Create bonds
        for off_bond in self.topology_bonds():
            oe_mol.NewBond(oe_atoms[bond.atom1], oe_atoms[bond.atom2],
                           bond.bond_order)
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
            oechem.OEThrow.Erorr(
                "OEMol has an unexpected number of bonds:: "
                "Molecule has {} bonds, while OEMol has {} bonds".format(
                    self.n_bond, oe_mol.NumBonds()))

        if positions is not None:
            # Set the OEMol positions
            particle_indices = [
                atom.particle_index for atom in self.topology_atoms
            ]  # get particle indices
            pos = positions[particle_indices].value_in_units_of(unit.angstrom)
            pos = list(itertools.chain.from_iterable(pos))
            oe_mol.SetCoords(pos)
            oechem.OESetDimensionFromCoords(oe_mol)

        return oe_mol

    def get_bond_between(self, i, j):
        """Returns the bond between two atoms

        Parameters
        ----------
        i, j : int or TopologyAtom
            Atoms or atom indices to check

        Returns
        -------
        bond : TopologyBond
            The bond between i and j.

        """
        if (type(i) is int) and (type(j) is int):
            atomi = self.atom(i)
            atomj = self.atom(j)
        elif (type(i) is TopologyAtom) and (type(j) is TopologyAtom):
            atomi = i
            atomj = j
        else:
            raise Exception(
                "Invalid input passed to is_bonded(). Expected ints or TopologyAtoms, "
                "got {} and {}".format(i, j))

        for top_bond in atomi.topology_bonds:
            for top_atom in top_bond.atoms:
                if top_atom == atomi:
                    continue
                if top_atom == atomj:
                    return top_bond

        raise NotBondedError('No bond between atom {} and {}'.format(i, j))


    def is_bonded(self, i, j):
        """Returns True if the two atoms are bonded

        Parameters
        ----------
        i, j : int or TopologyAtom
            Atoms or atom indices to check

        Returns
        -------
        is_bonded : bool
            True if atoms are bonded, False otherwise.

        """
        try:
            self.get_bond_between(i, j)
            return True
        except NotBondedError:
            return False

    def atom(self, atom_topology_index):
        """
        Get the TopologyAtom at a given Topology atom index.

        Parameters
        ----------
        atom_topology_index : int
             The index of the TopologyAtom in this Topology

        Returns
        -------
        An openforcefield.topology.TopologyAtom
        """
        assert type(atom_topology_index) is int
        assert 0 <= atom_topology_index < self.n_topology_atoms
        this_molecule_start_index = 0
        next_molecule_start_index = 0
        for topology_molecule in self._topology_molecules:
            next_molecule_start_index += topology_molecule.n_atoms
            if next_molecule_start_index > atom_topology_index:
                atom_molecule_index = atom_topology_index - this_molecule_start_index
                # NOTE: the index here should still be in the topology index order, NOT the reference molecule's
                return topology_molecule.atom(atom_molecule_index)
            this_molecule_start_index += topology_molecule.n_atoms

        # Potentially more computationally efficient lookup ( O(largest_molecule_natoms)? )
        # start_index_2_top_mol is an ordered dict of [starting_atom_index] --> [topology_molecule]
        # search_range = range(atom_topology_index - largest_molecule_natoms, atom_topology_index)
        # search_index = atom_topology_index
        # while not(search_index in start_index_2_top_mol.keys()): # Only efficient if start_index_2_top_mol.keys() is a set (constant time lookups)
        #     search_index -= 1
        # topology_molecule = start_index_2_top_mol(search_index)
        # atom_molecule_index = atom_topology_index - search_index
        # return topology_molecule.atom(atom_molecule_index)

    def virtual_site(self, vsite_topology_index):
        """
        Get the TopologyAtom at a given Topology atom index.

        Parameters
        ----------
        vsite_topology_index : int
             The index of the TopologyVirtualSite in this Topology

        Returns
        -------
        An openforcefield.topology.TopologyVirtualSite

        """
        assert type(vsite_topology_index) is int
        assert 0 <= vsite_topology_index < self.n_topology_virtual_sites
        this_molecule_start_index = 0
        next_molecule_start_index = 0
        for topology_molecule in self._topology_molecules:
            next_molecule_start_index += topology_molecule.n_virtual_sites
            if next_molecule_start_index > vsite_topology_index:
                vsite_molecule_index = vsite_topology_index - this_molecule_start_index
                return topology_molecule.virtual_site(vsite_molecule_index)
            this_molecule_start_index += topology_molecule.n_virtual_sites

    def bond(self, bond_topology_index):
        """
        Get the TopologyBond at a given Topology bond index.

        Parameters
        ----------
        bond_topology_index : int
             The index of the TopologyBond in this Topology

        Returns
        -------
        An openforcefield.topology.TopologyBond
        """
        assert type(bond_topology_index) is int
        assert 0 <= bond_topology_index < self.n_topology_bonds
        this_molecule_start_index = 0
        next_molecule_start_index = 0
        for topology_molecule in self._topology_molecules:
            next_molecule_start_index += topology_molecule.n_bonds
            if next_molecule_start_index > bond_topology_index:
                bond_molecule_index = bond_topology_index - this_molecule_start_index
                return topology_molecule.bond(bond_molecule_index)
            this_molecule_start_index += topology_molecule.n_bonds

    def add_particle(self, particle):
        """Add a Particle to the Topology.

        Parameters
        ----------
        particle : Particle
            The Particle to be added.
            The Topology will take ownership of the Particle.

        """
        pass

    def add_molecule(self, molecule, local_topology_to_reference_index=None):
        """Add a Molecule to the Topology.

        Parameters
        ----------
        molecule : Molecule
            The Molecule to be added.
        local_topology_to_reference_index: dict, optional, default = None
            Dictionary of {TopologyMolecule_atom_index : Molecule_atom_index} for the TopologyMolecule that will be built

        Returns
        -------
        index : int
            The index of this molecule in the topology
        """
        from openforcefield.topology.molecule import FrozenMolecule
        #molecule.set_aromaticity_model(self._aromaticity_model)

        mol_smiles = molecule.to_smiles()
        reference_molecule = None
        for potential_ref_mol in self._reference_molecule_to_topology_molecules.keys(
        ):
            if mol_smiles == potential_ref_mol.to_smiles():
                # If the molecule is already in the Topology.reference_molecules, add another reference to it in
                # Topology.molecules
                reference_molecule = potential_ref_mol
                break
        if reference_molecule is None:
            # If it's a new unique molecule, make and store an immutable copy of it
            reference_molecule = FrozenMolecule(molecule)
            self._reference_molecule_to_topology_molecules[
                reference_molecule] = list()

        topology_molecule = TopologyMolecule(
            reference_molecule, self, local_topology_to_reference_index)
        self._topology_molecules.append(topology_molecule)
        self._reference_molecule_to_topology_molecules[
            reference_molecule].append(self._topology_molecules[-1])

        index = len(self._topology_molecules) - 1
        return index

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
        if (iatom, jatom) in self._constrained_atom_pairs:
            existing_distance = self._constrained_atom_pairs[(iatom, jatom)]
            if unit.is_quantity(existing_distance) and (distance is True):
                raise Exception(
                    'Atoms (%d,%d) already constrained with distance %s but attempting to override with unspecified distance'
                    % (iatom, jatom, existing_distance))
            if (existing_distance is True) and (distance is True):
                raise Exception(
                    'Atoms (%d,%d) already constrained with unspecified distance but attempting to override with unspecified distance'
                    % (iatom, jatom))
            if distance is False:
                del self._constrained_atom_pairs[(iatom, jatom)]
                del self._constrained_atom_pairs[(jatom, iatom)]
                return

        self._constrained_atom_pairs[(iatom, jatom)] = distance
        self._constrained_atom_pairs[(jatom, iatom)] = distance

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
        if (iatom, jatom) in self._constrained_atom_pairs:
            return self._constrained_atom_pairs[(iatom, jatom)]
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
