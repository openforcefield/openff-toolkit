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

import copy
import itertools

import collections

import numpy as np

from simtk import openmm, unit
from simtk.openmm.app import element as elem
from simtk.openmm import app

#from openforcefield.utils import get_data_filename
from openforcefield.typing.chemistry import ChemicalEnvironment, SMIRKSParsingError
from openforcefield.topology.molecule import DEFAULT_AROMATICITY_MODEL, DEFAULT_FRACTIONAL_BONDORDER_MODEL, DEFAULT_CHARGE_MODEL

from openforcefield.utils.serialization import Serializable

#=============================================================================================
# GLOBAL PARAMETERS
#=============================================================================================

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

class _TransformedDict(collections.MutableMapping):
    """A dictionary that applies an arbitrary key-altering
       function before accessing the keys"""

    def __init__(self, *args, **kwargs):
        self.store = collections.OrderedDict()
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

class ValenceDict(_TransformedDict):
    """Enforce uniqueness in atom indices"""
    def __keytransform__(self, key):
        """Reverse tuple if first element is larger than last element."""
        # Ensure key is a tuple.
        key = tuple(key)
        # Reverse the key if the first element is bigger than the last.
        if key[0] > key[-1]:
            key = tuple(reversed(key))
        return key

class ImproperDict(_TransformedDict):
    """Symmetrize improper torsions"""
    def __keytransform__(self,key):
        """Reorder tuple in numerical order except for element[1] which is the central atom; it retains its position."""
        # Ensure key is a tuple
        key = tuple(key)
        # Retrieve connected atoms
        connectedatoms = [key[0], key[2], key[3]]
        # Sort connected atoms
        connectedatoms.sort()
        # Re-store connected atoms
        key = tuple( [connectedatoms[0], key[1], connectedatoms[1], connectedatoms[2]])
        return(key)

#=============================================================================================
# TOPOLOGY OBJECTS
#=============================================================================================

#=============================================================================================
# Particle
#=============================================================================================

class Particle(object):
    """
    Base class for all particles in a molecule.

    A particle object could be an ``Atom`` or a ``VirtualSite``.

    """
    def __init__(self, name):
        """
        Create a particle.
        """
        self._name = name # the particle name
        self._topology = None # the Topology object this Particle belongs to

    @property
    def topology(self):
        """
        The Topology object that owns this particle, or None.
        """
        return self._topology

    @property
    def name(self):
        """
        An arbitrary label assigned to the particle.

        """
        return self._name

    @property
    def particle_index(self):
        """
        Index of this particle within the ``Topology`` or corresponding OpenMM ``System`` object.

        .. todo::

           Should ``atom.particle_index`` just be called ``index``, or does that risk confusion within
           the index within ``topology.atoms``, which will differ if the system has virtual sites?

        """
        if self._topology is None:
            raise Exception('This particle does not belong to a Topology')
        # Return index of this particle within the Topology
        # TODO: This will be slow; can we cache this and update it only when needed?
        #       Deleting atoms/molecules in the Topology would have to invalidate the cached index.
        return self._topology.particles.index(self)

    def __repr__(self):
        pass

    def __str__(self):
        pass

#=============================================================================================
# Atom
#=============================================================================================

class Atom(Particle):
    """
    A particle representing a chemical atom.

    Note that non-chemical virtual sites are represented by the ``VirtualSite`` object.

    .. todo::
        * Should ``Atom`` objects be immutable or mutable?
        * Should an ``Atom`` be able to belong to more than one ``Topology`` object?
        * Do we want to support the addition of arbitrary additional properties,
          such as floating point quantities (e.g. ``charge``), integral quantities (such as ``id`` or ``serial`` index in a PDB file),
          or string labels (such as Lennard-Jones types)?
        * Should we be able to create ``Atom`` objects on their own, or only in the context of a ``Topology`` object they belong to?

    """
    def __init__(self, name, element, topology=None):
        """
        Create an Atom object.

        Parameters
        ----------
        name : str
            A unique name for this atom
        element : str
            The element name

        """
        super(Atom, self).__init__(name)
        self._element = element # TODO: Validate and store Element

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

        """
        if self._topology is None:
            raise ValueError('This Atom does not belong to a Topology object')
        # TODO: This will be slow; can we cache this and update it only when needed?
        #       Deleting atoms/molecules in the Topology would have to invalidate the cached index.
        return self._topology.atoms.index(self)

    def __repr__(self):
        # TODO: Also include particle_index and which topology this atom belongs to?
        return "Atom(name={}, element={})".format(self.name, self.element)

    def __str__(self):
        # TODO: Also include particle_index and which topology this atom belongs to?
        return "<Atom name='{}' element='{}'>".format(self.name, self.element)

#=============================================================================================
# VirtualSite
#=============================================================================================

class VirtualSite(Particle):
    """
    A particle representing a virtual site whose position is defined in terms of ``Atom`` positions.

    Note that chemical atoms are represented by the ``Atom``.

    .. todo::
        * Should a virtual site be able to belong to more than one Topology?
        * Should virtual sites be immutable or mutable?

    """

    # TODO: This will need to be generalized for virtual sites to allow out-of-plane sites.
    # TODO: How do we want users to specify virtual site type?
    def __init__(self, name, sitetype, weights, atoms):
        """
        Create a virtual site whose position is defined by a linear combination of multiple Atoms.

        Parameters
        ----------
        name : str
            The name of this virtual site
        sitetype : str
            The virtual site type.
        weights : list of floats of shape [N]
            weights[index] is the weight of particles[index] contributing to the position of the virtual site.
        atoms : list of Atom of shape [N]
            atoms[index] is the corresponding Atom for weights[index]
        virtual_site_type : str
            Virtual site type.
            TODO: What types are allowed?

        """
        self._name = name
        self._type = sitetype # TODO: Validate site types against allowed values
        self._weights = np.array(weights) # make a copy and convert to array internally
        self._atoms = [ atom for atom in atoms ] # create a list of Particles

    @property
    def virtual_site_index(self):
        """
        The index of this VirtualSite within the list of virtual sites within ``Topology``
        Note that this can be different from ``particle_index``.

        """
        if self._topology is None:
            raise ValueError('This VirtualSite does not belong to a Topology object')
        # TODO: This will be slow; can we cache this and update it only when needed?
        #       Deleting atoms/molecules in the Topology would have to invalidate the cached index.
        return self._topology.virtual_sites.index(self)

    @property
    def atoms(self):
        """
        Atoms on whose position this VirtualSite depends.
        """
        for atom in self._atoms:
            yield atom

    def __repr__(self):
        # TODO: Also include particle_index, which topology this atom belongs to, and which atoms/weights it is defined by?
        return "VirtualSite(name={}, type={}, weights={}, atoms={})".format(self.name, self.type, self.weights, self.atoms)

    def __str__(self):
        # TODO: Also include particle_index, which topology this atom belongs to, and which atoms/weights it is defined by?
        return "<VirtualSite name={} type={} weights={}, atoms={}>".format(self.name, self.type, self.weights, self.atoms)

#=============================================================================================
# Bond
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

#=============================================================================================
# TOPOLOGY
#=============================================================================================

class Topology(Serializable):
    """
    A Topology is a chemical representation of a system containing one or more molecules appearing in a specified order.

    Attributes
    ----------
    chains : list of Chain
        Iterate over all Chain objects in the topology
    molecules : list of Molecule
        Iterate over all Molecule objects in the system in the topology
    unique_molecules : list of Molecule
        Iterate over all unique Molecule objects in the topology

    n_chains : int
        Number of chains in the topology
    n_molecules : int
        Number of molecules in the topology
    n_unique_molecules : int
        Number of unique molecules in the topology

    Examples
    --------

    Import some utilities

    >>> from simtk.openmm import app
    >>> from openforcefield.tests.utils.utils import get_monomer_mol2file, get_packmol_pdbfile
    >>> pdb_filename = get_packmol_pdbfile('cyclohexane_ethanol_0.4_0.6.pdb')
    >>> monomer_names = ('cyclohexane', 'ethanol')

    Create a Topology object from a PDB file and mol2 files defining the molecular contents

    >>> pdbfile = app.PDBFile(pdb_filename)
    >>> mol2_filenames = [ get_monomer_mol2file(name) for name in monomer_names ]
    >>> unique_molecules = [ Molecule.from_file(mol2_filename) for mol2_filename in mol2_filenames ]
    >>> topology = Topology.from_openmm(pdbfile.topology, unique_molecules=unique_molecules)

    Create a Topology object from a PDB file and IUPAC names of the molecular contents

    >>> pdbfile = app.PDBFile(pdb_filename)
    >>> unique_molecules = [ Molecule.from_iupac(name) for name in monomer_names ]
    >>> topology = Topology.from_openmm(pdbfile.topology, unique_molecules=unique_molecules)

    Create an empty Topology object and add a few copies of a single benzene molecule

    >>> topology = Topology()
    >>> molecule = Molecule.from_iupac('benzene')
    >>> [ topology.add_molecule(molecule) for index in range(10) ]

    Create a deep copy of the Topology and its contents

    >>> topology_copy = Topology(topology)

    Create a Topology from an OpenEye Molecule, including perception of chains and residues
    (requires the OpenEye toolkit)

    >>> oemol = oechem.oemolistream('input.pdb')
    >>> topology = Topology.from_openeye(oemol)

    .. todo ::

       Should the :class:`Topology` object be able to have optional positions and box vectors?
       If so, this would make the creation of input files for other molecular simulation packages much easier.

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
        # Assign cheminformatics models
        self._aromaticity_model = DEFAULT_AROMATICITY_MODEL
        self._fractional_bondorder_model = DEFAULT_FRACTIONAL_BONDORDER_MODEL
        self._charge_model = DEFAULT_CHARGE_MODEL

        # Initialize storage
        self._constrained_atom_pairs = dict()

        # TODO: Try to construct Topology copy from `other` if specified
        pass

    @staticmethod
    def from_molecules(molecules):
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
            some_object_iterator = iter(files)
        except TypeError as te:
            # Make iterable object
            files = [files]

        # Create Topology and populate it with specified molecules
        topology = Topology()
        for molecule in molecules:
            topology.add_molecule(molecule)

        return topology

    def assert_bonded(atom1, atom2):
        """
        Raise an exception if the specified atoms are not bonded in the topology.

        Parameters
        ----------
        atom1, atom2 : openforcefield.topology.Atom
            The atoms to check to ensure they are bonded

        .. todo :: Should we move this to ForceField?

        """
        # TODO: Should atom1 and atom2 be int or Atom objects?
        assert self.is_bonded(atom1, atom2), 'Atoms {} and {} are not bonded in topology'.format(atom1, atom2)

    def set_aromaticity_model(self, aromaticity_model):
        """
        Set the aromaticity model applied to all molecules in the topology.

        Parameters
        ----------
        aromaticity_model : str
            Aromaticity model to use. One of: ['MDL']

        """
        if not aromaticity_model in ALLOWED_AROMATICITY_MODELS:
            raise ValueError("Aromaticity model must be one of {}; specified '{}'".format(ALLOWED_AROMATICITY_MODELS, aromaticity_model))
        self._aromaticity_model = aromaticity_model

    def get_aromaticity_model(self):
        """
        Get the aromaticity model applied to all molecules in the topology.

        Returns
        -------
        aromaticity_model : str
            Aromaticity model in use.

        """
        return self._aromaticity_model

    def set_fractional_bondorder_model(self, fractional_bondorder_model):
        """
        Set the fractional bond order model applied to all molecules in the topology.

        Parameters
        ----------
        fractional_bondorder_model : str
            Fractional bond order model to use. One of: ['Wiberg']

        """
        if not fractional_bondorder_model in ALLOWED_FRACTIONAL_BONDORDER_MODELS:
            raise ValueError("Fractional bond order model must be one of {}; specified '{}'".format(ALLOWED_FRACTIONAL_BONDORDER_MODELS, fractional_bondorder_model))
        self._fractional_bondorder_model = fractional_bondorder_model

    def get_fractional_bond_order(self):
        """
        Get the fractional bond order model for the Topology.

        Returns
        -------
        fractional_bondorder_model : str
            Fractional bond order model in use.

        """
        return self._fractional_bondorder_model

    def set_charge_model(self, charge_model):
        """
        Set the fractional bond order model applied to all molecules in the topology.

        Parameters
        ----------
        charge_model : str
            Charge model to use for all molecules in the Topology.
            Allowed values: ['AM1-CM2', 'AM1-BCC', 'Mulliken']
            * ``AM1-CM2``: AM1 wavefunction with CM2 population analysis
            * ``AM1-BCC``: Canonical AM1-BCC scheme
            * ``Mulliken``: Mulliken charges
        """
        if not charge_model in ALLOWED_CHARGE_MODELS:
            raise ValueError("Charge model must be one of {}; specified '{}'".format(ALLOWED_CHARGE_MODELS, charge_model))
        self._charge_model = charge_model

    def chemical_environment_matches(self, query, aromaticity_model='MDL'):
        """Retrieve all matches for a given chemical environment query.

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
        matches : list of Atom tuples
            A list of all matching Atom tuples

        """
        # Render the query to a SMARTS string
        if type(query) is str:
            smarts = query
        elif type(query) is ChemicalEnvironment:
            smarts = query.as_smarts()
        else:
            raise ValueError("Don't know how to convert query '%s' into SMARTS string" % query)

        # Perform matching on each unique molecule, unrolling the matches to all matching copies of that molecule in the Topology object.
        matches = list()
        for molecule in self.unique_molecules:
            # Find all atomsets that match this definition in the reference molecule
            # This will automatically attempt to match chemically identical atoms in a canonical order within the Topology
            refmol_matches = molecule.chemical_environment_matches(smarts)

            # Loop over matches
            for reference_match in refmol_matches:
                # Unroll corresponding atom indices over all instances of this molecule
                for reference_to_topology_atom_mapping in self._reference_to_topology_atom_mappings[reference_molecule]:
                    # Create match.
                    match = tuple([ reference_to_topology_atom_mapping[atom] for atom in reference_match ])
                    matches.append(match)

        return matches

    def chemical_environment_matches(self, query):
        """Retrieve all matches for a given chemical environment query.

        .. todo ::

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
            refmol_matches = molecule.chemical_environment_matches(query)

            # Loop over matches
            for reference_atom_indices in refmol_matches:
                # Unroll corresponding atom indices over all instances of this molecule
                # TODO: This is now handled through MoleculeImages
                for reference_to_topology_atom_mapping in self._reference_to_topology_atom_mappings[reference_molecule]:
                    # Create match.
                    atom_indices = tuple([ reference_to_topology_atom_mapping[atom_index] for atom_index in reference_atom_indices ])
                    matches.append(atom_indices)

        return matches

    # TODO: Overhaul this function so that we identify molecules as they are added to the Topology
    # TODO: We also need to ensure the order of atoms is matched the same way for each unique_molecule if possible
    def _identify_molecules(self):
        """Identify all unique reference molecules and atom mappings to all instances in the Topology.
        """
        import networkx as nx

        # Generate list of topology atoms.
        atoms = [ atom for atom in self.atoms() ]

        # Generate graphs for reference molecules.
        self._reference_molecule_graphs = list()
        for reference_molecule in self._reference_molecules:
            # Generate Topology
            reference_molecule_topology = generateTopologyFromOEMol(reference_molecule)
            # Generate Graph
            reference_molecule_graph = reference_molecule_topology.to_networkx()
            self._reference_molecule_graphs.append(reference_molecule_graph)

        # Generate a graph for the current topology.
        G = self.to_networkx()

        # Extract molecules (as connected component subgraphs).
        from networkx.algorithms import isomorphism
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

    @staticmethod
    def from_openmm(openmm_topology, unique_molecules=None):
        """
        Construct an openforcefield Topology object from an OpenMM Topology object.

        Parameters
        ----------
        openmm_topology : simtk.openmm.app.Topology
            An OpenMM Topology object
        unique_molecules : iterable of objects that can be used to construct unique Molecule objects
            All unique molecules mult be provided, in any order, though multiple copies of each molecule are allowed.
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
        raise NotImplementedError

    def to_openmm(self):
        """
        Create an OpenMM Topology object.

        Parameters
        ----------
        openmm_topology : simtk.openmm.app.Topology
            An OpenMM Topology object
        """
        raise NotImplementedError

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
        return Topology.from_openmm(mdtraj_topology.to_openmm(), unique_molecules=unique_molecules)

    def to_mdtraj(self):
        """
        Create an MDTraj Topology object.

        Returns
        ----------
        mdtraj_topology : mdtraj.Topology
            An MDTraj Topology object
        """
        return md.Topology.from_openmm(self.to_openmm())

    @staticmethod
    def from_parmed(parmed_structure, unique_molecules=None):
        """
        Construct an openforcefield Topology object from a ParmEd Structure object.

        Parameters
        ----------
        mdtraj_topology : mdtraj.Topology
            An MDTraj Topology object
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
        Create a ParmEd Structure object.

        Returns
        ----------
        parmed_structure : parmed.Structure
            A ParmEd Structure objecft
        """
        import parmed
        # TODO: Implement functionality
        raise NotImplementedError


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
        hv = oechem.OEHierView(oemol, oechem.OEAssumption_BondedResidue +
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

    def to_openeye(self, positions=None, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
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
        if (iatom, jatom) in self._constrained_atom_pairs:
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

    @property
    def is_periodic(self):
        """
        ``True`` if the topology represents a periodic system; ``False`` otherwise
        """
        return self._is_periodic
