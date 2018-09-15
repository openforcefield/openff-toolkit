#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Molecular chemical entity representation and routines to interface with cheminformatics toolkits

.. todo::

   * Our main philosophy here is to keep the object contents of topology objects easily serializable/deserializable

   * Have ``Molecule`` raise an exception if loading/creating molecules with unspecified stereochemistry?
   * Create ``FrozenMolecule`` to represent immutable molecule
   * Make ``Atom`` and ``Bond`` an inner class of Molecule?
   * Add ``Molecule.from_smarts()`` or ``.from_tagged_smiles()`` to allow a tagged SMARTS string
     (where tags are zero-indexed atom indices) to be used to create a molecule with the given atom numbering.
   * How can we make the ``Molecule`` API more useful to codes like perses that modify molecules on the fly?
   * Use `attrs <http://www.attrs.org/>`_ for convenient class initialization?
   * JSON/BSON representations of objects?
   * Generalize Molecule infrastructure to provide "plug-in" support for cheminformatics toolkits
   * Do we need a way to write a bunch of molecules to a file, or serialize a set of molecules to a file?
     We currently don't have a way to do that through the ``Molecule`` API, even though there is a way to
     read multiple molecules via ``Molecules.from_file()``.
   * Should we allow the removal of atoms too?
   * Should invalidation of cached properties be handled via something like a tracked list?
   * Refactor toolkit encapsulation to generalize and provide only a few major toolkit methods and toolkit objects that can be queried for features
   * Speed up overall import time by putting non-global imports only where they are needed

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import numpy as np
from copy import deepcopy

from distutils.spawn import find_executable

from simtk import unit
from simtk.openmm.app import element

import openforcefield
from openforcefield.utils.toolkits import OPENEYE_INSTALLED, RDKIT_INSTALLED, TOOLKIT_PRECEDENCE, SUPPORTED_FILE_FORMATS
from openforcefield.utils.toolkits import requires_rdkit, requires_openeye
from openforcefield.typing.chemistry import ChemicalEnvironment, SMIRKSParsingError

# TODO: Do we need these?
from openforcefield.utils.toolkits import DEFAULT_AROMATICITY_MODEL, ALLOWED_AROMATICITY_MODELS
from openforcefield.utils.toolkits import DEFAULT_FRACTIONAL_BONDORDER_MODEL, ALLOWED_FRACTIONAL_BONDORDER_MODELS
from openforcefield.utils.toolkits import DEFAULT_CHARGE_MODEL, ALLOWED_CHARGE_MODELS

from openforcefield.utils.serialization import Serializable

#=============================================================================================
# GLOBAL PARAMETERS
#=============================================================================================

# TODO: Move these to utils.toolkits?

# TODO: Can we have the `ALLOWED_*_MODELS` list automatically appear in the docstrings below?
# TODO: Should `ALLOWED_*_MODELS` be objects instead of strings?
# TODO: Should these be imported from `openforcefield.cheminformatics.aromaticity_models` and `.bondorder_models`?

# TODO: Allow all OpenEye aromaticity models to be used with OpenEye names?
#       Only support OEAroModel_MDL in RDKit version?

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

#=============================================================================================
# Particle
#=============================================================================================

class Particle(Serializable):
    """
    Base class for all particles in a molecule.

    .. warning :: This API experimental and subject to change.

    A particle object could be an ``Atom`` or a ``VirtualSite``.

    """
    pass

#=============================================================================================
# Atom
#=============================================================================================

class Atom(Particle):
    """
    A particle representing a chemical atom.

    .. warning :: This API experimental and subject to change.

    Note that non-chemical virtual sites are represented by the ``VirtualSite`` object.

    .. todo::

       * Should ``Atom`` objects be immutable or mutable?
       * Should an ``Atom`` be able to belong to more than one ``Topology`` object?
       * Do we want to support the addition of arbitrary additional properties,
         such as floating point quantities (e.g. ``charge``), integral quantities (such as ``id`` or ``serial`` index in a PDB file),
         or string labels (such as Lennard-Jones types)?
       * Should we be able to create ``Atom`` objects on their own, or only in the context of a ``Topology`` object they belong to?

    .. todo :: Allow atoms to have associated properties.

    """
    def __init__(self, atomic_number, formal_charge, is_aromatic, stereochemistry=None, name=None):
        """
        Create an immutable Atom object.

        Object is serializable and immutable.

        .. todo :: Use attrs to validate?

        .. todo :: We can add setters if we need to.

        Parameters
        ----------
        atomic_number : int
            Atomic number of the atom
        formal_charge : int
            Formal charge of the atom
        is_aromatic : bool
            If True, atom is aromatic; if False, not aromatic
        stereochemistry : str, optional, default=None
            Either 'R' or 'S' for specified stereochemistry, or None for ambiguous stereochemistry
        name : str, optional, default=None
            An optional name to be associated with the atom

        Examples
        --------

        Create a non-aromatic carbon atom

        >>> atom = Atom(6, 0, False)

        Create a chiral carbon atom

        >>> atom = Atom(6, 0, False, stereochemistry='R', name='CT')

        """
        self._atomic_number = atomic_number
        self._formal_charge = formal_charge
        self._is_aromatic = is_aromatic
        self._stereochemistry = stereochemistry
        self._name = name

    def to_dict(self):
        """Return a dict representation of the atom."""
        # TODO
        return self.__dict__

    @staticmethod
    def from_dict(atom_dict):
        """Create an Atom from a dict representation."""
        # TODO
        return Atom(**atom_dict)

    @property
    def element(self):
        """
        The element name

        """
        return element.Element.getByAtomicNumber(self._atomic_number)

    @property
    def atomic_number(self):
        """
        The integer atomic number of the atom.

        """
        return self._atomic_number

    @property
    def mass(self):
        """
        The standard atomic weight (abundance-weighted isotopic mass) of the atomic site.

        .. todo :: Should we discriminate between standard atomic weight and most abundant isotopic mass?

        """
        return self.element.mass

    # TODO: How are we keeping track of bonds, angles, etc?

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
        # TODO: Also include particle_index and which molecule this atom belongs to?
        return "Atom(name={}, element={})".format(self.name, self.element)

    def __str__(self):
        # TODO: Also include particle_index and which molecule this atom belongs to?
        return "<Atom name='{}' element='{}'>".format(self.name, self.element)

#=============================================================================================
# VirtualSite
#=============================================================================================

class VirtualSite(Particle):
    """
    A particle representing a virtual site whose position is defined in terms of ``Atom`` positions.

    .. warning :: This API experimental and subject to change.

    Note that chemical atoms are represented by the ``Atom``.

    .. todo::

       * Should a virtual site be able to belong to more than one Topology?
       * Should virtual sites be immutable or mutable?

    """

    def __init__(self, name, sitetype, weights, atoms):
        """
        Create a virtual site whose position is defined by a linear combination of multiple Atoms.

        .. todo ::

           * This will need to be gneeralized for virtual sites to allow out-of-plane sites, which are not simply a linear combination of atomic positions
           * Add checking for allowed virtual site types
           * How do we want users to specify virtual site types?

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

        """
        self._name = name
        self._type = sitetype # TODO: Validate site types against allowed values
        self._weights = np.array(weights) # make a copy and convert to array internally
        self._atoms = [ atom for atom in atoms ] # create a list of Particles

    def to_dict(self):
        """Return a dict representation of the atom."""
        # TODO
        return self.__dict__

    @staticmethod
    def from_dict(vsite_dict):
        """Create an Atom from a dict representation."""
        # TODO
        return VirtualSite(**vsite_dict)

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
        # TODO: Also include particle_index, which molecule this atom belongs to?
        return "VirtualSite(name={}, type={}, weights={}, atoms={})".format(self.name, self.type, self.weights, self.atoms)

    def __str__(self):
        # TODO: Also include particle_index, which molecule this atom belongs to?
        return "<VirtualSite name={} type={} weights={}, atoms={}>".format(self.name, self.type, self.weights, self.atoms)

#=============================================================================================
# Bond
#=============================================================================================

class Bond(Serializable):
    """
    Chemical bond representation.

    .. warning :: This API experimental and subject to change.

    .. todo :: Allow bonds to have associated properties.

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
        assert type(atom1) == Atom
        assert type(atom2) == Atom
        self._atom1 = atom1
        self._atom2 = atom2
        # TODO: Check bondtype and fractional_bondorder are valid?
        self._type = bondtype
        self._fractional_bondorder = fractional_bondorder

    def to_dict(self):
        """Return a dict representation of the bond."""
        # TODO
        return self.__dict__

    @staticmethod
    def from_dict(d):
        """Create an Bond from a dict representation."""
        # TODO
        return Bond(**d)

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
# Molecule
#=============================================================================================

# TODO: Make Molecule immutable (by default)

# TODO: How do we automatically trigger invalidation of cached properties if an ``Atom``, ``Bond``, or ``VirtualSite`` is modified,
#       rather than added/deleted via the API? The simplest resolution is simply to make them immutable.

class Molecule(Serializable):
    """
    Chemical representation of a molecule, such as a small molecule or biopolymer.

    .. todo :: What other API calls would be useful for supporting biopolymers as small molecules? Perhaps iterating over chains and residues?

    Examples
    --------

    Create a molecule from a mol2 file

    >>> molecule = Molecule.from_file('molecule.mol2')

    Create a molecule from an OpenEye molecule

    >>> molecule = Molecule.from_openeye(oemol)

    Create a molecule from an RDKit molecule

    >>> molecule = Molecule.from_rdkit(rdmol)

    Create a molecule from IUPAC name (requires the OpenEye toolkit)

    >>> molecule = Molecule.from_iupac('imatinib')

    Create a molecule from SMILES

    >>> molecule = Molecule.from_smiles('Cc1ccccc1')

    """
    def __init__(self, other=None):
        """
        Create a new Molecule object

        .. todo ::

           * If a filename or file-like object is specified but the file contains more than one molecule, what is the proper behavior?
           Read just the first molecule, or raise an exception if more than one molecule is found?

           * Should we also support SMILES strings or IUPAC names for ``other``?

        Parameters
        ----------
        other : optional, default=None
            If specified, attempt to construct a copy of the Molecule from the specified object.
            This can be any one of the following:

            * a :class:`Molecule` object
            * a file that can be used to construct a :class:`Molecule` object
            * an ``openeye.oechem.OEMol``
            * an ``rdkit.Chem.rdchem.Mol``
            * a serialized :class:`Molecule` object

        Examples
        --------

        Create an empty molecule:

        >>> empty_molecule = Molecule()

        Create a molecule from another molecule:

        >>> molecule_copy = Molecule(molecule_copy)

        Create a molecule from a file that can be used to construct a molecule,
        using either a filename or file-like object:

        >>> molecule = Molecule('molecule.mol2')
        >>> molecule = Molecule(open('molecule.mol2', 'r'))
        >>> molecule = Molecule(gzip.GzipFile('molecule.mol2.gz', 'r'))

        Create a molecule from an OpenEye molecule:

        >>> molecule = Molecule(oemol)

        Create a molecule from an RDKit molecule:

        >>> molecule = Molecule(rdmol)

        Create a molecule from a serialized molecule object:

        >>> serialized_molecule = molecule.__getstate__()
        >>> molecule_copy = Molecule(serialized_molecule)

        """
        if other is None:
            self._initialize()
        else:
            # TODO: Can we check interface compliance (in a try..except) instead of checking instances?
            if isinstance(other, openforcefield.topology.Molecule):
                self._copy_initializer(other)
            elif isinstance(other, str):
                self.__setstate__(other)
            elif OPENEYE_INSTALLED and issubclass(other, openeye.oechem.OEMolBase):
                mol = Molecule.from_openeye(other)
                self._copy_initializer(mol)
            elif RDKIT_INSTALLED and isinstance(other, rdkit.Chem.rdchem.Mol):
                mol = Molecule.from_rdkit(other)
                self._copy_initializer(mol)
            elif isinstance(other, str) or hasattr(other, 'read'):
                mol = Molecule.from_file(other) # returns a list only if multiple molecules are found
                if type(mol) == list:
                    raise ValueError('Specified file or file-like object must contain exactly one molecule')
                self._copy_initializer(molecule)
            else:
                msg = 'Cannot construct openforcefield.topology.Molecule from {}\n'.format(other)
                raise Exception(msg)

    ####################################################################################################
    # Safe serialization
    ####################################################################################################

    def to_dict(self):
        """
        Return a dictionary representation of the molecule.

        .. todo ::

           * Document the representation standard.
           * How do we do version control with this standard?

        Returns
        -------
        molecule_dict : OrderedDict
            A dictionary representation of the molecule.

        """
        molecule_dict = OrderedDict()
        molecule_dict['name'] = self.name
        molecule_dict['atoms'] = [ atom.to_dict() for atom in self.atoms ]
        molecule_dict['virtual_sites'] = [ vsite.to_dict() for vsite in self.virtual_sites ]
        molecule_dict['bonds'] = [ bond.to_dict() for bond in self.bonds ] # TODO: How do we make sure bonds are serializable?
        # TODO: Charges
        # TODO: Properties
        # TODO: Conformers
        return molecule_dict

    @staticmethod
    def from_dict(molecule_dict):
        """
        Create a Molecule from a dictionary representation

        Parameters
        ----------
        molecule_dict : OrderedDict
            A dictionary representation of the molecule.

        Returns
        -------
        molecule : Molecule
            A Molecule created from the dictionary representation

        """
        # TODO: Provide useful exception messages if there are any failures
        molecule = Molecule(name=molecule_dict['name'])
        for atom in molecule_dict['atoms']:
            molecule.add_atom(*atom)
        for vsite in molecule_dict['virtual_sites']:
            molecule.add_atom(*vsite)
        for bond in molecule_dict['bonds']:
            molecule.add_bond(*bond) # TODO: How is this correctly handled? using indices?
        # TODO: Charges
        # TODO: Properties
        # TODO: Conformers
        return molecule

    def __getstate__(self):
        return self.to_dict()

    def __setstate__(self, state):
        return self.from_dict(state)

    def _initialize(self):
        """
        Clear the contents of the current molecule.
        """
        self._name = None # TODO: Should we keep a name, or just store that in _properties?
        self._particles = list() # List of particles (atoms or virtual sites) # TODO: Should this be a dict?
        self._bonds = list() # List of bonds between Atom objects # TODO: Should this be a dict?
        self._charges = None # TODO: Storage charges
        self._properties = None # Attached properties to be preserved
        self._cached_properties = None # Cached properties (such as partial charges) can be recomputed as needed
        self._conformers = None # Optional conformers

    def _copy_initializer(self, other):
        """
        Copy contents of the specified molecule

        .. todo :: Should this be a ``@staticmethod`` where we have an explicit copy constructor?

        Parameters
        ----------
        other : optional
            Overwrite the state of this Molecule with the specified Molecule object.
            A deep copy is made.

        """
        assert isinstance(other, type(self)), "can only copy instances of {}".format(typse(self))
        self.__dict__ = deepcopy(other.__dict__)

    def _eq_(self, other):
        """Test two molecules for equality to see if they are the chemical species, but do not check other annotated properties.

        .. note ::

           Note that this method simply tests whether two molecules are identical chemical species using equivalence of their canonical isomeric SMILES.
           No effort is made to ensure that the atoms are in the same order or that any annotated properties are preserved.

        """
        return self.to_smiles() == other.to_smiles()

    def to_smiles(self):
        """
        Return a canonical isomeric SMILES representation of the current molecule

        .. note :: RDKit and OpenEye versions will not necessarily return the same representation.

        .. todo :: Can we ensure RDKit and OpenEye versions return the same representation?

        Returns
        -------
        smiles : str
            Canonical isomeric explicit-hydrogen SMILES

        Examples
        --------

        >>> smiles = molecule.to_smiles()

        """
        # TODO: Rework this toolkit usage to use a standard ToolkitWrapper interface
        for toolkit in REGISTERED_TOOLKITS:
            try:
                return toolkit.to_smiles(self)
            except NotImplementedError as e:
                pass
        raise NotImplementedError('No registered toolkits can provide this capability.')

        # TODO: Is there a simmpler scheme that could streamline this?
        return toolkit_handler.call('to_smiles', self)

        # Old code that does not use a toolkit wrapper
        toolkit = TOOLKIT_PRECEDENCE[0]
        if toolkit == 'oechem':
            oemol = self.to_openeye()
            return oechem.OEMolToSmiles(oemol)
        elif toolkit == 'rdkit':
            rdmol = self.to_rdkit()
            return Chem.MolToSmiles(rdmol, isomericSmiles=True)
        else:
            raise Exception('Unknown toolkit {}'.format(toolkit))

    @staticmethod
    def from_smiles(smiles):
        """
        Construct a Molecule from a SMILES representation

        Parameters
        ----------
        smiles : str
            The SMILES representation of the molecule.

        Returns
        -------
        molecule : Molecule
            The molecule

        Examples
        --------

        >>> molecule = Molecule.from_smiles('Cc1ccccc1')

        """
        toolkit = TOOLKIT_PRECEDENCE[0]
        if toolkit == 'oechem':
            oemol = oechem.OEMol()
            if not oechem.OESmilesToMol(oemol, smiles):
                raise ValueError("Could not parse SMILES string: {}".format(smiles))
            return Molecule.from_openeye(oemol)
        elif toolkit == 'rdkit':
            rdmol = Chem.MolFromSmiles(smiles)
            if rdmol is None:
                raise ValueError("Could not parse SMILES string: {}".format(smiles))
            return Molecule.from_rdkit(rdmol)
        else:
            raise Exception('Unknown toolkit {}'.format(toolkit))

    def _invalidate_cached_properties(self):
        """
        Indicate that the chemical entity has been altered.
        """
        if hasattr(self, '_cached_properties'):
            delattr(self, '_cached_properties')

    def to_networkx(self):
        """Geneate a NetworkX undirected graph from the Topology.

        Nodes are Atoms labeled with particle indices and atomic elements (via the ``element`` node atrribute).
        Edges denote chemical bonds between Atoms.
        Virtual sites are not included, since they lack a concept of chemical connectivity.

        .. todo ::

           * Do we need a ``from_networkx()`` method? If so, what would the Graph be required to provide?
           * Should edges be labeled with discrete bond types in some aromaticity model?
           * Should edges be labeled with fractional bond order if a method is specified?
           * Should we add other per-atom and per-bond properties (e.g. partial charges) if present?

        Returns
        -------
        graph : networkx.Graph
            The resulting graph, with nodes labeled with atom indices and elements

        Examples
        --------
        Retrieve the bond graph for imatinib (OpenEye toolkit required)

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> nxgraph = molecule.to_networkx()

        """
        import networkx as nx
        G = nx.Graph()
        for atom in topology.atoms():
            G.add_node(atom.particle_index, element=atom.element)
        for (atom1, atom2) in topology.bonds():
            G.add_edge(atom1.index, atom2.index)

        return G

    def add_atom(atomic_number, formal_charge, is_aromatic, stereochemistry=None, name=None):
        """
        Add an atom

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        atomic_number : int
            Atomic number of the atom
        formal_charge : int
            Formal charge of the atom
        is_aromatic : bool
            If True, atom is aromatic; if False, not aromatic
        stereochemistry : str, optional, default=None
            Either 'R' or 'S' for specified stereochemistry, or None if stereochemistry is irrelevant
        name : str, optional, default=None
            An optional name for the atom

        Returns
        -------
        index : int
            The index of the atom in the molecule

        Examples
        --------

        Define a methane molecule

        >>> molecule = Molecule(name='methane')
        >>> C = molecule.add_atom(6, 0, False)
        >>> H1 = molecule.add_atom(1, 0, False)
        >>> H2 = molecule.add_atom(1, 0, False)
        >>> H3 = molecule.add_atom(1, 0, False)
        >>> H4 = molecule.add_atom(1, 0, False)
        >>> molecule.add_bond(C, H1, False, 1)
        >>> molecule.add_bond(C, H2, False, 1)
        >>> molecule.add_bond(C, H3, False, 1)
        >>> molecule.add_bond(C, H4, False, 1)

        """
        # Create an atom
        atom = Atom(atomic_number=atomic_number, formal_charge=formal_charge, is_aromatic=is_aromatic, stereochemistry=stereochemistry)
        self._atoms.append(atom)
        self._invalidate_cached_properties()

    def add_bond(self, atom1_index, atom2_index, is_aromatic, order, stereochemistry=None):
        """
        Add a bond between two specified atom indices

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        atom1_index : int
            Index of first atom
        atom2_index : int
            Index of second atom
        order : int
            Integral bond order of Kekulized form
        is_aromatic : bool
            True if this bond is aromatic, False otherwise
        stereochemistry : str, optional, default=None
            Either 'E' or 'Z' for specified stereochemistry, or None if stereochemistry is irrelevant

        """
        # TODO: Check to make sure bond does not already exist
        self._bonds.append(atom1, atom2)
        self._invalidate_cached_properties()

    def add_virtual_site(self, virtual_site):
        """
        Add a Virtual Site.

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        virtual_site : VirtualSite
            The VirtualSite to add.

        """
        # Make sure that all Atoms referenced in the virtual site are already present in the entity.
        for atom in virtual_site.atoms:
            if atom not in self._particles:
                raise Exception("{} depends on {}, which is not a mamber of the molecule".format(virtual_site, atom))
        self._particles.append(virtual_site)
        self._invalidate_cached_properties()

    @property
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
        for proper in self.propers:
            yield property
        for improper in self.impropers:
            yield improper

    @property
    def total_charge(self):
        """Return the total charge on the molecule"""
        return sum([atom.formal_charge for atom in self.atoms])

    def chemical_environment_matches(self, query):
        """Retrieve all matches for a given chemical environment query.

        .. todo ::

           * Do we want to generalize ``query`` to allow other kinds of queries, such as mdtraj DSL, pymol selections, atom index slices, etc?
             We could call it ``topology.matches(query)`` instead of ``chemical_environment_matches``

        Parameters
        ----------
        query : str or ChemicalEnvironment
            SMARTS string (with one or more tagged atoms) or ``ChemicalEnvironment`` query
            Query will internally be resolved to SMIRKS using ``query.asSMIRKS()`` if it has an ``.asSMIRKS`` method.

        Returns
        -------
        matches : list of Atom tuples
            A list of all matching Atom tuples

        Examples
        --------
        Retrieve all the carbon-carbon bond matches in a molecule

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> matches = molecule.chemical_environment_matches('[#6X3:1]~[#6X3:2]')

        """
        # Resolve to SMIRKS if needed
        # TODO: Update this to use updated ChemicalEnvironment API
        if hasattr(query, 'asSMIRKS'):
            smirks = query.asSMIRKS()
        elif type(query) == str:
            smirks = query
        else:
            raise ValueError("'query' must be either a string or a ChemicalEnvironment")

        # Use specified cheminformatics toolkit to determine matches with specified aromaticity model
        # TODO: Move this to toolkit handling
        toolkit = TOOLKIT_PRECEDENCE[0]
        if toolkit == 'oechem':
            matches = _openye_smirks_matches(self.to_openeye(), smirks, aromaticity_model=self._aromaticity_model)
        elif toolkit == 'rdkit':
            matches = _openeye_smirks_matches(self.to_rdkit(), smirks, aromaticity_model=self._aromaticity_model)
        else:
            raise Exception('Unknown toolkit {}'.format(toolkit))

        return matches

    # TODO: Move this to ToolkitWrapper
    @staticmethod
    @requires_rdkit()
    def _rdkit_smirks_matches(rdmol, smirks, aromaticity_model='OEAroModel_MDL'):
        """Find all sets of atoms in the provided RDKit molecule that match the provided SMARTS string.

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        rdmol : rdkit.Chem.Mol
            rdmol to process with the SMIRKS in order to find matches
        smarts : str
            SMARTS string with any number of sequentially tagged atoms.
            If there are N tagged atoms numbered 1..N, the resulting matches will be N-tuples of atoms that match the corresponding tagged atoms.
        aromaticity_model : str, optional, default='OEAroModel_MDL'
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

           * Raises ``ValueError`` if ``smarts`` query is malformed

        """
        from rdkit import Chem

        # Make a copy of the molecule
        rdmol = Chem.Mol(rdmol)
        # Use designated aromaticity model
        if aromaticity_model == 'OEAroModel_MDL':
            Chem.SanitizeMol(mol, Chem.SANITIZE_ALL^Chem.SANITIZE_SETAROMATICITY)
            Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        else:
            raise ValueError('Unknown aromaticity model: {}'.aromaticity_models)

        # Set up query.
        qmol = Chem.MolFromSmarts(smirks)   #cannot catch the error
        if qmol is None:
            raise SMIRKSParsingError('RDKit could not parse the SMIRKS string "{}"'.format(smirks))

        # Create atom mapping for query molecule
        index_map = dict()
        for atom in qmol.GetAtoms():
             smirks_index = atom.GetAtomMapNum()
             if smirks_index != 0:
                ind_map[smirks_index - 1] = atom.GetIdx()
        map_list = [ index_map[x] for x in sorted(index_map) ]

        # Perform matching
        # TODO: The MoleculeImage mapping should preserve ordering of template molecule for equivalent atoms
        #       and speed matching for larger molecules.
        matches = list()
        for match in rdmol.GetSubstructMatches(qmol, uniquify=False):
            mas = [ match[x] for x in map_list ]
            matches.append(tuple(mas))

        return matches

    # TODO: Move this to ToolkitWrapper
    @staticmethod
    @requires_openeye('oechem')
    def _oechem_smirks_matches(oemol, smirks):
        """Find all sets of atoms in the provided OpenEye molecule that match the provided SMARTS string.

        .. warning :: This API experimental and subject to change.

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
        # TODO: The MoleculeImage mapping should preserve ordering of template molecule for equivalent atoms
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

    @property
    def name(self):
        """The name (or title) of the molecule
        """
        return _name

    @staticmethod
    @requires_openeye('oechem', 'oeiupac')
    def from_iupac(iupac_name):
        """Generate a molecule from IUPAC or common name

        Parameters
        ----------
        iupac_name : str
            IUPAC name of molecule to be generated

        Returns
        -------
        molecule : Molecule
            The resulting molecule with position

        .. note :: This method requires the OpenEye toolkit to be installed.

        Examples
        --------

        Create a molecule from a common name

        >>> molecule = Molecule.from_iupac('4-[(4-methylpiperazin-1-yl)methyl]-N-(4-methyl-3-{[4-(pyridin-3-yl)pyrimidin-2-yl]amino}phenyl)benzamide')

        Create a molecule from a common name

        >>> molecule = Molecule.from_iupac('imatinib')

        """
        from openeye import oechem, oeiupac
        oemol = oechem.OEMol()
        oeiupac.OEParseIUPACName(oemol, iupac_name)
        oechem.OETriposAtomNames(oemol)
        return Molecule.from_openeye(oemol)

    @requires_openeye('oechem', 'oeiupac')
    def to_iupac(self):
        """Generate IUPAC name from Molecule

        Returns
        -------
        iupac_name : str
            IUPAC name of the molecule

        .. note :: This method requires the OpenEye toolkit to be installed.

        Examples
        --------

        >>> iupac_name = molecule.to_iupac()

        """
        from openeye import oeiupac
        return oeiupac.OECreateIUPACName(self.to_openeye())

    @staticmethod
    def from_topology(topology):
        """Return a Molecule representation of an openforcefield Topology containing a single Molecule object.

        Parameters
        ----------
        topology : openforcefield.topology.Topology
            The :class:`Topology` object containing a single :class:`Molecule` object.
            Note that OpenMM and MDTraj ``Topology`` objects are not supported.

        Returns
        -------
        molecule : openforcefield.topology.Molecule
            The Molecule object in the topology

        Raises
        ------
        ValueError
            If the topology does not contain exactly one molecule.

        Examples
        --------

        Create a molecule from a Topology object that contains exactly one molecule

        >>> molecule = Molecule.from_topology(topology)

        """
        # TODO: Ensure we are dealing with an openforcefield Topology object
        if topology.n_molecules != 1:
            raise ValueError('Topology must contain exactly one molecule')
        molecule = topology.unique_molecules.next()
        return Molecule(molecule)

    def to_topology(self):
        """
        Return an openforcefield Topology representation containing one copy of this molecule

        Returns
        -------
        topology : openforcefield.topology.Topology
            A Topology representation of this molecule

        Examples
        --------

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> topology = molecule.to_topology()

        """
        return Topology.from_molecules(self)

    @staticmethod
    def from_file(filename):
        """
        Create one or more molecules from a file

        .. todo::

           * Extend this to also include some form of .offmol Open Force Field Molecule format?
           * Generalize this to also include file-like objects?

        Parameters
        ----------
        filename : str
            The name of the file to stream one or more molecules from.

        Returns
        -------
        molecules : Molecule or list of Molecules
            If there is a single molecule in the file, a Molecule is returned;
            otherwise, a list of Molecule objects is returned.

        Examples
        --------
        >>> from openforcefield.tests.utils.utils import get_monomer_mol2file
        >>> mol2_filename = get_monomer_mol2file('cyclohexane')
        >>> molecule = Molecule.from_file(mol2_filename)

        """
        # Use highest-precendence toolkit
        toolkit = TOOLKIT_PRECEDENCE[0]
        mols = list()
        if toolkit == 'openeye':
            # Read molecules from an OpenEye-supported file, converting them one by one
            from openeye import oechem
            oemol = oechem.OEGraphMol()
            ifs = oechem.oemolistream(filename)
            while oechem.OEReadMolecule(ifs, oemol):
                mol = Molecule.from_openeye(oemol)
                mols.append(mol)
        elif toolkit == 'rdkit':
            from rdkit import Chem
            for rdmol in Chem.SupplierFromFilename(filename):
                mol = Molecule.from_rdkit(rdmol)
                mols.append(mol)
        else:
            raise Exception('Toolkit {} unsupported.'.format(toolkit))
        return mols

    def to_file(self, outfile, format):
        """Write the current molecule to a file or file-like object

        Parameters
        ----------
        outfile : str or file-like object
            A file-like object or the filename of the file to be written to
        format : str
            Format specifier, one of ['MOL2', 'MOL2H', 'SDF', 'PDB', 'SMI', 'CAN', 'TDT']
            Note that not all toolkits support all formats

        Raises
        ------
        ValueError
            If the requested format is not supported by one of the installed cheminformatics toolkits

        Examples
        --------

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> molecule.to_file('imatinib.mol2', format='mol2')
        >>> molecule.to_file('imatinib.sdf', format='sdf')
        >>> molecule.to_file('imatinib.pdb', format='pdb')

        """
        # Determine which formats are supported
        toolkit = None
        for query_toolkit in TOOLKIT_PRECEDENCE:
            if format in SUPPORTED_FILE_FORMATS[query_toolkit]:
                toolkit = query_toolkit
                break

        # Raise an exception if no toolkit was found to provide the requested format
        if toolkit == None:
            supported_formats = set()
            for toolkit in TOOLKIT_PRECEDENCE:
                supported_formats.add(SUPPORTED_FILE_FORMATS[toolkit])
            raise ValueError('The requested file format ({}) is not available from any of the installed toolkits (supported formats: {})'.format(format, supported_formats))

        # Write file
        if type(outfile) == str:
            # Open file for writing
            outfile = open(outfile, 'w')
            close_file_on_return = True
        else:
            close_file_on_return = False

        if toolkit == 'openeye':
            from openeye import oechem
            oemol = self.to_openeye()
            ofs = oechem.oemolostream(outfile)
            openeye_formats = getattr(oechem, 'OEFormat_' + format)
            ofs.SetFormat(openeye_formats[format])
            oechem.OEWriteMolecule(ofs, mol)
            ofs.close()
        elif toolkit == 'rdkit':
            from rdkit import Chem
            rdmol = self.to_rdkit()
            rdkit_writers = { 'SDF' : Chem.SDWriter, 'PDB' : Chem.PDBWriter, 'SMI' : Chem.SmilesWriter, 'TDT' : Chem.TDTWriter }
            writer = rdkit_writers[format](outfile)
            writer.write(rdmol)
            writer.close()

        if close_file_on_return:
            outfile.close()

    @staticmethod
    @requires_rdkit()
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

        Examples
        --------

        Create a molecule from an RDKit molecule

        >>> molecule = Molecule.from_rdkit(rdmol)

        """
        from rdkit import Chem

        # Create a new openforcefield Molecule
        mol = Molecule()

        # If RDMol has a title save it
        if rdmol.HasProp("_Name"):
            self.name == rdmol.GetProp("_Name")

        # Store all properties
        # TODO: Should Title or _Name be a special property?
        # TODO: Should there be an API point for storing properties?
        properties = rdmol.GetPropsAsDict()
        mol.properties = properties

        # We store bond orders as integers regardless of aromaticity.
        # In order to properly extract these, we need to have the "Kekulized" version of the rdkit mol
        kekul_mol = Chem.Mol(rdmol)
        Chem.Kekulize(kekul_mol, True)

        # setting chirality in openeye requires using neighbor atoms
        # therefore we can't do it until after the atoms and bonds are all added
        chiral_atoms = dict() # {rd_idx: openeye chirality}
        for rda in rdmol.GetAtoms():
            rd_idx = rda.GetIdx()

            # create a new atom
            atomic_number = oemol.NewAtom(rda.GetAtomicNum())
            formal_charge = rda.GetFormalCharge()
            is_aromatic = rda.GetIsAromatic()

            # If chiral, store the chirality to be set later
            stereochemistry = None
            tag = rda.GetChiralTag()
            if tag == Chem.CHI_TETRAHEDRAL_CCW:
                stereochemistry = 'R'
            if tag == Chem.CHI_TETRAHEDRAL_CW:
                stereochemistry = 'S'

            atom_index = mol.add_atom(atomic_number=atomic_number, formal_charge=formal_charge, is_aromatic=is_aromatic, stereochemistry=stereochemistry)
            map_atoms[rd_idx] = atom_index

        # Similar to chirality, stereochemistry of bonds in OE is set relative to their neighbors
        stereo_bonds = list()
        # stereo_bonds stores tuples in the form (oe_bond, rd_idx1, rd_idx2, OE stereo specification)
        # where rd_idx1 and 2 are the atoms on the outside of the bond
        # i.e. Cl and F in the example above
        aro_bond = 0
        for rdb in rdmol.GetBonds():
            a1 = rdb.GetBeginAtomIdx()
            a2 = rdb.GetEndAtomIdx()

            # Determine bond aromaticity and Kekulized bond order
            is_aromatic = False
            order = rdb.GetBondTypeAsDouble()
            if order == 1.5:
                # get the bond order for this bond in the kekulized molecule
                order = kekul_mol.GetBondWithIdx(rdb.GetIdx()).GetBondTypeAsDouble()
                is_aromatic = True
            # Convert floating-point bond order to integral bond order
            order = int(order)

            # determine if stereochemistry is needed
            stereochemistry = None
            tag = rdb.GetStereo()
            if tag == Chem.BondStereo.STEREOCIS or tag == Chem.BondStereo.STEREOZ:
                stereochemistry = 'Z'
            if tag == Chem.BondStereo.STEREOTRANS or tag == Chem.BondStereo.STEREOE:
                stereochemistry = 'E'

            # create a new bond
            bond_index = mol.add_bond(map_atoms[a1], map_atoms[a2], is_aromatic=is_aromatic, order=order, stereochemistry=stereochemistry)

        # TODO: Save conformer(s), if present
        # If the rdmol has a conformer, store its coordinates
        # TODO: Note, this currently only adds the first conformer, it will need to be adjusted if the
        # you wanted to convert multiple sets of coordinates
        if rdmol.GetConformers():
            conf = rdmol.GetConformer()
            # TODO: Store conformers
            #for rd_idx, oeatom in map_atoms.items():
            #    coords = conf.GetAtomPosition(rd_idx)
            #    oemol.SetCoords(oeatom, oechem.OEFloatArray(coords))

        return mol

    @requires_rdkit()
    def to_rdkit(self, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an RDKit molecule

        Requires the RDKit to be installed.

        Parameters
        ----------
        aromaticity_model : str, optional, default=DEFAULT_AROMATICITY_MODEL
            The aromaticity model to use

        Returns
        -------
        rdmol : rkit.RDMol
            An RDKit molecule

        Examples
        --------

        Convert a molecule to RDKit

        >>> rdmol = molecule.to_rdkit()

        """
        from rdkit import Chem, Geometry

        # Create an editable RDKit molecule
        rdmol = Chem.RWMol()

        # Set name
        # TODO: What is the best practice for how this should be named?
        rdmol.SetProp('_Name', self.name)

        # TODO: Set other properties
        for name, value in self.properties.items():
            if type(value) == str:
                rdmol.SetProp(name, value)
            elif type(value) == int:
                rdmol.SetIntProp(name, value)
            elif type(value) == float:
                rdmol.SetDoubleProp(name, value)
            elif type(value) == bool:
                rdmol.SetBoolProp(name, value)
            else:
                # Shove everything else into a string
                rdmol.SetProp(name, str(value))

        _bondtypes = {1: Chem.BondType.SINGLE,
                      1.5: Chem.BondType.AROMATIC,
                      2: Chem.BondType.DOUBLE,
                      3: Chem.BondType.TRIPLE,
                      4: Chem.BondType.QUADRUPLE,
                      5: Chem.BondType.QUINTUPLE,
                      6: Chem.BondType.HEXTUPLE,
                      7: Chem.BondType.ONEANDAHALF,}

        # atom map lets you find atoms again
        map_atoms = dict() # { molecule index : rdkit index }
        for index, atom in enumerate(self.atoms):
            rdatom = Chem.Atom(atom.atomic_number)
            rdatom.SetFormalCharge(atom.formal_charge)
            rdatom.SetIsAromatic(atom.is_aromatic)

            if atom.stereochemistry == 'S':
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)
            elif atom.stereochemistry == 'R':
                rdatom.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)

            map_atoms[oe_idx] = rdmol.AddAtom(rdatom)

        for bond in self.bonds:
            rdatom1 = map_atoms[bond.atom1_index]
            rdatom2 = map_atoms[bond.atom2_index]
            rdmol.AddBond(rdatom1, rdatom2)
            rdbond = rdmol.GetBondBetweenAtoms(rdatom1, rdatom2)

            # Assign bond type, which is based on order unless it is aromatic
            if bond.is_aromatic:
                rdbond.SetBondType(_bondtypes[1.5])
                rdbond.SetIsAromatic(True)
            else:
                rdbond.SetBondType(_bondtypes[bond.order])
                rdbond.SetIsAromatic(False)

        # Assign bond stereochemistry
        for bond in self.bonds:
            if bond.stereochemistry:
                # Determine neighbors
                # TODO: This API needs to be created
                n1 = [n.index for n in bond.atom1.bonded_atoms if n != bond.atom2][0]
                n2 = [n.index for n in bond.atom2.bonded_atoms if n != bond.atom1][0]
                # Get rdmol bonds
                bond1 = rdmol.GetBondBetweenAtoms(map_atoms[n1], map_atoms[bond.atom1_index])
                bond2 = rdmol.GetBondBetweenAtoms(map_atoms[bond.atom1_index], map_atoms[bond.atom2_index])
                bond3 = rdmol.GetBondBetweenAtoms(map_atoms[bond.atom2_index], map_atoms[n2])
                # Set arbitrary stereochemistry
                # Since this is relative, the first bond always goes up
                # as explained above these names come from SMILES slashes so UP/UP is Trans and Up/Down is cis
                bond1.SetBondDir(Chem.BondDir.ENDUPRIGHT)
                bond3.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
                # Flip the stereochemistry if it is incorrect
                # TODO: Clean up _CIPCode atom and bond properties
                Chem.AssignStereochemistry(rdmol, cleanIt=True, force=True)
                if rdmol.GetProp('_CIPCode') != bond.stereochemistry:
                    # Flip it
                    bond3.SetBondDir(Chem.BondDir.ENDUPRIGHT)
                    # Validate we have the right stereochemistry as a sanity check
                    Chem.AssignStereochemistry(rdmol, cleanIt=True, force=True)
                    if rdmol.GetProp('_CIPCode') != bond.stereochemistry:
                        raise Exception('Programming error with assumptions about RDKit stereochemistry model')

        # Set coordinates if we have them
        # TODO: Fix this once conformer API is defined
        if self._conformers:
            for conformer in self._conformers:
                rdmol_conformer = Chem.Conformer()
                for index, rd_idx in map_atoms.items():
                    (x,y,z) = conformer[index,:]
                    rdmol_conformer.SetAtomPosition(rd_idx, Geometry.Point3D(x,y,z))
                rdmol.AddConformer(rdmol_conformer)

        # Cleanup the rdmol
        # Note I copied UpdatePropertyCache and GetSSSR from Shuzhe's code to convert oemol to rdmol here:
        rdmol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(rdmol)
        # I added AssignStereochemistry which takes the directions of the bond set
        # and assigns the stereochemistry tags on the double bonds
        Chem.AssignStereochemistry(rdmol, force=False)

        # Return non-editable version
        return rdkit.Mol(rdmol)

    @staticmethod
    @requires_openeye('oechem')
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
        molecule : openforcefield.topology.Molecule
            An openforcefield molecule

        Examples
        --------

        Create a Molecule from an OpenEye OEMol

        >>> molecule = Molecule.from_openeye(oemol)

        """
        return toolkit_registry.call('from_openeye', oemol)

    def to_openeye(self, aromaticity_model=DEFAULT_AROMATICITY_MODEL):
        """
        Create an OpenEye molecule

        Requires the OpenEye toolkit to be installed.

        .. todo ::

           * Use stored conformer positions instead of an argument.
           * Should the aromaticity model be specified in some other way?

        Parameters
        ----------
        aromaticity_model : str, optional, default=DEFAULT_AROMATICITY_MODEL
            The aromaticity model to use

        Returns
        -------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule

        Examples
        --------

        Create an OpenEye molecule from a Molecule

        >>> molecule = Molecule.from_smiles('CC')
        >>> oemol = molecule.to_openeye()

        """
        return toolkit_registry.call('to_openeye', self, aromaticity_model=aromaticity_model)

    # TODO: We have to distinguish between retrieving user-specified partial charges and providing a generalized semiempirical/pop analysis/BCC scheme according to the new SMIRNOFF spec
    def get_partial_charges(self, method='AM1-BCC'):
        """
        Retrieve partial atomic charges.

        .. warning :: This API experimental and subject to change.

        .. todo::
            * Generalize to allow specification of both QM method and population analysis method
            * Should we return the charges or store the charges as atom properties?
            * Refine API for this method to better correspond to new SMIRNOFF 1.0 spec
            * Is it OK that the ``Molecule`` object does not store geometry, but will create it using ``openeye.omega`` or ``rdkit``?
            * Should this method assign charges to the ``Atom``s in the molecule, a separate ``charges`` molecule property,
              or just return the charge array? Is it OK that the Topology is modified?
            * How do we add enough flexibility to specify the toolkit and optional parameters, such as:
              ``oequacpac.OEAssignPartialCharges(charged_copy, getattr(oequacpac, 'OECharges_AM1BCCSym'), False, False)``
            * What will we do about virtual sites, since this only assigns partial atomic charges?

        Parameters
        ----------
        method : str, optional, default='AM1-BCC'
            The name of the charge method to use.
            Options are:
            * `AM1-BCC` : symmetrized AM1 charges with BCC (no ELF)

        Returns
        -------
        charges : numpy.array of shape (natoms) of type float
            The partial charges

        Examples
        --------

        Get AM1-BCC charges for imatinib

        >>> molecule = Molecule.from_iupac('imatininb')
        >>> charges = molecule.get_partial_charges(method='AM1-BCC')

        """
        # TODO: Use memoization to speed up subsequent calls; use decorator?
        charges = toolkit_registry.call('compute_partial_charges', method=method)
        return charges

    def get_fractional_bond_orders(self, method='Wiberg'):
        """Get fractional bond orders.

        .. warning :: This API experimental and subject to change.

        .. todo::
            * Is it OK that the ``Molecule`` object does not store geometry, but will create it using ``openeye.omega`` or ``rdkit``?
            * Should this method assign fractional bond orders to the ``Bond``s in the molecule, a separate ``bond_orders`` molecule property,
              or just return the array of bond orders?
            * How do we add enough flexibility to specify the toolkit and optional parameters, such as:
              ``oequacpac.OEAssignPartialCharges(charged_copy, getattr(oequacpac, 'OECharges_AM1BCCSym'), False, False)``
            * Generalize to allow user to specify both QM method and bond order computation approach (e.g. ``AM1`` and ``Wiberg``)

        method : str, optional, default='Wiberg'
            The name of the charge method to use.
            Options are:
            * 'Wiberg' : Wiberg bond order

        Examples
        --------

        Get fractional Wiberg bond orders

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> fractional_bond_orders = molecule.get_fractional_bond_orders(method='Wiberg')

        """
        # TODO: Use memoization to speed up subsequent calls; use decorator?
        fractional_bond_orders = toolkit_registry.call('compute_fractional_bond_orders', method=method)
        return fractional_bond_orders

    # TODO: Compute terms for each unique molecule, then use mapping to molecules to enumerate all terms
    @property
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
    @property
    def torsions(self):
        """
        Get an iterator over all i-j-k-l torsions.
        Note that i-j-k-i torsions (cycles) are excluded.

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
