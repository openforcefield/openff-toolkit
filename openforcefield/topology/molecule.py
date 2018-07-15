#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Molecular chemical entity representation and routines to interface with cheminformatics toolkits

.. todo::

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
from simtk.openmm.app import element as elem

import openforcefield
from openforcefield.utils.toolkits import OPENEYE_INSTALLED, RDKIT_INSTALLED, TOOLKIT_PRECEDENCE, SUPPORTED_FILE_FORMATS
from openforcefield.utils.toolkits import requires_rdkit, requires_openeye
from openforcefield.typing.chemistry import ChemicalEnvironment, SMIRKSParsingError

#=============================================================================================
# GLOBAL PARAMETERS
#=============================================================================================

# TODO: Move these to utils.toolkits?

# TODO: Can we have the `ALLOWED_*_MODELS` list automatically appear in the docstrings below?
# TODO: Should `ALLOWED_*_MODELS` be objects instead of strings?
# TODO: Should these be imported from `openforcefield.cheminformatics.aromaticity_models` and `.bondorder_models`?

# TODO: Allow all OpenEye aromaticity models to be used with OpenEye names?
#       Only support OEAroModel_MDL in RDKit version?

DEFAULT_AROMATICITY_MODEL = 'OEAroModel_MDL' # TODO: Is there a more specific name and reference for the aromaticity model?
ALLOWED_AROMATICITY_MODELS = ['OEAroModel_MDL']

DEFAULT_FRACTIONAL_BONDORDER_MODEL = 'Wiberg' # TODO: Is there a more specific name and reference for the aromatciity model?
ALLOWED_FRACTIONAL_BONDORDER_MODELS = ['Wiberg']

DEFAULT_CHARGE_MODEL = 'AM1-BCC' # TODO: Should this be `AM1-BCC`, or should we encode BCCs explicitly via AM1-CM2 preprocessing?
ALLOWED_CHARGE_MODELS = ['AM1-BCC'] # TODO: Which models do we want to support?

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

#=============================================================================================
# Particle
#=============================================================================================

class Particle(object):
    """
    Base class for all particles in a molecule.

    .. warning :: This API experimental and subject to change.

    A particle object could be an ``Atom`` or a ``VirtualSite``.

    """
    def __init__(self):
        """
        Create a particle.
        """
        # TODO: Don't allow a Particle to be constructed directly; only an Atom or VirtualSite
        raise Exception('A Particle cannot be created directly; only an Atom or VirtualSite can be created directly')
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
        Create an Atom object.

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

        .. todo :: Use attrs to validate?
        """
        super(Atom, self).__init__()


        self._atomic_number = atomic_number
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

        .. todo ::

           * This will need to be gneeralized for virtual sites to allow out-of-plane sites, which are not simply a linear combination of atomic positions
           * Add checking for allowed virtual site types
           * How do we want users to specify virtual site types?

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
        # TODO: Also include particle_index, which molecule this atom belongs to?
        return "VirtualSite(name={}, type={}, weights={}, atoms={})".format(self.name, self.type, self.weights, self.atoms)

    def __str__(self):
        # TODO: Also include particle_index, which molecule this atom belongs to?
        return "<VirtualSite name={} type={} weights={}, atoms={}>".format(self.name, self.type, self.weights, self.atoms)

#=============================================================================================
# Bond
#=============================================================================================

class Bond(object):
    """
    Chemical bond representation.

    .. warning :: This API experimental and subject to change.

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

    .. todo :: Allow bonds to have associated properties.

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
# Molecule
#=============================================================================================

# TODO: Make Molecule immutable (by default)

class Molecule(object):
    """
    Chemical representation of a molecule.

    Examples
    --------

    Create a molecule from a mol2 file

    >>> molecule = Molecule.from_file('molecule.mol2')

    Create a molecule from an OpenEye molecule

    >>> molecule = Molecule.from_openeye(oemol)

    Create a molecule from an RDKit molecule

    >>> molecule = Molecule.from_rdkit(rdmol)

    Create a molecule from IUPAC name (OpenEye toolkit required)

    >>> molecule = Molecule.from_iupac('imatinib')

    .. todo ::

       * How do we automatically trigger invalidation of cached properties if an ``Atom``, ``Bond``, or ``VirtualSite`` is modified,
         rather than added/deleted via the API? The simplest resolution is simply to make them immutable.

    """
    def __init__(self, other=None):
        """
        Create a new Molecule object

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

        .. todo ::

           * If a filename or file-like object is specified but the file contains more than one molecule, what is the proper behavior?
           Read just the first molecule, or raise an exception if more than one molecule is found?

           * Should we also support SMILES strings for ``other``?

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

    def _initialize(self):
        """
        Clear the contents of the current molecule.
        """
        self._particles = list()
        self._bonds = None
        self._name = None # Set the name of the molecule
        self._charges = None # TODO: Storage charges
        self._properties = None # Attached properties
        self._cached_properties = None # Cached properties can be recomputed as needed
        self._conformers = None

    def _copy_initializer(self, other):
        """
        Copy contents of the specified molecule

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

        .. todo :: Should this be a ``@staticmethod`` where we have an explicit copy constructor?

        """
        assert isinstance(other, type(self)), "can only copy instances of {}".format(typse(self))
        self.__dict__ = deepcopy(other.__dict__)

    def _eq_(self, other):
        """Test two molecules for equality.

        .. note ::

           Note that this method simply tests whether two molecules are identical chemical species.
           No effort is made to ensure that the atoms are in the same order or that any annotated properties are preserved.

        """
        return self.to_smiles() == other.to_smiles()

    def to_smiles(self):
        """
        Return a canonical isomeric SMILES representation of the current molecule

        .. note :: RDKit and OpenEye versions will not necessarily return the same representation.

        Returns
        -------
        smiles : str
            Canonical isomeric explicit-hydrogen SMILES

        .. todo :: Can we ensure RDKit and OpenEye versions return the same representation?

        """
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

        Returns
        -------
        graph : networkx.Graph
            The resulting graph, with nodes labeled with atom indices and elements

        .. todo ::

           * Do we need a ``from_networkx()`` method? If so, what would the Graph be required to provide?
           * Should edges be labeled with discrete bond types in some aromaticity model?
           * Should edges be labeled with fractional bond order if a method is specified?
           * Should we add other per-atom and per-bond properties (e.g. partial charges) if present?

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

        Examples
        --------
        Retrieve all the carbon-carbon bond matches in a molecule

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> matches = molecule.chemical_environment_matches('[#6X3:1]~[#6X3:2]')

        """
        # Resolve to SMIRKS if needed
        # TODO: Udpdate
        if hasattr(query, 'asSMIRKS'):
            smirks = query.asSMIRKS()
        else:
            smirks = query

        # Use specified cheminformatics toolkit to determine matches with specified aromaticity model
        toolkit = TOOLKIT_PRECEDENCE[0]
        if toolkit == 'oechem':
            matches = _openye_smirks_matches(self.to_openeye(), smirks, aromaticity_model=self._aromaticity_model)
        elif toolkit == 'rdkit':
            matches = _openeye_smirks_matches(self.to_rdkit(), smirks, aromaticity_model=self._aromaticity_model)
        else:
            raise Exception('Unknown toolkit {}'.format(toolkit))

        return matches

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

    # TODO: Implement specialized JSON/BSON deserialization, if needed.
    #def __setstate__(self, state):
    #    pass

    # TODO: Implement specialized JSON/BSON serialization, if needed
    #def __getstate__(self):
    #    pass

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

        .. todo::

           * Extend this to also include some form of .offmol Open Force Field Molecule format?
           * Generalize this to also include file-like objects?

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

        Returns
        -------
        rdmol : rkit.RDMol
            An RDKit molecule
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

        """
        from openeye import oechem
        from openforcefield.utils.toolkits import openeye_cip_atom_stereochemistry, openeye_cip_bond_stereochemistry

        molecule = Molecule()

        # TODO: What other information should we preserve besides name?
        # TODO: How should we preserve the name?
        molecule.name = oemol.GetTitle()

        # Copy any attached SD tag information
        # TODO: Should we use an API for this?
        molecule.properties = dict()
        for dp in oechem.OEGetSDDataPairs(oemol):
            molecule[dp.GetTag()] = dp.GetValue()

        map_atoms = dict() # {oemol_idx: molecule_idx}
        for oeatom in oemol.GetAtoms():
            oe_idx = oeatom.GetIdx()
            atomic_number = oeatom.GetAtomicNum()
            formal_charge = oeatom.GetFormalCharge()
            is_aromatic = oeatom.IsAromatic()
            stereochemistry = openeye_cip_atom_stereochemistry(oemol, oeatom)
            atom_index = molecule.add_atom(atomic_number=atomic_number, formal_charge=formal_charge, is_aromatic=is_aromatic, stereochemistry=stereochemistry)
            map_atoms[oe_idx] = atom_index # store for mapping oeatom to molecule atom indices below

        for oebond in oemol.GetBonds():
            atom1_index = map_atoms[oebond.GetBgnIdx()]
            atom2_index = map_atoms[oebond.GetEndIdx()]
            order = oeb.GetOrder()
            is_aromatic = oeb.IsAromatic()
            stereochemistry = openeye_cip_bond_stereochemistry(oemol, oebond)
            molecule.add_bond(atom1_index, atom2_index, order=order, is_aromatic=is_aromatic, stereochemistry=stereochemistry)

        # TODO: Copy conformations, if present

        return molecule

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
        from openeye import oechem
        from openforcefield.utils.toolkits import openeye_cip_atom_stereochemistry, openeye_cip_bond_stereochemistry

        oemol = oechem.OEMol()

        # Add atoms
        oemol_atoms = list() # list of corresponding oemol atoms
        for atom in self.atoms:
            oeatom = oemol.NewAtom(atom.atomic_number)
            oeatom.SetFormalCharge(atom.formal_charge)
            oeatom.SetAromatic(atom.is_aromatic)
            oemol_atoms.append(oeatom)

        # Add bonds
        oemol_bonds = list() # list of corresponding oemol bonds
        for bond in self.bonds:
            oebond = oemol.NewBond(oemol_atoms[bond.atom1_index], oemol_atoms[bond.atom2_index])
            newbond.SetOrder(bond.order)
            newbond.SetAromatic(bond.is_aromatic)
            oemol_bonds.append(oebond)

        # Set atom stereochemistry now that all connectivity is in place
        for atom, oeatom in zip(self.atoms, oemol_atoms):
            if not atom.stereochemistry:
                continue

            # Set arbitrary initial stereochemistry
            neighs = [n for n in oeatom.GetAtoms()]
            oeatom.SetStereo(neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Right)

            # Flip chirality if stereochemistry is incorrect
            oeatom_stereochemistry = openeye_cip_atom_stereochemistry(oemol, oeatom)
            if oeatom_stereochemistry != atom.sterechemistry:
                # Flip the stereochemistry
                oea.SetStereo(neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Left)
                # Verify it matches now as a sanity check
                oeatom_stereochemistry = openeye_cip_atom_stereochemistry(oemol, oeatom)
                if oeatom_stereochemistry != atom.stereochemistry:
                    raise Exception('Programming error: OpenEye atom stereochemistry assumptions failed.')

        # Set bond stereochemistry
        for bond, oebond in zip(self.atoms, oemol_bonds):
            if not bond.stereochemistry:
                continue

            # Set arbitrary initial stereochemistry
            oeatom1, oeatom2 = oemol_atoms[bond.atom1_index], oemol_atoms[bond.atom2_index]
            oeatom1_neighbor = [n for n in oeatom1.GetAtoms()][0]
            oeatom2_neighbor = [n for n in oeatom2.GetAtoms()][0]
            oebond.SetStereo([oeatom1, oeatom2], oechem.OEBondStereo_CisTrans, oechem.OEBondStereo_Cis)

            # Flip stereochemistry if incorrect
            oebond_stereochemistry = openeye_cip_bond_stereochemistry(oemol, oebond)
            if oebond_stereochemistry != bond.sterechemistry:
                # Flip the stereochemistry
                oebond.SetStereo([oeatom1, oeatom2], oechem.OEBondStereo_CisTrans, oechem.OEBondStereo_Trans)
                # Verify it matches now as a sanity check
                oebond_stereochemistry = openeye_cip_bond_stereochemistry(oemol, oebond)
                if oebond_stereochemistry != bond.stereochemistry:
                    raise Exception('Programming error: OpenEye bond stereochemistry assumptions failed.')

        # TODO: Save conformations, if present

        # TODO: Save name and properties, if present

        # Clean Up phase
        # The only feature of a molecule that wasn't perceived above seemed to be ring connectivity, better to run it
        # here then for someone to inquire about ring sizes and get 0 when it shouldn't be
        oechem.OEFindRingAtomsAndBonds(oemol)

        return oemol

    # TODO: We have to distinguish between retrieving user-specified partial charges and providing a generalized semiempirical/pop analysis/BCC scheme according to the new SMIRNOFF spec
    def get_partial_charges(self, method='AM1-BCC'):
        """
        Retrieve partial atomic charges.

        .. warning :: This API experimental and subject to change.

        .. todo::
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

        """
        # TODO: Cache charges to speed things up

        # TODO: Support methods beyond AM1-BCC
        if method != 'AM1-BCC':
            raise NotImplementedError()

        if OPENEYE_INSTALLED:
            # TODO: Translate charge model and check compatibility
            charge_model = 'am1bcc'
            charges = _assign_partial_charges_using_quacpac(charge_model)
        elif RDKIT_INSTALLED:
            # TODO: Translate charge model and check compatibility
            charge_model = 'bcc'
            charges = _assign_partial_charges_using_antechamber(charge_model)

        # TODO: Add units before returning charges, either here or in private methods below
        return charges

    def _assign_partial_charges_using_antechamber(self, charge_model="bcc"):
        """
        Assign partial charges with AmberTools using antechamber/sqm

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        charge_model : str, optional, default='bcc'
            The charge model to use. One of ['gas', 'mul', 'cm1', 'cm2', 'bcc']

        Returns
        -------
        charges : numpy.array of shape (natoms) of type float
            The partial charges

        Notes
        -----
        Currently only sdf file supported as input and mol2 as output
        https://github.com/choderalab/openmoltools/blob/master/openmoltools/packmol.py

        .. todo ::

           * Do we want to also allow ESP/RESP charges?

        """
        # Check that the requested charge method is supported
        SUPPORTED_ANTECHAMBER_CHARGE_METHODS = ['gas', 'mul', 'cm1', 'cm2', 'bcc']
        if charge_method not in SUPPORTED_ANTECHAMBER_CHARGE_METHODS:
            raise ValueError('Requested charge method {} not among supported charge methods {}'.format(charge_method, SUPPORTED_ANTECHAMBER_CHARGE_METHODS))

        # Find the path to antechamber
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            raise(IOError("Antechamber not found, cannot run charge_mol()"))

        # Compute charges
        from openmmtools.utils import temporary_directory, temporary_cd
        with temporary_directory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = self.net_charge()
                # Write out molecule in SDF format
                self.to_file('molecule.sdf', format='SDF')
                # Compute desired charges
                # TODO: Add error handling if antechamber chokes
                os.system("antechamber -i molecule.sdf -fi sdf -o charged.mol2 -fo mol2 -pf yes -c {} -nc {}".format(charge_method, net_charge))
                # Write out just charges
                os.system("antechamber -i charges.mol2 -fi mol2 -o charges2.mol2 -fo mol2 -c wc -cf charges.txt -pf yes")
                # Read the charges
                with open('charges.txt', 'r') as infile:
                    contents = infile.read()
                text_charges = contents.split()
                charges = np.zeros([self.n_atoms], np.float64)
                for index, token in enumerate(text_charges):
                    charges[index] = float(token)

        return charges

    def _assign_partial_charges_using_quacpac(molecule, charge_model="am1bcc"):
        """
        Assign partial charges with OpenEye quacpac

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        molecule : Molecule
            Molecule for which partial charges are to be computed
        charge_model : str, optional, default='bcc'
            The charge model to use. One of ['noop', 'mmff', 'mmff94', 'am1bcc', 'am1bccnosymspt', 'amber', 'amberff94', 'am1bccelf10']

        Returns
        -------
        charges : numpy.array of shape (natoms) of type float
            The partial charges

        .. todo ::

           * Should the default be ELF?
           * Can we expose more charge models?

        """
        oemol = molecule.to_openeye()

        result = False
        if name == "noop":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEChargeEngineNoOp())
        elif name == "mmff" or name == "mmff94":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEMMFF94Charges())
        elif name == "am1bcc":
            result = oequacpac.OEAssignCharges(oemol, oequacpac.OEAM1BCCCharges())
        elif name == "am1bccnosymspt":
            optimize = True
            symmetrize = True
            result = oequacpac.OEAssignCharges(mol, oequacpac.OEAM1BCCCharges(not optimize, not symmetrize))
        elif name == "amber" or name == "amberff94":
            result = oequacpac.OEAssignCharges(mol, oequacpac.OEAmberFF94Charges())
        elif name == "am1bccelf10":
            result = oequacpac.OEAssignCharges(mol, oequacpac.OEAM1BCCELF10Charges())
        else:
            raise ValueError('charge_model {} unknown'.format(charge_model))

        if result is False:
            raise Exception('Unable to assign charges')

        # Extract and return charges
        charges = np.zeros([oemol.NumAtoms()], np.float64)
        for index, atom in enumerate(oemol.GetAtoms()):
            charges[index] = atom.GetPartialCharge()
        return charges

    def get_fractional_bond_orders(self, method='Wiberg', toolkit=None, **kwargs):
        """Get fractional bond orders.

        .. warning :: This API experimental and subject to change.

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
