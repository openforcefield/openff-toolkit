#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Molecular chemical entity representation and routines to interface with cheminformatics toolkits

.. todo::

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

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import numpy as np
from copy import deepcopy

from distutils.spawn import find_executable

from simtk import unit
from simtk.openmm.app import element as elem

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

DEFAULT_CHARGE_MODEL = 'AM1-CM2' # TODO: Should this be `AM1-BCC`, or should we encode BCCs explicitly?
ALLOWED_CHARGE_MODELS = ['AM1-CM2', 'AM1-BCC', 'Mulliken'] # TODO: Which models do we want to support?

#=============================================================================================
# PRIVATE SUBROUTINES
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
        # TODO: Don't allow a Particle to be constructed directly; only an Atom or VirtualSite

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
        self._cached_properties = None # clear cached properties

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
                raise ValueError("Could not process SMILES string: {}".format(smiles))
            return Molecule.from_openeye(oemol)
        elif toolkit == 'rdkit':
            rdmol = Chem.MolFromSmiles(smiles)
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

    def add_atom(atomic_number, formal_charge, is_aromatic, stereochemistry=None):
        """
        Add an atom

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
        >>> molecule.add_bond(C, H1)
        >>> molecule.add_bond(C, H2)
        >>> molecule.add_bond(C, H3)
        >>> molecule.add_bond(C, H4)

        .. warning :: This API experimental and subject to change.

        """
        # Create an atom
        atom = Atom(atomic_number=atomic_number, formal_charge=formal_charge, is_aromatic=is_aromatic, stereochemistry=stereochemistry)
        self._atoms.append(atom)
        self._invalidate_cached_properties()

    # TODO: Should we allow the removal of atoms too?
    # TODO: Should invalidation of cached properties be handled via something like a tracked list?

    def add_bond(self, atom1_index, atom2_index, ):
        """
        Add a bond between two specified atom indices

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
            If None, stereochemistry is ambiguous; otherwise 'E' or 'Z'

        .. warning :: This API experimental and subject to change.

        """
        # TODO: Check to make sure bond does not already exist
        self._bonds.append(atom1, atom2)
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
        pass

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

        .. warning :: This API experimental and subject to change.

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

        .. warning :: This API experimental and subject to change.

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
        if toolkit == 'openeye':
            # Read molecules from an OpenEye-supported file, converting them one by one
            from openeye import oechem
            oemol = oechem.OEGraphMol()
            ifs = oechem.oemolistream(filename)
            while oechem.OEReadMolecule(ifs, mol):
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
        from rdkit import Chem

        # Create a new openforcefield Molecule
        rdmol = Chem.EditableMol()
        #
        # # Set name
        # # TODO: What is the best practice for how this should be named?
        # rdmol.SetProp('Name', self.name)
        #
        # # Add atoms
        # atom_map = dict() # atom_map[offmol_atom] is the rdkit atom index of openforcefield Atom offmol_atom
        # for atom in self.atoms:
        #     rdatom_index = rdmol.AddAtom(Chem.Atom(atom.element_index))
        #     rdatom = rdmol.GetAtomWithIdx(rdmol, rdatom_index)
        #     rdatom.SetProp('_TriposAtomName', atom.name) # TODO: Should we use a different property name to store the atom name?
        #     atom_map[atom] = rdatom_index
        #
        # # Add bonds
        # from rdkit.Chem import rdchem
        # # Mapping from RDKit BondType to openforcefield (type, int_order)
        # # TODO: This mapping may be woefully incomplete
        # BOND_TYPE_MAP = {
        #     rdchem.BondType.SINGLE : ('single', 1),
        #     rdchem.BondType.DOUBLE : ('double', 2),
        #     rdchem.BondType.TRIPLE : ('triple', 3),
        #     rdchem.BondType.AROMATIC : ('aromatic', 1),
        # }
        # for bond in self.bonds:
        #     if bond.type == 'single':
        #         rdkit_bond_type = rdchem.BondType.SINGLE
        #     elif bond.type == 'double':
        #         rdkit_bond_type = rdchem.BondType.DOUBLE
        #     elif bond.type == 'triple':
        #         rdkit_bond_type = rdchem.BondType.TRIPLE
        #     elif bond.type == 'aromatic':
        #         rdkit_bond_type = rdchem.BondType.AROMATIC
        #     else:
        #         raise ValueError('bond type {} unknown'.format(bond.type))
        #     rdmol.AddBond(atom_map[bond.atom1], atom_map[bond.atom2], rdkit_bond_type)
        #
        # # TODO: Preserve atom and bond stereochemistry

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
        for openeye_atom in oemol.GetAtoms():
            oe_idx = openeye_atom.GetIdx()
            atomic_number = oea.GetAtomicNum()
            formal_charge = oea.GetFormalCharge()
            is_aromatic = oea.IsAromatic()

            # Retrieve stereochemistry
            # TODO: The stereochemistry doesn't really make sense except in the context of other atoms; is there a better way to handle this?
            cip = oechem.OEPerceiveCIPStereo(oemol, oea)
            if cip == oechem.OECIPAtomStereo_S:
                stereochemistry = 'S'
            if cip == oechem.OECIPAtomStereo_R:
                stereochemistry = 'R'

            atom_index = molecule.add_atom(atomic_number=atomic_number, formal_charge=formal_charge, is_aromatic=is_aromatic, stereochemistry=stereochemistry)
            map_atoms[oe_idx] = atom_index

        aro_bond = 0
        for oeb in oemol.GetBonds():
            # get neighboring rd atoms
            atom1_index = map_atoms[oeb.GetBgnIdx()]
            atom2_index = map_atoms[oeb.GetEndIdx()]

            order = oeb.GetOrder()
            is_aromatic = oeb.IsAromatic()

            stereochemistry = None
            # If the bond has specified stereo add the required information to stereo_bonds
            # TODO: Convert this into E/Z
            if oeb.HasStereoSpecified(oechem.OEBondStereo_CisTrans):
                # OpenEye determines stereo based on neighboring atoms so get two outside atoms
                n1 = [n for n in oeb.GetBgn().GetAtoms() if n != oeb.GetEnd()][0]
                n2 = [n for n in oeb.GetEnd().GetAtoms() if n != oeb.GetBgn()][0]
                stereo = oeb.GetStereo([n1,n2], oechem.OEBondStereo_CisTrans)
                if stereo == oechem.OEBondStereo_Cis:
                    stereochemistry = 'cis'
                elif stereo == oechem.OEBondStereo_Trans:
                    stereochemistry = 'trans'

            molecule.add_bond(atom1_index, atom2_index, order=order, is_aromatic=is_aromatic, stereochemistry=stereochemistry)

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

            # Set default stereochemistry of OpenEye atom
            neighs = [n for n in oeatom.GetAtoms()]
            oeatom.SetStereo(neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Right)

            # Flip chirality if stereochemistry is incorrecft
            cip = oechem.OEPerceiveCIPStereo(oemol, oeatom)
            # TODO: Is there a chance the CIP stereo could be other than R and S in corner cases?
            oeatom_stereochemistry = 'R' if (cip == oechem.OECIPAtomStereo_R) else 'S'

            # Flip chirality if needed
            if oeatom_stereochemistry != atom.sterechemistry:
                # Flip the stereochemistry
                oea.SetStereo(neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Left)
                # Verify it matches now
                cip = oechem.OEPerceiveCIPStereo(oemol, oeatom)
                # TODO: Is there a chance the CIP stereo could be other than R and S in corner cases?
                oeatom_stereochemistry = 'R' if (cip == oechem.OECIPAtomStereo_R) else 'S'
                if cip != chirality:
                    # Note, I haven't seen this happen yet, but it shouldn't be a problem since there
                    # is only 2 directions for handedness and we're only running this for chiral atoms
                    raise Exception('Programming error: OpenEye atom stereochemistry assumptions failed.')

        # Set bond stereochemistry
        for bond, oebond in zip(self.atoms, oemol_bonds):
            if not bond.stereochemistry:
                continue

            oeatom1, oeatom2 = oemol_atoms[bond.atom1_index], oemol_atoms[bond.atom2_index]
            # Get two neighboring atoms
            oeatom1_neighbor = [n for n in oeatom1.GetAtoms()][0]
            oeatom2_neighbor = [n for n in oeatom2.GetAtoms()][0]
            # Define stereochemistry



            oestereo = oechem.OEBondStereo_Cis
            oebond.SetStereo([oeatom1, oeatom2], oechem.OEBondStereo_CisTrans, oestereo)

            # determine if stereochemistry is needed
            tag = rdb.GetStereo()
            if tag == Chem.BondStereo.STEREOCIS or tag == Chem.BondStereo.STEREOZ:
                stereo_atoms = rdb.GetStereoAtoms()
                stereo_bonds.append((newbond, stereo_atoms[0], stereo_atoms[1], oechem.OEBondStereo_Cis))

                bond2 = rdmol.GetBondBetweenAtoms(stereo_atoms[0], a1)
                bond4 = rdmol.GetBondBetweenAtoms(stereo_atoms[1], a2)
                print(tag, bond2.GetBondDir(), bond4.GetBondDir())
            if tag == Chem.BondStereo.STEREOTRANS or tag == Chem.BondStereo.STEREOE:
                stereo_atoms = rdb.GetStereoAtoms()
                stereo_bonds.append((newbond, stereo_atoms[0], stereo_atoms[1], oechem.OEBondStereo_Trans))
                bond2 = rdmol.GetBondBetweenAtoms(stereo_atoms[0], a1)
                bond4 = rdmol.GetBondBetweenAtoms(stereo_atoms[1], a2)
                print(tag, bond2.GetBondDir(), bond4.GetBondDir())

        # Set stereochemistry using the reference atoms extracted above
        for oeb, idx1, idx2, oestereo in stereo_bonds:
            oeb.SetStereo([map_atoms[idx1], map_atoms[idx2]], oechem.OEBondStereo_CisTrans, oestereo)

        # If the rdmol has a conformer, add its coordinates to the oemol
        # Note, this currently only adds the first conformer, it will need to be adjusted if the
        # you wanted to convert multiple sets of coordinates
        if rdmol.GetConformers():
            print("found an rdmol conformer")
            conf = rdmol.GetConformer()
            for rd_idx, oeatom in map_atoms.items():
                coords = conf.GetAtomPosition(rd_idx)
                oemol.SetCoords(oeatom, oechem.OEFloatArray(coords))

        # If RDMol has a title save it
        if rdmol.HasProp("_Name"):
            oemol.SetTitle(rdmol.GetProp("_Name"))

        # Clean Up phase
        # The only feature of a molecule that wasn't perceived above seemed to be ring connectivity, better to run it
        # here then for someone to inquire about ring sizes and get 0 when it shouldn't be
        oechem.OEFindRingAtomsAndBonds(oemol)

        print('Final Molecule: ', oechem.OEMolToSmiles(oemol))
        return oemol


    # TODO: We have to distinguish between retrieving user-specified partial charges and providing a generalized semiempirical/pop analysis/BCC scheme according to the new SMIRNOFF spec
    def get_partial_charges(self, method='AM1-BCC'):
        """Retrieve partial atomic charges.

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

        .. warning :: This API experimental and subject to change.

        """
        pass

    # TODO: Rework this.
    def _assign_partial_charges_using_sqm(sdf_filename, output_filename, charge_model="bcc"):
        """
        Assign partial charges with AmberTools using antechamber/sqm.

        Notes
        -----
        Currently only sdf file supported as input and mol2 as output
        https://github.com/choderalab/openmoltools/blob/master/openmoltools/packmol.py

        .. warning :: This API experimental and subject to change.

        """
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            raise(IOError("Antechamber not found, cannot run charge_mol()"))

        with temporary_directory() as tmp_dir:
            # TODO: Remove sqm files after completion
            # output_filename = tempfile.mktemp(suffix = ".mol2", dir = tmp_dir)
            os.system("antechamber -i {} -fi sdf -o {} -fo mol2 -pf y -c {}".format(sdf_filename, output_filename, charge_model))
            os.system("rm sqm.*")

    # TODO: Rework this.
    def _assign_partial_charges_using_quacpac(oemol, charge_model="bcc"):
        """
        Assign partial charges with OpenEye quacpac.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol
            OpenEye molecule for which partial charges are to be computed
        charge_model : str, optional, default='bcc'
            The charge model to use. One of ['AM1', 'AM1BCC', 'AM1BCCELF10', 'Gasteiger', 'MMFF94']

        Returns
        -------
        charges : numpy.array of shape (natoms) of type float
            The partial charges

        .. warning :: This API experimental and subject to change.

        """
        raise NotImplementedError()

    # TODO: Ditch this in favor of just-in-time computation of partial charges by get_partial_charges()
    @property
    def has_partial_charges(self):
        """Return True if any atom has nonzero charges; False otherwise.

        .. warning :: This API experimental and subject to change.

        """
        if (self.charges is None) or np.all(self.charges == 0.0):
            return False
        else:
            return True

    def get_fractional_bond_orders(self, method='Wiberg', toolkit=None, **kwargs):
        """Get fractional bond orders.

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

        .. warning :: This API experimental and subject to change.

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
