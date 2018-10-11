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
from collections import OrderedDict

from simtk import unit
from simtk.openmm.app import element

import openforcefield
from openforcefield.utils.toolkits import OPENEYE_AVAILABLE, RDKIT_AVAILABLE, AMBERTOOLS_AVAILABLE, GLOBAL_TOOLKIT_REGISTRY, SUPPORTED_FILE_FORMATS
#from openforcefield.utils.toolkits import requires_rdkit, requires_openeye
from openforcefield.typing.chemistry import ChemicalEnvironment, SMIRKSParsingError

# TODO: Do we need these?
from openforcefield.utils.toolkits import DEFAULT_AROMATICITY_MODEL, ALLOWED_AROMATICITY_MODELS
from openforcefield.utils.toolkits import DEFAULT_FRACTIONAL_BONDORDER_MODEL, ALLOWED_FRACTIONAL_BONDORDER_MODELS
from openforcefield.utils.toolkits import DEFAULT_CHARGE_MODEL, ALLOWED_CHARGE_MODELS

from openforcefield.utils.serialization import Serializable

from openforcefield.utils.toolkits import ToolkitRegistry, ToolkitWrapper, RDKitToolkitWrapper, OpenEyeToolkitWrapper

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
    ## From Jeff: This version of the API will push forward with immutable atoms and molecules
    def __init__(self, atomic_number, formal_charge, is_aromatic, stereochemistry=None,
                 name=None, molecule=None):
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
        self._molecule = molecule
        ## From Jeff: I'm going to assume that this is implicit in the parent Molecule's ordering of atoms
        #self._molecule_atom_index = molecule_atom_index
        self._bonds = list()


    def add_bond(self, bond):
        """Adds a bond that this atom is involved in
        .. todo :: Is this how we want to keep records?

        Parameters
        ----------
        bond: an openforcefield.topology.molecule.Bond
            A bond involving this atom

        """
        #TODO: I forgot if python will hold a pointer to the molecule, or a copy here
        self._bonds.append(bond)
        

    def to_dict(self):
        """Return a dict representation of the atom."""
        # TODO
        atom_dict = OrderedDict()
        atom_dict['atomic_number'] = self._atomic_number
        atom_dict['formal_charge'] = self._formal_charge
        atom_dict['is_aromatic'] = self._is_aromatic
        atom_dict['stereochemistry'] = self._stereochemistry
        # TODO: Should we let atoms have names?
        atom_dict['name'] = self._name
        # TODO: Should this be implicit in the atom ordering when saved?
        #atom_dict['molecule_atom_index'] = self._molecule_atom_index
        return atom_dict


            
    @classmethod
    def from_dict(cls, atom_dict):
        """Create an Atom from a dict representation."""
        ## TODO: classmethod or static method? Classmethod is needed for Bond, so it have
        ## its _molecule set and then look up the Atom on each side of it by ID
        return cls.__init__(*atom_dict)


    @property
    def formal_charge(self):
        """
        The atom's formal charge
        """
        return self._formal_charge

    @property
    def is_aromatic(self):
        """
        The atom's is_aromatic flag 
        """
        return self._is_aromatic

    @property
    def stereochemistry(self):
        """
        The atom's stereochemistry (if defined, otherwise None) 
        """
        return self._stereochemistry

    
    
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
        TODO (from jeff): Are there atoms that have different chemical properties based on their isotopes?
        
        """
        return self.element.mass

    # TODO: How are we keeping track of bonds, angles, etc?

    @property
    def bonds(self):
        """
        The list of ``Bond`` objects this atom is involved in.

        """
        for bond in self._bonds:
            yield bond

    @property
    def bonded_to(self):
        """
        The list of ``Atom`` objects this atom is involved in

        """
        for bond in self._bonds:
            for atom in bond.atoms:
                if not(atom == self):
                    # TODO: This seems dangerous. Ask John for a better way
                    yield atom

    @property
    def molecule(self):
        """
        The ``Molecule`` this atom is part of.

        .. todo::

           * Should we have a single unique ``Molecule`` for each molecule type in the system,
           or if we have multiple copies of the same molecule, should we have multiple ``Molecule``s?
        """
        return self._molecule


    
    @molecule.setter
    def molecule(self, molecule):
        """
        Set the atom's molecule pointer. Note that this will only work if the atom currently 
        doesn't have a molecule
        """
        #TODO: I forgot if python will hold a pointer to the molecule, or some sort of bulky copy here
        #TODO: Should _molecule be a property?
        #TODO: Add informative exception here
        assert self._molecule == None
        self._molecule = molecule

    
    @property
    def molecule_atom_index(self):
        """
        The index of this Atom within the the list of atoms in ``Molecules``.
        Note that this can be different from ``particle_index``.

        """
        if self._molecule is None:
            raise ValueError('This Atom does not belong to a Molecule object')
        return self._molecule.atoms.index(atom)

    ## From Jeff: Not sure if we actually need this
    @property
    def topology_atom_index(self):
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
    #def __init__(self, atom1, atom2, bondtype, fractional_bondorder=None):
    def __init__(self, atom1, atom2, bond_order, is_aromatic, stereochemistry=None):
        """
        Create a new chemical bond.
        """
        assert type(atom1) == Atom
        assert type(atom2) == Atom
        ## From Jeff: For now, I'm assuming each atom can only belong to one molecule
        # TODO: The molecule equality test just tests for identical SMILESes right now,
        # need to implement test for object identity
        assert atom1.molecule == atom2.molecule
        assert isinstance(atom1.molecule, FrozenMolecule)
        #if molecule == None:
        self._molecule = atom1.molecule
        ## From Jeff: Don't do this -- always perceive the molecule from the Atoms
        #elif type(molecule) == openforcefield.topology.molecule.Molecule:
        #    self._molecule = molecule
        #else:
        #    raise Exception("Bond created between atoms without a molecule specified, and no molecule provided")
            ## TODO: It might be interesting to think about some sort of architecture where
        ## creating a bond between two Molecules could turn them into one molecule...
        
        self._atom1 = atom1
        self._atom2 = atom2
        
        atom1.add_bond(self)
        atom2.add_bond(self)
        # TODO: Check bondtype and fractional_bondorder are valid?
        #self._type = bondtype
        #self._fractional_bondorder = fractional_bondorder
        self._bond_order = bond_order
        self._is_aromatic = is_aromatic
        self._stereochemistry = stereochemistry
        #self._molecule = molecule


        
    def to_dict(self):
        """Return a dict representation of the bond."""
        bond_dict = OrderedDict()
        bond_dict['atom1'] = self._molecule.atoms.index(self._atom1)
        bond_dict['atom2'] = self._molecule.atoms.index(self._atom2)
        bond_dict['bond_order'] = self._bond_order
        bond_dict['is_aromatic'] = self._is_aromatic
        bond_dict['stereochemistry'] = self._stereochemistry
        return bond_dict

    @classmethod
    def from_dict(cls, molecule, d):
        """Create a Bond from a dict representation."""
        # TODO
        d['molecule'] = molecule
        d['atom1'] = molecule.atoms[d['atom1']]
        d['atom2'] = molecule.atoms[d['atom2']]
        return cls(*d)        
        #self.__init__(atom1, atom2, bondorder, 
        #self._molecule = molecule
        #atom1_index = d['atom1']
        #atom2_index = d['atom2']
        #atom1 = self._molecule.atoms(atom1_index)
        #atom2 = self._molecule.atoms(atom2_index)
        #bond_order = d['bond_order']
        #is_aromatic = d['is_aromatic']
        #stereochemistry = d['stereochemistry']

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    @property
    def atoms(self):
        return (self._atom1, self._atom2)

    #def type(self):
    #    return self._type

    @property
    def bond_order(self):
        return self._bond_order

    @bond_order.setter
    def bond_order(self, value):
        self._bond_order = value

    @property
    def stereochemistry(self):
        return self._stereochemistry
    @property
    def is_aromatic(self):
        return self._is_aromatic

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, value):
        """
        Sets the Bond's parent molecule. Can not be changed after assignment
        """
        assert self._molecule == None
        self._molecule = value

    @property
    def molecule_bond_index(self):
        """
        The index of this Bond within the the list of bonds in ``Molecules``.

        """
        if self._molecule is None:
            raise ValueError('This Atom does not belong to a Molecule object')
        return self._molecule.bonds.index(bond)
    
    #@property
    #def fractional_bondorder(self):
    #    return self._fractional_bondorder

    #@fractional_bondorder.setter
    #def fractional_bondorder(self, value):
    #    self._fractional_bondorder = value

#=============================================================================================
# Molecule
#=============================================================================================

# TODO: How do we automatically trigger invalidation of cached properties if an ``Atom``, ``Bond``, or ``VirtualSite`` is modified,
#       rather than added/deleted via the API? The simplest resolution is simply to make them immutable.

class FrozenMolecule(Serializable):
    """
    Immutable chemical representation of a molecule, such as a small molecule or biopolymer.

    .. todo :: What other API calls would be useful for supporting biopolymers as small molecules? Perhaps iterating over chains and residues?

    Examples
    --------

    Create a molecule from a mol2 file

    >>> molecule = FrozenMolecule.from_file('molecule.mol2')

    Create a molecule from an OpenEye molecule

    >>> molecule = FrozenMolecule.from_openeye(oemol)

    Create a molecule from an RDKit molecule

    >>> molecule = FrozenMolecule.from_rdkit(rdmol)

    Create a molecule from IUPAC name (requires the OpenEye toolkit)

    >>> molecule = FrozenMolecule.from_iupac('imatinib')

    Create a molecule from SMILES

    >>> molecule = FrozenMolecule.from_smiles('Cc1ccccc1')

    """
    def __init__(self, other=None):
        """
        Create a new FrozenMolecule object

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

        >>> empty_molecule = FrozenMolecule()

        Create a molecule from another molecule:

        >>> molecule_copy = FrozenMolecule(molecule_copy)

        Create a molecule from a file that can be used to construct a molecule,
        using either a filename or file-like object:

        >>> molecule = FrozenMolecule('molecule.mol2')
        >>> molecule = FrozenMolecule(open('molecule.mol2', 'r'))
        >>> molecule = FrozenMolecule(gzip.GzipFile('molecule.mol2.gz', 'r'))

        Create a molecule from an OpenEye molecule:

        >>> molecule = FrozenMolecule(oemol)

        Create a molecule from an RDKit molecule:

        >>> molecule = FrozenMolecule(rdmol)

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
            elif OPENEYE_AVAILABLE and issubclass(other, openeye.oechem.OEMolBase):
                mol = Molecule.from_openeye(other)
                self._copy_initializer(mol)
            elif RDKIT_AVAILABLE and isinstance(other, rdkit.Chem.rdchem.Mol):
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
        molecule_dict['name'] = self._name
        ## From Jeff: If we go the properties-as-dict route, then _properties should, at
        ## the top level, be a dict. Should we go through recursively and ensure all values are dicts too?
        molecule_dict['atoms'] = [ atom.to_dict() for atom in self._atoms ]
        molecule_dict['virtual_sites'] = [ vsite.to_dict() for vsite in self._virtual_sites ]
        molecule_dict['bonds'] = [ bond.to_dict() for bond in self._bonds ]
        # TODO: How do we make sure bonds are serializable?
        ## From Jeff: Maybe each molecule could have its own atomic ID scheme?
        # TODO: Charges
        # TODO: Properties
        ## From Jeff: We could have the onerous requirement that all "properties" have to_dict() functions.
        ## Or we could restrict properties to simple stuff (ints, strings, floats, and the like)
        ## Or pickle anything unusual
        ## Or not allow user-defined properties at all (just use our internal _cached_properties)
        #molecule_dict['properties'] = dict([(key, value._to_dict()) for key.value in self._properties])
        # TODO: Assuming "simple stuff" properties right now, figure out a better standard
        molecule_dict['properties'] = self._properties
        if hasattr(self, '_cached_properties'):
            molecule_dict['cached_properties'] = self._cached_properties
        # TODO: Conformers
        return molecule_dict

    @classmethod
    def from_dict(cls, molecule_dict):
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
        molecule = cls(name=molecule_dict['name'])
        for atom_dict in molecule_dict['atoms']:
            #atom = Atom.from_dict(atom_dict)
            ## From Jeff: This is really ugly
            #atom.set_molecule(self)
            #self._atoms.append(atom)
            molecule.add_atom(*atom_dict)
        #for vsite in molecule_dict['virtual_sites']:
        #    molecule._add_virtual_site(*vsite)
        for bond_dict in molecule_dict['bonds']:
            ## TODO: There's probably a more graceful way to do this than overwriting a bond_dict entry
            bond_dict['atom1'] = molecule.atoms[bond_dict['atom1']]
            bond_dict['atom2'] = molecule.atoms[bond_dict['atom2']]
            molecule.add_bond(*bond_dict)
            #bond = Bond.from_dict(bond_dict)
            #bond.set_molecule(
            #self._bonds.append(bond)
            
            #molecule._add_bond(self, bond_dict = bond) # TODO: How is this correctly handled? using indices?
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
        ## From Jeff:
        ## Pro: More informative error messages for users when there are parameterization problems
        ## Pro: Easy way to find molecules of interest in a Topology or other collection of Molecules
        ## Pro: More informative when manually inspecting a serialized file
        ## Con: We'd need to make rules about duplicate names in a Topology
        ## Con: If we have a Topology.select_atoms() function, people would expect names to be searchable and we'd have to implement that
        ## Opinion: We should store the name in _properties
        self._atoms = list()
        self._virtual_sites = list()
        #self._particles = list() # List of particles (atoms or virtual sites) # TODO: Should this be a dict?
        self._bonds = list() # List of bonds between Atom objects # TODO: Should this be a dict?
        #self._charges = None # TODO: Storage charges
        self._properties = None # Attached properties to be preserved
        self._cached_properties = None # Cached properties (such as partial charges) can be recomputed as needed
        #self._conformers = None # Optional conformers
        ## From Jeff: Why is it necessary to know mutliple conformers for the same molecule?
        ## I get the sense that Molecule objects should just deal with connectivity, not geometry.
        ## Having a _conformers field would make sense if we're expecting cases where we assign
        ## different parameters to molecules based on their conformation, but then we might need
        ## a whole geometry implementation in here.

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

    def __eq__(self, other):
        """Test two molecules for equality to see if they are the chemical species, but do not check other annotated properties.

        .. note ::

           Note that this method simply tests whether two molecules are identical chemical species using equivalence of their canonical isomeric SMILES.
           No effort is made to ensure that the atoms are in the same order or that any annotated properties are preserved.

        """
        return self.to_smiles() == other.to_smiles()

    def to_smiles(self, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Return a canonical isomeric SMILES representation of the current molecule

        .. note :: RDKit and OpenEye versions will not necessarily return the same representation.

        .. todo :: Can we ensure RDKit and OpenEye versions return the same representation?

        Parameters
        ----------
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES conversion

        Returns
        -------
        smiles : str
            Canonical isomeric explicit-hydrogen SMILES

        Examples
        --------

        >>> smiles = molecule.to_smiles()

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call('to_smiles', self)
        elif isinstance(toolkit_registry, ToolkitWrapper):
            return toolkit.to_smiles(self)
        #return toolkit_registry.call('to_smiles', self)

    @staticmethod
    def from_smiles(smiles, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Construct a Molecule from a SMILES representation

        Parameters
        ----------
        smiles : str
            The SMILES representation of the molecule.

        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion

        Returns
        -------
        molecule : Molecule
            The molecule

        Examples
        --------

        >>> molecule = Molecule.from_smiles('Cc1ccccc1')

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call('from_smiles', smiles)
        elif isinstance(toolkit_registry, ToolkitWrapper):
            return toolkit_registry.from_smiles(self, smiles)
        #return tookit_registry.call('from_smiles', smiles)

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

    def _add_atom(self, atomic_number, formal_charge, is_aromatic, stereochemistry=None,
                  name=None):
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
        atom = Atom(atomic_number, formal_charge, is_aromatic, stereochemistry=stereochemistry, name=name, molecule=self)
        self._atoms.append(atom)
        #self._particles.append(atom)
        self._invalidate_cached_properties()
        return self._atoms.index(atom)

    def _add_bond(self, atom1, atom2, bond_order, is_aromatic, stereochemistry=None):
        """
        Add a bond between two specified atom indices

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        atom1 : int or openforcefield.topology.molecule.Atom
            Index of first atom or first atom
        atom2_index : int or openforcefield.topology.molecule.Atom
            Index of second atom or second atom
        bond_order : int
            Integral bond order of Kekulized form
        is_aromatic : bool
            True if this bond is aromatic, False otherwise
        stereochemistry : str, optional, default=None
            Either 'E' or 'Z' for specified stereochemistry, or None if stereochemistry is irrelevant

        Returns
        -------
        index : int
            The index of the bond in the molecule

        """
        if isinstance(atom1, int) and isinstance(atom2, int):
            atom1_atom = self.atoms[atom1]
            atom2_atom = self.atoms[atom2]
        elif isinstance(atom1, Atom) and isinstance(atom2, Atom):
            atom1_atom = atom1
            atom2_atom = atom2
        else:
            raise Exception('Invalid inputs to molecule._add_bond. Expected ints or Atoms. Received {} and {} '.format(type(atom1), type(atom2), atom1))
        # TODO: Check to make sure bond does not already exist
        bond = Bond(atom1_atom, atom2_atom, bond_order, is_aromatic, stereochemistry=stereochemistry)
        # TODO: This is a bad way to get bond index
        #bond_index = len(self._bonds)
        #bond.set_molecule(self, bond_index)
        #atom1.add_bond(bond)
        #atom2.add_bond(bond)
        self._bonds.append(bond)
        self._invalidate_cached_properties()
        return self._bonds.index(bond)

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
        return self._atoms

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
        return self._bonds

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

    def chemical_environment_matches(self, query, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """Retrieve all matches for a given chemical environment query.

        .. todo ::

           * Do we want to generalize ``query`` to allow other kinds of queries, such as mdtraj DSL, pymol selections, atom index slices, etc?
             We could call it ``topology.matches(query)`` instead of ``chemical_environment_matches``

        Parameters
        ----------
        query : str or ChemicalEnvironment
            SMARTS string (with one or more tagged atoms) or ``ChemicalEnvironment`` query
            Query will internally be resolved to SMIRKS using ``query.asSMIRKS()`` if it has an ``.asSMIRKS`` method.
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for chemical environment matches
            

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
        # TODO: Simplify this by requiring a toolkit registry for the molecule?
        # TODO: Do we have to pass along an aromaticity model?
        if isinstance(toolkit, ToolkitRegistry):
            matches = toolkit.call('find_smarts_matches', self, smirks)
        elif isinstance(toolkit, ToolkitWrapper):
            matches = toolkit.find_smarts_matches(self, smirks)
        else:
            raise ValueError("'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper")
            #else:
        #    matches = GLOBAL_TOOLKIT_REGISTRY.find_smarts_matches(self, smirks)

        return matches

    @property
    def name(self):
        """The name (or title) of the molecule
        """
        return self._name

    @staticmethod    
    #@requires_openeye('oechem', 'oeiupac')
    @OpenEyeToolkitWrapper.requires_toolkit()
    #@OpenEyeToolkitWrapper.requires_toolkit(OpenEyeToolkitWrapper)
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

    #@requires_openeye('oechem', 'oeiupac')
    @OpenEyeToolkitWrapper.requires_toolkit()
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
    def from_file(filename, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Create one or more molecules from a file

        .. todo::

           * Extend this to also include some form of .offmol Open Force Field Molecule format?
           * Generalize this to also include file-like objects?

        Parameters
        ----------
        filename : str
            The name of the file to stream one or more molecules from.

        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for file loading. If a Toolkit is passed, only the highest-precedence toolkit is used
            
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
        # TODO: This needs to be cleaned up to use the new ToolkitRegistry and ToolkitWrappers
        # TODO: Check file extensions for compatible toolkit wrappers

        # Use highest-precendence toolkit
        #toolkit = TOOLKIT_PRECEDENCE[0]
        if isinstance(toolkit, ToolkitRegistry):
            toolkit = toolkit_registry.registered_toolkits[0]
        elif isinstance(toolkit, ToolkitWrapper):
            toolkit = toolkit_registry
        else:
            raise ValueError("'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper")
        #toolkit = 
        mols = list()
        
        if type(toolkit) is OpenEyeToolkitWrapper:
            # Read molecules from an OpenEye-supported file, converting them one by one
            from openeye import oechem
            oemol = oechem.OEGraphMol()
            ifs = oechem.oemolistream(filename)
            while oechem.OEReadMolecule(ifs, oemol):
                mol = Molecule.from_openeye(oemol)
                mols.append(mol)
        elif type(toolkit) is RDKitToolkitWrapper:
            from rdkit import Chem
            for rdmol in Chem.SupplierFromFilename(filename):
                mol = Molecule.from_rdkit(rdmol)
                mols.append(mol)
        else:
            raise Exception('Toolkit {} unsupported.'.format(toolkit))
        return mols

    def to_file(self, outfile, format, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """Write the current molecule to a file or file-like object

        Parameters
        ----------
        outfile : str or file-like object
            A file-like object or the filename of the file to be written to
        format : str
            Format specifier, one of ['MOL2', 'MOL2H', 'SDF', 'PDB', 'SMI', 'CAN', 'TDT']
            Note that not all toolkits support all formats
        toolkit_registry : openforcefield.utils.toolkits.ToolRegistry or openforcefield.utils.toolkits.ToolkitWrapper, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for file writing. If a Toolkit is passed, only the highest-precedence toolkit is used

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
        # TODO: This needs to be cleaned up to use the new ToolkitRegistry and ToolkitWrappers

        # Determine which formats are supported
        toolkit = None
        for query_toolkit in toolkit_registry:
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
    @RDKitToolkitWrapper.requires_toolkit()
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
        toolkit = RDKitToolkitWrapper()
        return toolkit.from_rdkit(rdmol)

    @RDKitToolkitWrapper.requires_toolkit()
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
        toolkit = RDKitToolkitWrapper()
        return toolkit.to_rdkit(self, aromaticity_model=aromaticity_model)

    @staticmethod
    @OpenEyeToolkitWrapper.requires_toolkit()
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
        toolkit = OpenEyeToolkitWrapper()
        return toolkit.from_openeye(oemol)

    #@OpenEyeToolkitWrapper.requires_toolkit(OpenEyeToolkitWrapper)
    @OpenEyeToolkitWrapper.requires_toolkit()
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
        toolkit = OpenEyeToolkitWrapper()
        return toolkit.to_openeye(self, aromaticity_model=aromaticity_model)

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


class Molecule(FrozenMolecule):
    """
    Mutable chemical representation of a molecule, such as a small molecule or biopolymer.

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
    def __init__(self, *args, **kwargs):
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
        #super(self, Molecule).__init__(*args, **kwargs)
        super(Molecule, self).__init__(*args, **kwargs)

    def add_atom(self, atomic_number, formal_charge, is_aromatic, stereochemistry=None,
                 name=None):
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
        atom_index = self._add_atom(atomic_number, formal_charge, is_aromatic, stereochemistry=stereochemistry, name=name)
        return atom_index
        
    def add_bond(self, atom1_index, atom2_index, bond_order, is_aromatic,
                 stereochemistry=None):
        """
        Add a bond between two specified atom indices

        .. warning :: This API experimental and subject to change.

        Parameters
        ----------
        atom1_index : int
            Index of first atom
        atom2_index : int
            Index of second atom
        bond_order : int
            Integral bond order of Kekulized form
        is_aromatic : bool
            True if this bond is aromatic, False otherwise
        stereochemistry : str, optional, default=None
            Either 'E' or 'Z' for specified stereochemistry, or None if stereochemistry is irrelevant

        """
        bond_index = self._add_bond(atom1_index, atom2_index, bond_order, is_aromatic, stereochemistry=stereochemistry)
