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
   * Refactor toolkit encapsulation to generalize and provide only a few major toolkit methods and toolkit objects
        that can be queried for features
   * Speed up overall import time by putting non-global imports only where they are needed

"""
import json
import operator
import warnings
from collections import OrderedDict, UserDict
from copy import deepcopy
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    Generator,
    List,
    Optional,
    Set,
    TextIO,
    Tuple,
    Union,
)

import networkx as nx
import numpy as np
from openff.units import unit
from openff.units.elements import MASSES, SYMBOLS
from openff.utilities.exceptions import MissingOptionalDependencyError
from packaging import version

import openff.toolkit
from openff.toolkit.utils.exceptions import (
    HierarchySchemeNotFoundException,
    HierarchySchemeWithIteratorNameAlreadyRegisteredException,
    IncompatibleUnitError,
    InvalidAtomMetadataError,
    InvalidConformerError,
    SmilesParsingError,
    UnsupportedFileTypeError,
)
from openff.toolkit.utils.serialization import Serializable
from openff.toolkit.utils.toolkits import (
    DEFAULT_AROMATICITY_MODEL,
    GLOBAL_TOOLKIT_REGISTRY,
    InvalidToolkitRegistryError,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
    ToolkitRegistry,
    ToolkitWrapper,
    UndefinedStereochemistryError,
)
from openff.toolkit.utils.utils import get_data_file_path, requires_package

if TYPE_CHECKING:
    from openff.units.unit import Quantity

    from openff.toolkit.topology._mm_molecule import _SimpleAtom, _SimpleMolecule

# TODO: Can we have the `ALLOWED_*_MODELS` list automatically appear in the docstrings below?
# TODO: Should `ALLOWED_*_MODELS` be objects instead of strings?
# TODO: Should these be imported from `openff.toolkit.cheminformatics.aromaticity_models` and `.bondorder_models`?

# TODO: Allow all OpenEye aromaticity models to be used with OpenEye names?
#       Only support OEAroModel_MDL in RDKit version?


def _molecule_deprecation(old_method, new_method):
    warnings.warn(
        f"Molecule.{old_method} is deprecated. Use Molecule.{new_method} instead.",
        MoleculeDeprecationWarning,
    )


class MoleculeDeprecationWarning(UserWarning):
    """Warning for deprecated portions of the Molecule API."""


class Particle(Serializable):
    """
    Base class for all particles in a molecule.

    A particle object could be an ``Atom`` or similar.

    .. warning :: This API is experimental and subject to change.
    """

    @property
    def molecule(self):
        r"""
        The ``Molecule`` this particle is part of.

        .. todo::

            * Should we have a single unique ``Molecule`` for each molecule
              type in the system, or if we have multiple copies of the same
              molecule, should we have multiple ``Molecule``\ s?

        """
        return self._molecule

    @molecule.setter
    def molecule(self, molecule):
        """
        Set the particle's molecule pointer. Note that this will only work if the particle currently
        doesn't have a molecule
        """
        err = f"{type(self).__name__} already has an associated molecule"
        assert self._molecule is None, err
        self._molecule = molecule

    @property
    def molecule_particle_index(self):
        """
        Returns the index of this particle in its molecule
        """
        return self._molecule.particles.index(self)

    @property
    def name(self):
        """
        The name of the particle
        """
        return self._name

    def to_dict(self):
        """Convert to dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO

    @classmethod
    def from_dict(cls, d):
        """Static constructor from dictionary representation."""
        # Implement abstract method Serializable.to_dict()
        raise NotImplementedError()  # TODO


class AtomMetadataDict(UserDict):
    def __init__(self, *args, **kwargs):
        self.data = {}
        self.update(dict(*args, **kwargs))

    def __setitem__(self, key, value):
        if not isinstance(key, str):
            raise InvalidAtomMetadataError(
                f"Attempted to set atom metadata with a non-string key. (key: {key}"
            )
        if not isinstance(value, (str, int)):
            raise InvalidAtomMetadataError(
                f"Attempted to set atom metadata with a non-string or integer "
                f"value. (value: {value})"
            )
        super().__setitem__(key, value)


class Atom(Particle):
    """
    A chemical atom.

    .. todo::

       * Should ``Atom`` objects be immutable or mutable?
       * Do we want to support the addition of arbitrary additional properties,
         such as floating point quantities (e.g. ``charge``), integral
         quantities (such as ``id`` or ``serial`` index in a PDB file),
         or string labels (such as Lennard-Jones types)?

    .. todo :: Allow atoms to have associated properties.

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(
        self,
        atomic_number,
        formal_charge,
        is_aromatic,
        name=None,
        molecule=None,
        stereochemistry=None,
        metadata=None,
    ):
        """
        Create an immutable Atom object.

        Object is serializable and immutable.

        .. todo :: Use attrs to validate?

        .. todo :: We can add setters if we need to.

        Parameters
        ----------
        atomic_number : int
            Atomic number of the atom
        formal_charge : int or openff.units.unit.Quantity-wrapped int with dimension "charge"
            Formal charge of the atom
        is_aromatic : bool
            If True, atom is aromatic; if False, not aromatic
        stereochemistry : str, optional, default=None
            Either 'R' or 'S' for specified stereochemistry, or None for ambiguous stereochemistry
        name : str, optional, default=None
            An optional name to be associated with the atom
        metadata : dict[str: (int, str)], default=None
            An optional dictionary where keys are strings and values are strings or ints. This is intended
            to record atom-level information used to inform hierarchy definition and iteration, such as
            grouping atom by residue and chain.

        Examples
        --------

        Create a non-aromatic carbon atom

        >>> atom = Atom(6, 0, False)

        Create a chiral carbon atom

        >>> atom = Atom(6, 0, False, stereochemistry='R', name='CT')

        """
        self._atomic_number = atomic_number
        # Use the setter here, since it will handle either ints or Quantities
        if hasattr(formal_charge, "units"):
            # Faster check than ` == unit.dimensionless`
            if str(formal_charge.units) == "":
                raise Exception
        self.formal_charge = formal_charge
        self._is_aromatic = is_aromatic
        self._stereochemistry = stereochemistry
        if name is None:
            name = ""
        self._name = name
        self._molecule = molecule
        # From Jeff: I'm going to assume that this is implicit in the parent Molecule's ordering of atoms
        # self._molecule_atom_index = molecule_atom_index
        self._bonds = list()

        if metadata is None:
            self._metadata = AtomMetadataDict()
        else:
            self._metadata = AtomMetadataDict(metadata)

    # TODO: We can probably avoid an explicit call and determine this dynamically
    #   from self._molecule (maybe caching the result) to get rid of some bookkeeping.
    # TODO: Should stereochemistry be reset/cleared/recomputed upon addition of a bond?
    def add_bond(self, bond):
        """Adds a bond that this atom is involved in
        .. todo :: Is this how we want to keep records?

        Parameters
        ----------
        bond: an openff.toolkit.topology.molecule.Bond
            A bond involving this atom
        """

        self._bonds.append(bond)

    def to_dict(self):
        """Return a dict representation of the atom."""
        # TODO
        atom_dict = OrderedDict()
        atom_dict["atomic_number"] = self._atomic_number
        atom_dict["formal_charge"] = self._formal_charge.m_as(unit.elementary_charge)
        atom_dict["is_aromatic"] = self._is_aromatic
        atom_dict["stereochemistry"] = self._stereochemistry
        # TODO: Should we let atoms have names?
        atom_dict["name"] = self._name
        atom_dict["metadata"] = dict(self._metadata)
        # TODO: Should this be implicit in the atom ordering when saved?
        # atom_dict['molecule_atom_index'] = self._molecule_atom_index
        return atom_dict

    @classmethod
    def from_dict(cls, atom_dict):
        """Create an Atom from a dict representation."""
        return cls(**atom_dict)

    @property
    def metadata(self):
        """
        The atom's metadata dictionary
        """
        return self._metadata

    @property
    def formal_charge(self):
        """
        The atom's formal charge
        """
        return self._formal_charge

    @formal_charge.setter
    def formal_charge(self, other):
        """
        Set the atom's formal charge. Accepts either ints or unit-wrapped ints with units of charge.
        """
        if isinstance(other, int):
            self._formal_charge = unit.Quantity(other, unit.elementary_charge)
        elif isinstance(other, unit.Quantity):
            # Faster to check equality than convert, so short-circuit
            if other.units is unit.elementary_charge:
                self.formal_charge = other
            elif other.units in unit.elementary_charge.compatible_units():
                self._formal_charge = other
            else:
                raise IncompatibleUnitError(
                    f"Cannot set formal charge with a quantity with units {other.units}"
                )
        elif hasattr(other, "unit"):
            from openmm import unit as openmm_unit

            if not isinstance(other, openmm_unit.Quantity):
                raise IncompatibleUnitError(
                    "Unsupported type passed to formal_charge setter. "
                    "Found object of type {type(other)}."
                )

            from openff.units.openmm import from_openmm

            converted = from_openmm(other)
            if converted.units in unit.elementary_charge.compatible_units():
                self._formal_charge = converted
            else:
                raise IncompatibleUnitError(
                    f"Cannot set formal charge with a quantity with units {converted.units}"
                )
        else:
            raise ValueError

    @property
    def partial_charge(self):
        """
        The partial charge of the atom, if any.

        Returns
        -------
        unit-wrapped float with dimension of atomic charge, or None if no charge has been specified
        """
        if self._molecule._partial_charges is None:
            return None
        else:
            index = self.molecule_atom_index
            return self._molecule._partial_charges[index]

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

    @stereochemistry.setter
    def stereochemistry(self, value):
        """Set the atoms stereochemistry
        Parameters
        ----------
        value : str
            The stereochemistry around this atom, allowed values are "CW", "CCW", or None,
        """

        # if (value != 'CW') and (value != 'CCW') and not(value is None):
        #    raise Exception(
        #       "Atom stereochemistry setter expected 'CW', 'CCW', or None. ""
        #       "Received {} (type {})".format(value, type(value))"
        # )
        self._stereochemistry = value

    @property
    def atomic_number(self) -> int:
        """
        The integer atomic number of the atom.

        """
        return self._atomic_number

    @property
    def symbol(self) -> str:
        """
        Return the symbol implied by the atomic number of this atom

        """
        return SYMBOLS[self.atomic_number]

    @property
    def mass(self) -> "Quantity":
        """
        The standard atomic weight (abundance-weighted isotopic mass) of the atomic site.

        The mass is reported in units of Dalton.
        """
        # This is assumed elsewhere in the codebase to be in units of Dalton, which is what is
        # reported by MASSES as of openff-units v0.1.5. There may be performance implications if
        # other functions need to verify or convert units.
        # https://github.com/openforcefield/openff-toolkit/pull/1182#discussion_r802078273
        return MASSES[self.atomic_number]

    @property
    def name(self):
        """
        The name of this atom, if any
        """
        return self._name

    @name.setter
    def name(self, other):
        """

        Parameters
        ----------
        other : string
            The new name for this atom
        """
        if type(other) != str:
            raise Exception(
                f"In setting atom name. Expected str, received {other} (type {type(other)})."
            )
        self._name = other

    # TODO: How are we keeping track of bonds, angles, etc?

    @property
    def bonds(self):
        """
        The list of ``Bond`` objects this atom is involved in.

        """
        return self._bonds
        # for bond in self._bonds:
        #    yield bond

    @property
    # def bonded_to(self):
    def bonded_atoms(self):
        """
        The list of ``Atom`` objects this atom is involved in bonds with

        """
        for bond in self._bonds:
            for atom in bond.atoms:
                if atom is not self:
                    # TODO: This seems dangerous. Ask John for a better way
                    yield atom

    def is_bonded_to(self, atom2):
        """
        Determine whether this atom is bound to another atom

        Parameters
        ----------
        atom2: openff.toolkit.topology.molecule.Atom
            a different atom in the same molecule

        Returns
        -------
        bool
            Whether this atom is bound to atom2
        """
        # TODO: Sanity check (check for same molecule?)
        assert self != atom2
        for bond in self._bonds:
            for bonded_atom in bond.atoms:
                if atom2 == bonded_atom:
                    return True
        return False

    def is_in_ring(self, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY) -> bool:
        """
        Return whether or not this atom is in a ring(s) (of any size)

        This Atom is expected to be attached to a molecule (`Atom.molecule`).

        Parameters
        ----------
        toolkit_registry: openff.toolkit.utils.toolkits.ToolkitRegistry, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` to use to enumerate the tautomers.

        """
        _is_in_ring = toolkit_registry.call("atom_is_in_ring", self)

        return _is_in_ring

    @property
    def molecule_atom_index(self):
        """
        The index of this Atom within the the list of atoms in the parent ``Molecule``.
        """
        if self._molecule is None:
            raise ValueError("This Atom does not belong to a Molecule object")
        if "_molecule_atom_index" in self.__dict__:
            return self._molecule_atom_index
        self._molecule_atom_index = self._molecule.atoms.index(self)
        return self._molecule_atom_index

    def __repr__(self):
        # TODO: Also include which molecule this atom belongs to?
        return f"Atom(name={self._name}, atomic number={self._atomic_number})"

    def __str__(self):
        # TODO: Also include which molecule this atom belongs to?
        return "<Atom name='{}' atomic number='{}'>".format(
            self._name, self._atomic_number
        )


# =============================================================================================
# Bond Stereochemistry
# =============================================================================================

# class BondStereochemistry(Serializable):
# """
# Bond stereochemistry representation
# """
# def __init__(self, stereo_type, neighbor1, neighbor2):
#    """
#
#    Parameters
#    ----------
#    stereo_type
#    neighbor1
#    neighbor2
#    """
#    assert isinstance(neighbor1, Atom)
#    assert isinstance(neighbor2, Atom)
#    # Use stereo_type @setter to check stereo type is a permitted value
#    self.stereo_type = stereo_type
#    self._neighbor1 = neighbor1
#    self._neighbor2 = neighbor2

# def to_dict(self):
#    bs_dict = OrderedDict()
#    bs_dict['stereo_type'] = self._stereo_type
#    bs_dict['neighbor1_index'] = self._neighbor1.molecule_atom_index
#    bs_dict['neighbor2_index'] = self._neighbor2.molecule_atom_index
#    return bs_dict

# classmethod
# def from_dict(cls, molecule, bs_dict):
#    neighbor1 = molecule.atoms[bs_dict['neighbor1_index']]
#    neighbor2 = molecule.atoms[bs_dict['neighbor2_index']]
#    return cls.__init__(bs_dict['stereo_type'], neighbor1, neighbor2)

# @property
# def stereo_type(self):
#    return self._stereo_type

# @stereo_type.setter
# def stereo_type(self, value):
#    assert (value == 'CIS') or (value == 'TRANS') or (value is None)
#    self._stereo_type = value

# @property
# def neighbor1(self):
#    return self._neighbor1

# @property
# def neighbor2(self):
#    return self._neighbor2

# @property
# def neighbors(self):
#    return (self._neighbor1, self._neighbor2)


class Bond(Serializable):
    """
    Chemical bond representation.

    .. warning :: This API is experimental and subject to change.

    .. todo :: Allow bonds to have associated properties.

    Attributes
    ----------
    atom1, atom2 : openff.toolkit.topology.Atom
        Atoms involved in the bond
    bond_order : int
        The (integer) bond order of this bond.
    is_aromatic : bool
        Whether or not this bond is aromatic.
    fractional_bond_order : float, optional
        The fractional bond order, or partial bond order of this bond.
    stereochemstry : str, optional, default=None
        A string representing this stereochemistry of this bond.

    .. warning :: This API is experimental and subject to change.
    """

    def __init__(
        self,
        atom1,
        atom2,
        bond_order,
        is_aromatic,
        fractional_bond_order=None,
        stereochemistry=None,
    ):
        """
        Create a new chemical bond.

        """
        assert type(atom1) == Atom
        assert type(atom2) == Atom
        assert atom1.molecule is atom2.molecule
        assert isinstance(atom1.molecule, FrozenMolecule)
        self._molecule = atom1.molecule

        self._atom1 = atom1
        self._atom2 = atom2

        atom1.add_bond(self)
        atom2.add_bond(self)
        # TODO: Check bondtype and fractional_bond_order are valid?
        # TODO: Dative bonds
        self._fractional_bond_order = fractional_bond_order
        self._bond_order = bond_order
        self._is_aromatic = is_aromatic
        self._stereochemistry = stereochemistry

    def to_dict(self):
        """
        Return a dict representation of the bond.

        """
        bond_dict = OrderedDict()
        bond_dict["atom1"] = self.atom1.molecule_atom_index
        bond_dict["atom2"] = self.atom2.molecule_atom_index
        bond_dict["bond_order"] = self._bond_order
        bond_dict["is_aromatic"] = self._is_aromatic
        bond_dict["stereochemistry"] = self._stereochemistry
        bond_dict["fractional_bond_order"] = self._fractional_bond_order
        return bond_dict

    @classmethod
    def from_dict(cls, molecule, d):
        """Create a Bond from a dict representation."""
        # TODO
        d["molecule"] = molecule
        d["atom1"] = molecule.atoms[d["atom1"]]
        d["atom2"] = molecule.atoms[d["atom2"]]
        return cls(*d)

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    @property
    def atom1_index(self):
        return self.molecule.atoms.index(self._atom1)

    @property
    def atom2_index(self):
        return self.molecule.atoms.index(self._atom2)

    @property
    def atoms(self):
        return (self._atom1, self._atom2)

    @property
    def bond_order(self):
        return self._bond_order

    @bond_order.setter
    def bond_order(self, value):
        self._bond_order = value

    @property
    def fractional_bond_order(self):
        return self._fractional_bond_order

    @fractional_bond_order.setter
    def fractional_bond_order(self, value):
        self._fractional_bond_order = value

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
        assert self._molecule is None
        self._molecule = value

    @property
    def molecule_bond_index(self):
        """
        The index of this Bond within the the list of bonds in ``Molecules``.

        """
        if self._molecule is None:
            raise ValueError("This Atom does not belong to a Molecule object")
        return self._molecule.bonds.index(self)

    def is_in_ring(self, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY) -> bool:
        """
        Return whether or not this bond is in a ring(s) (of any size)

        This Bond is expected to be attached to a molecule (`Bond.molecule`).

        Note: Bonds containing atoms that are only in separate rings, i.e. the central bond in a biphenyl,
            are not considered to be bonded by this criteria.

        Parameters
        ----------
        toolkit_registry: openff.toolkit.utils.toolkits.ToolkitRegistry, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` to use to enumerate the tautomers.

        Returns
        -------
        is_in_ring: bool
            Whether or not this bond is in a ring.

        """
        _is_in_ring = toolkit_registry.call("bond_is_in_ring", self)

        return _is_in_ring

    def __repr__(self):
        return f"Bond(atom1 index={self.atom1_index}, atom2 index={self.atom2_index})"

    def __str__(self):
        return (
            f"<Bond atom1 index='{self.atom1_index}', atom2 index='{self.atom2_index}'>"
        )


# TODO: How do we automatically trigger invalidation of cached properties if an ``Atom`` or ``Bond`` is modified,
#       rather than added/deleted via the API? The simplest resolution is simply to make them immutable.


class FrozenMolecule(Serializable):
    """
    Immutable chemical representation of a molecule, such as a small molecule or biopolymer.

    .. todo :: What other API calls would be useful for supporting biopolymers
               as small molecules? Perhaps iterating over chains and residues?

    Examples
    --------

    Create a molecule from a sdf file

    >>> from openff.toolkit.utils import get_data_file_path
    >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
    >>> molecule = FrozenMolecule.from_file(sdf_filepath)

    Convert to OpenEye OEMol object

    >>> oemol = molecule.to_openeye()

    Create a molecule from an OpenEye molecule

    >>> molecule = FrozenMolecule.from_openeye(oemol)

    Convert to RDKit Mol object

    >>> rdmol = molecule.to_rdkit()

    Create a molecule from an RDKit molecule

    >>> molecule = FrozenMolecule.from_rdkit(rdmol)

    Create a molecule from IUPAC name (requires the OpenEye toolkit)

    >>> molecule = FrozenMolecule.from_iupac('imatinib')

    Create a molecule from SMILES

    >>> molecule = FrozenMolecule.from_smiles('Cc1ccccc1')

    .. warning :: This API is experimental and subject to change.


    """

    def __init__(
        self,
        other=None,
        file_format=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
    ):
        r"""
        Create a new FrozenMolecule object

        .. todo ::

           * If a filename or file-like object is specified but the file
             contains more than one molecule, what is the proper behavior?
             Read just the first molecule, or raise an exception if more
             than one molecule is found?

           * Should we also support SMILES strings or IUPAC names for
             ``other``\ ?

        Parameters
        ----------
        other : optional, default=None
            If specified, attempt to construct a copy of the molecule from
            the specified object. This can be any one of the following:

            * a :class:`Molecule` object
            * a file that can be used to construct a :class:`Molecule` object
            * an ``openeye.oechem.OEMol``
            * an ``rdkit.Chem.rdchem.Mol``
            * a serialized :class:`Molecule` object

        file_format : str, optional, default=None
            If providing a file-like object, you must specify the format
            of the data. If providing a file, the file format will attempt
            to be guessed from the suffix.
        toolkit_registry : a :class:`ToolkitRegistry` or
            :class:`ToolkitWrapper` object, optional,
            default=GLOBAL_TOOLKIT_REGISTRY :class:`ToolkitRegistry`
            or :class:`ToolkitWrapper` to use for I/O operations
        allow_undefined_stereo : bool, default=False
            If loaded from a file and ``False``, raises an exception if
            undefined stereochemistry is detected during the molecule's
            construction.

        Examples
        --------

        Create an empty molecule:

        >>> empty_molecule = FrozenMolecule()

        Create a molecule from a file that can be used to construct a molecule,
        using either a filename or file-like object:

        >>> from openff.toolkit.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = FrozenMolecule(sdf_filepath)
        >>> molecule = FrozenMolecule(open(sdf_filepath, 'r'), file_format='sdf')

        >>> import gzip
        >>> mol2_gz_filepath = get_data_file_path('molecules/toluene.mol2.gz')
        >>> molecule = FrozenMolecule(gzip.GzipFile(mol2_gz_filepath, 'r'), file_format='mol2')

        Create a molecule from another molecule:

        >>> molecule_copy = FrozenMolecule(molecule)

        Convert to OpenEye OEMol object

        >>> oemol = molecule.to_openeye()

        Create a molecule from an OpenEye molecule:

        >>> molecule = FrozenMolecule(oemol)

        Convert to RDKit Mol object

        >>> rdmol = molecule.to_rdkit()

        Create a molecule from an RDKit molecule:

        >>> molecule = FrozenMolecule(rdmol)

        Convert the molecule into a dictionary and back again:

        >>> serialized_molecule = molecule.to_dict()
        >>> molecule_copy = Molecule(serialized_molecule)

        """

        self._cached_smiles = None

        # Figure out if toolkit_registry is a whole registry, or just a single wrapper
        if isinstance(toolkit_registry, ToolkitRegistry):
            pass
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            toolkit_registry = ToolkitRegistry(toolkit_precedence=[])
            toolkit_registry.add_toolkit(toolkit)
        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        if other is None:
            self._initialize()
        else:
            loaded = False
            # Start a list of the ValueErrors the following logic encounters, so we can print it out
            # if there turned out to be no way to load this input
            value_errors = list()

            if isinstance(other, openff.toolkit.topology.FrozenMolecule) and not loaded:
                self._copy_initializer(other)
                loaded = True
            if isinstance(other, openff.toolkit.topology.Molecule) and not loaded:
                # TODO: This will need to be updated once FrozenMolecules and Molecules are significantly different
                self._copy_initializer(other)
                loaded = True
            if isinstance(other, dict) and not loaded:
                self._initialize_from_dict(other)
                loaded = True

            # Check through the toolkit registry to find a compatible wrapper for loading
            if not loaded:
                try:
                    # Each ToolkitWrapper may provide a from_object method, which turns some particular type(s)
                    # of object into OFFMols. For example, RDKitToolkitWrapper's from_object method will
                    # return an OFFMol if provided with an RDMol, or raise a ValueError if it is provided
                    # an OEMol (or anything else). This makes the assumption that any non-ValueError errors raised
                    # by the toolkit _really are_ bad and should be raised immediately, which may be a bad assumption.
                    result = toolkit_registry.call(
                        "from_object",
                        other,
                        allow_undefined_stereo=allow_undefined_stereo,
                        raise_exception_types=[UndefinedStereochemistryError],
                        _cls=self.__class__,
                    )
                # NotImplementedError should never be raised... Only from_file and from_file_obj are provided
                # in the base ToolkitWrapper class and require overwriting, so from_object should be excluded
                # except NotImplementedError as e:
                #    raise e
                # The toolkit registry will aggregate all errors except UndefinedStereochemistryErrors into a single
                # ValueError, which we should catch and and store that here.
                except ValueError as e:
                    value_errors.append(e)
                else:
                    self._copy_initializer(result)
                    loaded = True
            # TODO: Make this compatible with file-like objects (I couldn't figure out how to make an oemolistream
            # from a fileIO object)
            if isinstance(other, str) or hasattr(other, "read") and not loaded:
                try:
                    mol = Molecule.from_file(
                        other,
                        file_format=file_format,
                        toolkit_registry=toolkit_registry,
                        allow_undefined_stereo=allow_undefined_stereo,
                    )  # returns a list only if multiple molecules are found
                    if type(mol) == list:
                        raise ValueError(
                            "Specified file or file-like object must contain exactly one molecule"
                        )
                except ValueError as e:
                    value_errors.append(e)
                else:
                    self._copy_initializer(mol)
                    loaded = True

            # If none of the above methods worked, raise a ValueError summarizing the
            # errors from the different loading attempts

            if not loaded:
                msg = (
                    f"Cannot construct openff.toolkit.topology.Molecule from {other}\n"
                )
                for value_error in value_errors:
                    msg += str(value_error)
                raise ValueError(msg)

    @property
    def has_unique_atom_names(self) -> bool:
        """True if the molecule has unique atom names, False otherwise."""
        unique_atom_names = set([atom.name for atom in self.atoms])
        if len(unique_atom_names) < self.n_atoms:
            return False
        return True

    def generate_unique_atom_names(self):
        """
        Generate unique atom names using element name and number of times that element has occurred
        e.g. 'C1x', 'H1x', 'O1x', 'C2x', ...

        The character 'x' is appended to these generated names to reduce the odds that they clash with an atom name or
        type imported from another source.

        """
        from collections import defaultdict

        element_counts = defaultdict(int)
        for atom in self.atoms:
            symbol = atom.symbol
            element_counts[symbol] += 1
            # TODO: It may be worth exposing this as a user option, i.e. to avoid multiple ligands
            # parameterized with OpenFF clashing because they have atom names like O1x, H3x, etc.
            # i.e. an optional argument could enable a user to `generate_unique_atom_names(blah="y")
            # to have one ligand be O1y, etc.
            # https://github.com/openforcefield/openff-toolkit/pull/1096#pullrequestreview-767227391
            atom.name = symbol + str(element_counts[symbol]) + "x"

    def _validate(self):
        """
        Validate the molecule, ensuring it has unique atom names

        """
        if not self.has_unique_atom_names:
            self.generate_unique_atom_names()

    def strip_atom_stereochemistry(
        self, smarts, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY
    ):
        """Delete stereochemistry information for certain atoms, if it is present.
        This method can be used to "normalize" molecules imported from different cheminformatics
        toolkits, which differ in which atom centers are considered stereogenic.

        Parameters
        ----------
        smarts: str or ChemicalEnvironment
            Tagged SMARTS with a single atom with index 1. Any matches for this atom will have any assigned
            stereocheistry information removed.
        toolkit_registry : a :class:`ToolkitRegistry` or :class:`ToolkitWrapper` object, optional,
            default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for I/O operations

        """
        from openff.toolkit.typing.chemistry.environment import AtomChemicalEnvironment

        chem_env = AtomChemicalEnvironment(smarts)
        matches = self.chemical_environment_matches(
            chem_env, toolkit_registry=toolkit_registry
        )

        for match in set(matches):
            atom_idx = match[0]
            self.atoms[atom_idx].stereochemistry = None

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
        from openff.toolkit.utils.utils import serialize_numpy

        molecule_dict = OrderedDict()
        molecule_dict["name"] = self._name
        # From Jeff: If we go the properties-as-dict route, then _properties should, at
        # the top level, be a dict. Should we go through recursively and ensure all values are dicts too?
        molecule_dict["atoms"] = [atom.to_dict() for atom in self._atoms]
        molecule_dict["bonds"] = [bond.to_dict() for bond in self._bonds]
        # TODO: Charges
        # TODO: Properties
        # From Jeff: We could have the onerous requirement that all "properties" have to_dict() functions.
        # Or we could restrict properties to simple stuff (ints, strings, floats, and the like)
        # Or pickle anything unusual
        # Or not allow user-defined properties at all (just use our internal _cached_properties)
        # molecule_dict['properties'] = dict([(key, value._to_dict()) for key.value in self._properties])
        # TODO: Assuming "simple stuff" properties right now, figure out a better standard
        molecule_dict["properties"] = self._properties
        if hasattr(self, "_cached_properties"):
            molecule_dict["cached_properties"] = self._cached_properties
        # TODO: Conformers
        if self._conformers is None:
            molecule_dict["conformers"] = None
        else:
            molecule_dict["conformers"] = []
            molecule_dict[
                "conformers_unit"
            ] = "angstrom"  # Have this defined as a class variable?
            for conf in self._conformers:
                conf_unitless = conf.m_as(unit.angstrom)
                conf_serialized, conf_shape = serialize_numpy((conf_unitless))
                molecule_dict["conformers"].append(conf_serialized)
        if self._partial_charges is None:
            molecule_dict["partial_charges"] = None
            molecule_dict["partial_charges_unit"] = None

        else:
            charges_unitless = self._partial_charges.m_as(unit.elementary_charge)
            charges_serialized, charges_shape = serialize_numpy(charges_unitless)
            molecule_dict["partial_charges"] = charges_serialized
            molecule_dict["partial_charges_unit"] = "elementary_charge"

        molecule_dict["hierarchy_schemes"] = dict()
        for iter_name, hier_scheme in self._hierarchy_schemes.items():
            molecule_dict["hierarchy_schemes"][iter_name] = hier_scheme.to_dict()

        return molecule_dict

    def __hash__(self):
        """
        Returns a hash of this molecule. Used when checking molecule uniqueness in Topology creation.

        Returns
        -------
        string
        """
        return hash(self.to_smiles())

    # @cached_property
    def ordered_connection_table_hash(self):
        """Compute an ordered hash of the atoms and bonds in the molecule"""
        if self._ordered_connection_table_hash is not None:
            return self._ordered_connection_table_hash

        id = ""
        for atom in self.atoms:
            id += f"{atom.symbol}_{atom.formal_charge}_{atom.stereochemistry}__"
        for bond in self.bonds:
            id += f"{bond.bond_order}_{bond.stereochemistry}_{bond.atom1_index}_{bond.atom2_index}__"
        # return hash(id)
        self._ordered_connection_table_hash = hash(id)
        return self._ordered_connection_table_hash

    @classmethod
    def from_dict(cls, molecule_dict):
        """
        Create a new Molecule from a dictionary representation

        Parameters
        ----------
        molecule_dict : OrderedDict
            A dictionary representation of the molecule.

        Returns
        -------
        molecule : Molecule
            A Molecule created from the dictionary representation

        """
        # This implementation is a compromise to let this remain as a classmethod
        mol = cls()
        mol._initialize_from_dict(molecule_dict)
        return mol

    def _initialize_from_dict(self, molecule_dict):
        """
        Initialize the molecule from a dictionary representation

        Parameters
        ----------
        molecule_dict : OrderedDict
            A dictionary representation of the molecule.
        """
        # TODO: Provide useful exception messages if there are any failures
        from openff.toolkit.utils.utils import deserialize_numpy

        self._initialize()
        self.name = molecule_dict["name"]
        for atom_dict in molecule_dict["atoms"]:
            self._add_atom(**atom_dict)

        for bond_dict in molecule_dict["bonds"]:
            bond_dict["atom1"] = int(bond_dict["atom1"])
            bond_dict["atom2"] = int(bond_dict["atom2"])
            self._add_bond(**bond_dict)

        if molecule_dict["partial_charges"] is None:
            self._partial_charges = None
        else:
            charges_shape = (self.n_atoms,)
            partial_charges_unitless = deserialize_numpy(
                molecule_dict["partial_charges"], charges_shape
            )
            pc_unit = getattr(unit, molecule_dict["partial_charges_unit"])
            partial_charges = unit.Quantity(partial_charges_unitless, pc_unit)
            self._partial_charges = partial_charges

        if molecule_dict["conformers"] is None:
            self._conformers = None
        else:
            self._conformers = list()
            for ser_conf in molecule_dict["conformers"]:
                # TODO: Update to use string_to_quantity
                conformers_shape = (self.n_atoms, 3)
                conformer_unitless = deserialize_numpy(ser_conf, conformers_shape)
                c_unit = getattr(unit, molecule_dict["conformers_unit"])
                conformer = unit.Quantity(conformer_unitless, c_unit)
                self._conformers.append(conformer)

        self._properties = molecule_dict["properties"]

        for iter_name, hierarchy_scheme_dict in molecule_dict[
            "hierarchy_schemes"
        ].items():
            # It's important that we do NOT call `add_hierarchy_scheme` here, since we
            # need to deserialize these HierarchyElements exactly as they were serialized,
            # even if that conflicts with the current values in atom metadata.
            new_hier_scheme = HierarchyScheme(
                self,
                hierarchy_scheme_dict["uniqueness_criteria"],
                iter_name,
            )
            self._hierarchy_schemes[iter_name] = new_hier_scheme

            for element_dict in hierarchy_scheme_dict["hierarchy_elements"]:
                new_hier_scheme.add_hierarchy_element(
                    element_dict["identifier"], element_dict["atom_indices"]
                )

    def __repr__(self):
        """Return a summary of this molecule; SMILES if valid, Hill formula if not."""
        description = f"Molecule with name '{self.name}'"
        try:
            smiles = self.to_smiles()
        except Exception:
            hill = self.to_hill_formula()
            return description + f" with bad SMILES and Hill formula '{hill}'"
        return description + f" and SMILES '{smiles}'"

    def _initialize(self):
        """
        Clear the contents of the current molecule.
        """
        self._name = ""
        self._atoms = list()
        self._bonds = list()  # List of bonds between Atom objects
        self._properties = {}  # Attached properties to be preserved
        # self._cached_properties = None # Cached properties (such as partial charges) can be recomputed as needed
        self._partial_charges = None
        self._conformers = None  # Optional conformers
        self._hierarchy_schemes = dict()
        self._invalidate_cached_properties()

    def _copy_initializer(self, other):
        """
        Copy contents of the specified molecule

        .. todo :: Should this be a ``@staticmethod`` where we have an explicit copy constructor?

        Parameters
        ----------
        other : optional
            Overwrite the state of this FrozenMolecule with the specified FrozenMolecule object.
            A deep copy is made.

        """
        # assert isinstance(other, type(self)), "can only copy instances of {}".format(type(self))

        # Run a deepcopy here so that items that were _always_ dict (like other.properties) will
        # not have any references to the old molecule
        other_dict = deepcopy(other.to_dict())
        self._initialize_from_dict(other_dict)

    def __eq__(self, other):
        """
        Test two molecules for equality to see if they are the chemical species, but do not check other
        annotated properties.

        .. note ::

           Note that this method simply tests whether two molecules are identical chemical species using equivalence of
           their canonical isomeric SMILES.  No effort is made to ensure that the atoms are in the same order or that
           any annotated properties are preserved.

        """
        # updated to use the new isomorphic checking method, with full matching
        # TODO the doc string did not match the previous function what matching should this method do?
        return Molecule.are_isomorphic(self, other, return_atom_map=False)[0]

    def __deepcopy__(self, memo):
        cls = self.__class__
        return cls(self.to_dict())

    def add_default_hierarchy_schemes(self, overwrite_existing=True):
        """
        Adds ``chain`` and ``residue`` hierarchy schemes.

        The Open Force Field Toolkit has no native understanding of hierarchical
        atom organisation schemes common to other biomolecular software, such as
        "residues" or "chains" (see :ref:`userguide_hierarchy`). Hierarchy
        schemes allow iteration over groups of atoms according to their
        metadata. For more information, see
        :class:`~openff.toolkit.topology.molecule.HierarchyScheme`.

        If a ``Molecule`` with the default hierarchy schemes
        changes, :meth:`Molecule.update_hierarchy_schemes()` must be called before
        the residues or chains are iterated over again or else the iteration may
        be incorrect.

        Parameters
        ----------
        overwrite_existing : bool, default=True
            Whether to overwrite existing instances of the `residue` and `chain`
            hierarchy schemes. If this is ``False`` and either of the hierarchy
            schemes are already defined on this molecule, an exception will be
            raised.

        Raises
        ------
        HierarchySchemeWithIteratorNameAlreadyRegisteredException
            When ``overwrite_existing=False`` and either the ``chains`` or
            ``residues`` hierarchy scheme is already configured.
        """
        self._add_chain_hierarchy_scheme(overwrite_existing=overwrite_existing)
        self._add_residue_hierarchy_scheme(overwrite_existing=overwrite_existing)

    def _add_chain_hierarchy_scheme(self, overwrite_existing=True):
        """Add ``chain`` hierarchy scheme."""
        if overwrite_existing:
            if "chains" in self._hierarchy_schemes.keys():
                self.delete_hierarchy_scheme("chains")

        self.add_hierarchy_scheme(("chain_id",), "chains")

    def _add_residue_hierarchy_scheme(self, overwrite_existing=True):
        """Add ``residue`` hierarchy scheme."""
        if overwrite_existing:
            if "residues" in self._hierarchy_schemes.keys():
                self.delete_hierarchy_scheme("residues")

        self.add_hierarchy_scheme(
            ("chain_id", "residue_number", "residue_name"), "residues"
        )

    def add_hierarchy_scheme(
        self,
        uniqueness_criteria,
        iterator_name,
    ):
        """
        Use the molecule's metadata to facilitate iteration over its atoms.

        This method will add an attribute with the name given by the
        ``iterator_name`` argument that provides an iterator over groups of
        atoms. Atoms are grouped by the values in their ``atom.metadata``
        dictionary; any atoms with the same values for the keys given in the
        ``uniqueness_criteria`` argument will be in the same group. These groups
        have the type :class:`~openff.toolkit.topology.molecule.HierarchyElement`.

        Hierarchy schemes are not updated dynamically; if a ``Molecule`` with
        hierarchy schemes changes, :meth:`Molecule.update_hierarchy_schemes()` must
        be called before the scheme is iterated over again or else the grouping
        may be incorrect.

        Hierarchy schemes allow iteration over groups of atoms according to
        their metadata. For more information, see
        :class:`~openff.toolkit.topology.molecule.HierarchyScheme`.

        Parameters
        ----------
        uniqueness_criteria : tuple of str
            The names of ``Atom`` metadata entries that define this scheme. An
            atom belongs to a ``HierarchyElement`` only if its metadata has the
            same values for these criteria as the other atoms in the
            ``HierarchyElement``.

        iterator_name : str
            Name of the iterator that will be exposed to access the hierarchy
            elements generated by this scheme.

        Returns
        -------
        new_hier_scheme : openff.toolkit.topology.HierarchyScheme
            The newly created HierarchyScheme

        """
        if iterator_name in self._hierarchy_schemes:
            msg = (
                f'Can not add iterator with name "{iterator_name}" to this topology, as iterator '
                f"name is already used by {self._hierarchy_schemes[iterator_name]}"
            )
            raise HierarchySchemeWithIteratorNameAlreadyRegisteredException(msg)
        new_hier_scheme = HierarchyScheme(
            self,
            uniqueness_criteria,
            iterator_name,
        )
        self._hierarchy_schemes[iterator_name] = new_hier_scheme
        self.update_hierarchy_schemes([iterator_name])
        return new_hier_scheme

    @property
    def hierarchy_schemes(self) -> Dict[str, "HierarchyScheme"]:
        """
        The hierarchy schemes available on the molecule.

        Hierarchy schemes allow iteration over groups of atoms according to
        their metadata. For more information, see
        :class:`~openff.toolkit.topology.molecule.HierarchyScheme`.

        Returns
        -------
        A dict of the form {str: HierarchyScheme}
            The HierarchySchemes associated with the molecule.
        """
        return self._hierarchy_schemes

    def delete_hierarchy_scheme(self, iter_name):
        """
        Remove an existing ``HierarchyScheme`` specified by its iterator name.

        Hierarchy schemes allow iteration over groups of atoms according to
        their metadata. For more information, see
        :class:`~openff.toolkit.topology.molecule.HierarchyScheme`.

        Parameters
        ----------
        iter_name : str
        """
        if iter_name not in self._hierarchy_schemes:
            raise HierarchySchemeNotFoundException(
                f'Can not delete HierarchyScheme with name "{iter_name}" '
                f"because no HierarchyScheme with that iterator name exists"
            )
        self._hierarchy_schemes.pop(iter_name)

    def update_hierarchy_schemes(self, iter_names=None):
        """
        Infer a hierarchy from atom metadata according to the existing hierarchy
        schemes.

        Hierarchy schemes allow iteration over groups of atoms according to
        their metadata. For more information, see
        :class:`~openff.toolkit.topology.molecule.HierarchyScheme`.

        Parameters
        ----------
        iter_names : Iterable of str, Optional
            Only perceive hierarchy for HierarchySchemes that expose these
            iterator names. If not provided, all known hierarchies will be
            perceived, overwriting previous results if applicable.
        """
        if iter_names is None:
            iter_names = self._hierarchy_schemes.keys()

        for iter_name in iter_names:
            hierarchy_scheme = self._hierarchy_schemes[iter_name]
            hierarchy_scheme.perceive_hierarchy()

    def __getattr__(self, name: str):
        """If a requested attribute is not found, check the hierarchy schemes"""
        try:
            return self.__dict__["_hierarchy_schemes"][name].hierarchy_elements
        except KeyError:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute {name!r}"
            )

    def __dir__(self):
        """Add the hierarchy scheme iterator names to dir"""
        return list(self._hierarchy_schemes.keys()) + list(super().__dir__())

    def to_smiles(
        self,
        isomeric=True,
        explicit_hydrogens=True,
        mapped=False,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Return a canonical isomeric SMILES representation of the current molecule.
        A partially mapped smiles can also be generated for atoms of interest by supplying an `atom_map` to the
        properties dictionary.

        .. note :: RDKit and OpenEye versions will not necessarily return the same representation.

        Parameters
        ----------
        isomeric: bool optional, default= True
            return an isomeric smiles
        explicit_hydrogens: bool optional, default=True
            return a smiles string containing all hydrogens explicitly
        mapped: bool optional, default=False
            return a explicit hydrogen mapped smiles, the atoms to be mapped can be controlled by supplying an
            atom map into the properties dictionary. If no mapping is passed all atoms will be mapped in order, else
            an atom map dictionary from the current atom index to the map id should be supplied with no duplicates.
            The map ids (values) should start from 0 or 1.
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or
            openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES conversion

        Returns
        -------
        smiles : str
            Canonical isomeric explicit-hydrogen SMILES

        Examples
        --------

        >>> from openff.toolkit.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> smiles = molecule.to_smiles()

        """
        # Initialize cached_smiles dict for this molecule if none exists
        if self._cached_smiles is None:
            self._cached_smiles = {}

        # Figure out which toolkit should be used to create the SMILES
        if isinstance(toolkit_registry, ToolkitRegistry):
            to_smiles_method = toolkit_registry.resolve("to_smiles")
        elif isinstance(toolkit_registry, ToolkitWrapper):
            to_smiles_method = toolkit_registry.to_smiles
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to to_smiles. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        # Get a string representation of the function containing the toolkit name so we can check
        # if a SMILES was already cached for this molecule. This will return, for example
        # "RDKitToolkitWrapper.to_smiles"
        smiles_hash = (
            to_smiles_method.__qualname__
            + str(isomeric)
            + str(explicit_hydrogens)
            + str(mapped)
        )
        smiles_hash += str(self._properties.get("atom_map", None))
        # Check to see if a SMILES for this molecule was already cached using this method
        if smiles_hash in self._cached_smiles:
            return self._cached_smiles[smiles_hash]
        else:
            smiles = to_smiles_method(self, isomeric, explicit_hydrogens, mapped)
            self._cached_smiles[smiles_hash] = smiles
            return smiles

    @classmethod
    def from_inchi(
        cls,
        inchi,
        allow_undefined_stereo=False,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Construct a Molecule from a InChI representation

        Parameters
        ----------
        inchi : str
            The InChI representation of the molecule.

        allow_undefined_stereo : bool, default=False
            Whether to accept InChI with undefined stereochemistry. If False,
            an exception will be raised if a InChI with undefined stereochemistry
            is passed into this function.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for InChI-to-molecule conversion


        Returns
        -------
        molecule : openff.toolkit.topology.Molecule

        Examples
        --------
        Make cis-1,2-Dichloroethene:

        >>> molecule = Molecule.from_inchi('InChI=1S/C2H2Cl2/c3-1-2-4/h1-2H/b2-1-')
        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            molecule = toolkit_registry.call(
                "from_inchi",
                inchi,
                _cls=cls,
                allow_undefined_stereo=allow_undefined_stereo,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            molecule = toolkit.from_inchi(
                inchi, _cls=cls, allow_undefined_stereo=allow_undefined_stereo
            )
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to from_inchi. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        return molecule

    def to_inchi(self, fixed_hydrogens=False, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Create an InChI string for the molecule using the requested toolkit backend.
        InChI is a standardised representation that does not capture tautomers unless specified using the fixed
        hydrogen layer.

        For information on InChi see here https://iupac.org/who-we-are/divisions/division-details/inchi/

        Parameters
        ----------
        fixed_hydrogens: bool, default=False
            If a fixed hydrogen layer should be added to the InChI, if `True` this will produce a non standard
            specific InChI string of the molecule.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for molecule-to-InChI conversion

        Returns
        --------
        inchi: str
            The InChI string of the molecule.

        Raises
        -------
        InvalidToolkitRegistryError
             If an invalid object is passed as the toolkit_registry parameter
        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            inchi = toolkit_registry.call(
                "to_inchi", self, fixed_hydrogens=fixed_hydrogens
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            inchi = toolkit.to_inchi(self, fixed_hydrogens=fixed_hydrogens)
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to to_inchi. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        return inchi

    def to_inchikey(
        self, fixed_hydrogens=False, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY
    ):
        """
        Create an InChIKey for the molecule using the requested toolkit backend.
        InChIKey is a standardised representation that does not capture tautomers unless specified
        using the fixed hydrogen layer.

        For information on InChi see here https://iupac.org/who-we-are/divisions/division-details/inchi/

        Parameters
        ----------
        fixed_hydrogens: bool, default=False
            If a fixed hydrogen layer should be added to the InChI, if `True` this will produce a non standard specific
            InChI string of the molecule.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for molecule-to-InChIKey conversion

        Returns
        --------
        inchi_key: str
            The InChIKey representation of the molecule.

        Raises
        -------
        InvalidToolkitRegistryError
             If an invalid object is passed as the toolkit_registry parameter
        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            inchi_key = toolkit_registry.call(
                "to_inchikey", self, fixed_hydrogens=fixed_hydrogens
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            inchi_key = toolkit.to_inchikey(self, fixed_hydrogens=fixed_hydrogens)
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to to_inchikey. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        return inchi_key

    @classmethod
    def from_smiles(
        cls,
        smiles,
        hydrogens_are_explicit=False,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
    ):
        """
        Construct a Molecule from a SMILES representation

        Parameters
        ----------
        smiles : str
            The SMILES representation of the molecule.
        hydrogens_are_explicit : bool, default = False
            If False, the cheminformatics toolkit will perform hydrogen addition
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion
        allow_undefined_stereo : bool, default=False
            Whether to accept SMILES with undefined stereochemistry. If False,
            an exception will be raised if a SMILES with undefined stereochemistry
            is passed into this function.

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule

        Examples
        --------

        >>> molecule = Molecule.from_smiles('Cc1ccccc1')

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            molecule = toolkit_registry.call(
                "from_smiles",
                smiles,
                hydrogens_are_explicit=hydrogens_are_explicit,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            molecule = toolkit.from_smiles(
                smiles,
                hydrogens_are_explicit=hydrogens_are_explicit,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
            )
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to from_smiles. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        return molecule

    def _is_exactly_the_same_as(self, other):
        for atom1, atom2 in zip(self.atoms, other.atoms):
            if (
                (atom1.atomic_number != atom2.atomic_number)
                or (atom1.formal_charge != atom2.formal_charge)
                or (atom1.is_aromatic != atom2.is_aromatic)
                or (atom1.stereochemistry != atom2.stereochemistry)
            ):
                return False
        for bond1, bond2 in zip(self.bonds, other.bonds):
            if (
                (bond1.atom1_index != bond2.atom1_index)
                or (bond1.atom2_index != bond2.atom2_index)
                or (bond1.is_aromatic != bond2.is_aromatic)
                or (bond1.stereochemistry != bond2.stereochemistry)
            ):
                return False
        return True

    @staticmethod
    def are_isomorphic(
        mol1,
        mol2,
        return_atom_map=False,
        aromatic_matching=True,
        formal_charge_matching=True,
        bond_order_matching=True,
        atom_stereochemistry_matching=True,
        bond_stereochemistry_matching=True,
        strip_pyrimidal_n_atom_stereo=True,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Determine if ``mol1`` is isomorphic to ``mol2``.

        ``are_isomorphic()`` compares two molecule's graph representations and
        the chosen node/edge attributes. Connections and atomic numbers are
        always checked.

        If nx.Graphs() are given they must at least have ``atomic_number``
        attributes on nodes. Other attributes that ``are_isomorphic()`` can
        optionally check...

        -  ... in nodes are:

           -  ``is_aromatic``
           -  ``formal_charge``
           -  ``stereochemistry``

        -  ... in edges are:

           -  ``is_aromatic``
           -  ``bond_order``
           -  ``stereochemistry``

        By default, all attributes are checked, but stereochemistry around
        pyrimidal nitrogen is ignored.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        mol1 : an openff.toolkit.topology.molecule.FrozenMolecule or nx.Graph()
            The first molecule to test for isomorphism.

        mol2 : an openff.toolkit.topology.molecule.FrozenMolecule or nx.Graph()
            The second molecule to test for isomorphism.

        return_atom_map: bool, default=False, optional
            Return a ``dict`` containing the atomic mapping instead of a
            ``bool``.

        aromatic_matching: bool, default=True, optional
            If ``False``, aromaticity of graph nodes and edges are ignored for
            the purpose of determining isomorphism.

        formal_charge_matching: bool, default=True, optional
            If ``False``, formal charges of graph nodes are ignored for
            the purpose of determining isomorphism.

        bond_order_matching: bool, default=True, optional
            If ``False``, bond orders of graph edges are ignored for
            the purpose of determining isomorphism.

        atom_stereochemistry_matching : bool, default=True, optional
            If ``False``, atoms' stereochemistry is ignored for the
            purpose of determining isomorphism.

        bond_stereochemistry_matching : bool, default=True, optional
            If ``False``, bonds' stereochemistry is ignored for the
            purpose of determining isomorphism.

        strip_pyrimidal_n_atom_stereo: bool, default=True, optional
            If ``True``, any stereochemistry defined around pyrimidal
            nitrogen stereocenters will be disregarded in the isomorphism
            check.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for
            removing stereochemistry from pyrimidal nitrogens.

        Returns
        -------
        molecules_are_isomorphic : bool

        atom_map : default=None, Optional,
            [Dict[int,int]] ordered by mol1 indexing {mol1_index: mol2_index}
            If molecules are not isomorphic given input arguments, will return None instead of dict.
        """
        # Do a quick hill formula check first
        if Molecule._object_to_hill_formula(mol1) != Molecule._object_to_hill_formula(
            mol2
        ):
            return False, None

        # Do a quick check to see whether the inputs are totally identical (including being in the same atom order)
        if isinstance(mol1, FrozenMolecule) and isinstance(mol2, FrozenMolecule):
            if mol1._is_exactly_the_same_as(mol2):
                return True, {i: i for i in range(mol1.n_atoms)}

        # Build the user defined matching functions
        def node_match_func(x, y):
            # always match by atleast atomic number
            is_equal = x["atomic_number"] == y["atomic_number"]
            if aromatic_matching:
                is_equal &= x["is_aromatic"] == y["is_aromatic"]
            if formal_charge_matching:
                is_equal &= x["formal_charge"] == y["formal_charge"]
            if atom_stereochemistry_matching:
                is_equal &= x["stereochemistry"] == y["stereochemistry"]
            return is_equal

        # check if we want to do any bond matching if not the function is None
        if aromatic_matching or bond_order_matching or bond_stereochemistry_matching:

            def edge_match_func(x, y):
                # We don't need to check the exact bond order (which is 1 or 2)
                # if the bond is aromatic. This way we avoid missing a match only
                # if the alternate bond orders 1 and 2 are assigned differently.
                if aromatic_matching and bond_order_matching:
                    is_equal = (x["is_aromatic"] == y["is_aromatic"]) or (
                        x["bond_order"] == y["bond_order"]
                    )
                elif aromatic_matching:
                    is_equal = x["is_aromatic"] == y["is_aromatic"]
                elif bond_order_matching:
                    is_equal = x["bond_order"] == y["bond_order"]
                else:
                    is_equal = None
                if bond_stereochemistry_matching:
                    if is_equal is None:
                        is_equal = x["stereochemistry"] == y["stereochemistry"]
                    else:
                        is_equal &= x["stereochemistry"] == y["stereochemistry"]

                return is_equal

        else:
            edge_match_func = None

        # Here we should work out what data type we have, also deal with lists?
        def to_networkx(data):
            """For the given data type, return the networkx graph"""
            import networkx as nx

            if strip_pyrimidal_n_atom_stereo:
                SMARTS = "[N+0X3:1](-[*])(-[*])(-[*])"

            if isinstance(data, FrozenMolecule):
                # Molecule class instance
                if strip_pyrimidal_n_atom_stereo:
                    # Make a copy of the molecule so we don't modify the original
                    data = deepcopy(data)
                    data.strip_atom_stereochemistry(
                        SMARTS, toolkit_registry=toolkit_registry
                    )
                return data.to_networkx()

            elif isinstance(data, nx.Graph):
                return data

            else:
                raise NotImplementedError(
                    f"The input type {type(data)} is not supported,"
                    f"please supply an openff.toolkit.topology.molecule.Molecule "
                    f"or networkx.Graph representation of the molecule."
                )

        mol1_netx = to_networkx(mol1)
        mol2_netx = to_networkx(mol2)
        from networkx.algorithms.isomorphism import GraphMatcher  # type: ignore

        GM = GraphMatcher(
            mol1_netx, mol2_netx, node_match=node_match_func, edge_match=edge_match_func
        )
        isomorphic = GM.is_isomorphic()

        if isomorphic and return_atom_map:
            topology_atom_map = GM.mapping

            # reorder the mapping by keys
            sorted_mapping = {}
            for key in sorted(topology_atom_map.keys()):
                sorted_mapping[key] = topology_atom_map[key]

            return isomorphic, sorted_mapping

        else:
            return isomorphic, None

    def is_isomorphic_with(self, other, **kwargs):
        """
        Check if the molecule is isomorphic with the other molecule which can be an openff.toolkit.topology.Molecule
        or nx.Graph(). Full matching is done using the options described bellow.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        other: openff.toolkit.topology.Molecule or nx.Graph()

        aromatic_matching: bool, default=True, optional
        compare the aromatic attributes of bonds and atoms.

        formal_charge_matching: bool, default=True, optional
        compare the formal charges attributes of the atoms.

        bond_order_matching: bool, deafult=True, optional
        compare the bond order on attributes of the bonds.

        atom_stereochemistry_matching : bool, default=True, optional
            If ``False``, atoms' stereochemistry is ignored for the
            purpose of determining equality.

        bond_stereochemistry_matching : bool, default=True, optional
            If ``False``, bonds' stereochemistry is ignored for the
            purpose of determining equality.

        strip_pyrimidal_n_atom_stereo: bool, default=True, optional
            If ``True``, any stereochemistry defined around pyrimidal
            nitrogen stereocenters will be disregarded in the isomorphism
            check.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for
            removing stereochemistry from pyrimidal nitrogens.

        Returns
        -------
        isomorphic : bool
        """

        return Molecule.are_isomorphic(
            self,
            other,
            return_atom_map=False,
            aromatic_matching=kwargs.get("aromatic_matching", True),
            formal_charge_matching=kwargs.get("formal_charge_matching", True),
            bond_order_matching=kwargs.get("bond_order_matching", True),
            atom_stereochemistry_matching=kwargs.get(
                "atom_stereochemistry_matching", True
            ),
            bond_stereochemistry_matching=kwargs.get(
                "bond_stereochemistry_matching", True
            ),
            strip_pyrimidal_n_atom_stereo=kwargs.get(
                "strip_pyrimidal_n_atom_stereo", True
            ),
            toolkit_registry=kwargs.get("toolkit_registry", GLOBAL_TOOLKIT_REGISTRY),
        )[0]

    def generate_conformers(
        self,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        n_conformers=10,
        rms_cutoff=None,
        clear_existing=True,
        make_carboxylic_acids_cis=True,
    ):
        """
        Generate conformers for this molecule using an underlying toolkit.

        If ``n_conformers=0``, no toolkit wrapper will be called. If ``n_conformers=0``
        and ``clear_existing=True``, ``molecule.conformers`` will be set to ``None``.

        Parameters
        ----------
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or
            openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion
        n_conformers : int, default=1
            The maximum number of conformers to produce
        rms_cutoff : openmm.unit.Quantity-wrapped float, in units of distance, optional, default=None
            The minimum RMS value at which two conformers are considered redundant and one is deleted. Precise
            implementation of this cutoff may be toolkit-dependent. If ``None``, the cutoff is set to be the
            default value for each ``ToolkitWrapper`` (generally 1 Angstrom).
        clear_existing : bool, default=True
            Whether to overwrite existing conformers for the molecule
        make_carboxylic_acids_cis: bool, default=True
            Guarantee all conformers have exclusively cis carboxylic acid groups (COOH)
            by rotating the proton in any trans carboxylic acids 180 degrees around the
            C-O bond. Works around a bug in conformer generation by the OpenEye toolkit
            where trans COOH is much more common than it should be.

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.generate_conformers()

        Raises
        ------
        InvalidToolkitRegistryError
            If an invalid object is passed as the toolkit_registry parameter

        """
        # If no conformers are requested, do not call to a ToolkitWrapper at all
        if n_conformers == 0:
            if clear_existing:
                self._conformers = None
            return

        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call(
                "generate_conformers",
                self,
                n_conformers=n_conformers,
                rms_cutoff=rms_cutoff,
                clear_existing=clear_existing,
                raise_exception_types=[],
                make_carboxylic_acids_cis=make_carboxylic_acids_cis,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            return toolkit.generate_conformers(
                self,
                n_conformers=n_conformers,
                rms_cutoff=rms_cutoff,
                clear_existing=clear_existing,
                make_carboxylic_acids_cis=make_carboxylic_acids_cis,
            )
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to generate_conformers. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

    def _make_carboxylic_acids_cis(self, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Rotate dihedral angle of any conformers with trans COOH groups so they are cis

        Carboxylic acid groups almost always exist in nature in the cis conformation,
        with the hydrogen atom in between the two oxygen atoms::

                  O----H
                 /
                /
               /
            --C
               \\
                \\
                  O

        However, the OpenEye toolkit frequently produces carboxylic acid geometries
        in the unrealistic trans conformation::

             H----O
                 /
                /
               /
            --C
               \\
                \\
                  O

        This method converts all conformers in the molecule with the trans conformation
        into the corresponding cis conformer by rotating the OH bond around the CO bond
        by 180 degrees. Carboxylic acids that are already cis are unchanged. Carboxylic
        acid groups are considered cis if their O-C-O-H dihedral angle is acute.

        Parameters
        ----------
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion
        """

        # Return early if there are no conformers
        if not self._conformers:
            return

        # Convert all conformers into one big array
        conformers = np.asarray([q.m_as(unit.angstrom) for q in self._conformers])

        # Scan the molecule for carboxylic acids
        cooh_indices = self.chemical_environment_matches(
            "[C:2]([O:3][H:4])=[O:1]", toolkit_registry=toolkit_registry
        )
        n_conformers, n_cooh_groups = len(conformers), len(cooh_indices)
        # Exit early if there are no carboxylic acids
        if not n_cooh_groups:
            return

        # Pull out the coordinates of all carboxylic acid groups into cooh_xyz
        cooh_xyz = conformers[:, cooh_indices, :]
        assert cooh_xyz.shape == (n_conformers, n_cooh_groups, 4, 3)

        def dot(a, b):
            """Compute dot product along last axis of arrays"""
            return np.sum(a * b, axis=-1)[..., np.newaxis]

        def norm(a):
            """Compute norm along last axis of array"""
            return np.linalg.norm(a, axis=-1)[..., np.newaxis]

        def dihedral(a):
            """Compute dihedrals of array with shape (..., 4, 3)"""
            # Praxeolitic formula
            # 1 sqrt, 1 cross product
            # from https://stackoverflow.com/q/20305272
            p0 = a[..., 0, :]
            p1 = a[..., 1, :]
            p2 = a[..., 2, :]
            p3 = a[..., 3, :]

            b0 = -1.0 * (p1 - p0)
            b1 = p2 - p1
            b2 = p3 - p2

            # normalize b1 so that it does not influence magnitude of vector
            # rejections that come next
            b1 /= norm(b1)

            # vector rejections
            # v = projection of b0 onto plane perpendicular to b1
            #   = b0 minus component that aligns with b1
            # w = projection of b2 onto plane perpendicular to b1
            #   = b2 minus component that aligns with b1
            v = b0 - dot(b0, b1) * b1
            w = b2 - dot(b2, b1) * b1

            # angle between v and w in a plane is the torsion angle
            # v and w may not be normalized but that's fine since tan is y/x
            x = dot(v, w)
            y = dot(np.cross(b1, v), w)
            return np.arctan2(y, x)

        dihedrals = dihedral(cooh_xyz)
        assert dihedrals.shape == (n_conformers, n_cooh_groups, 1)
        dihedrals.shape = (n_conformers, n_cooh_groups, 1, 1)

        # Get indices of trans COOH groups
        trans_indices = np.logical_not(
            np.logical_and((-np.pi / 2) < dihedrals, dihedrals < (np.pi / 2))
        )
        # Expand array so it can be used to index cooh_xyz
        trans_indices = np.repeat(trans_indices, repeats=4, axis=2)
        trans_indices = np.repeat(trans_indices, repeats=3, axis=3)
        # Get indices of individual atoms in trans COOH groups (except terminal O)
        trans_indices_h = trans_indices.copy()
        trans_indices_h[:, :, (0, 1, 2), :] = False
        trans_indices_c = trans_indices.copy()
        trans_indices_c[:, :, (0, 2, 3), :] = False
        trans_indices_o = trans_indices.copy()
        trans_indices_o[:, :, (0, 1, 3), :] = False

        # Rotate OH around CO bond
        # We want to rotate H 180 degrees around the CO bond (b1)
        c = cooh_xyz[trans_indices_c].reshape(-1, 3)
        o = cooh_xyz[trans_indices_o].reshape(-1, 3)
        h = cooh_xyz[trans_indices_h].reshape(-1, 3)
        # Axis is defined as the line from the origin along a unit vector, so
        # move C to the origin and normalize
        point = h - c
        axis = o - c
        axis /= norm(axis)
        # Do the rotation
        # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
        rotated = axis * (dot(axis, point)) - np.cross(np.cross(axis, point), axis)
        # Move rotated point back to original coordinates
        rotated = rotated + c

        # Update the coordinates
        cooh_xyz[trans_indices_h] = rotated.reshape((-1))

        # Update conformers with rotated coordinates
        conformers[:, cooh_indices, :] = cooh_xyz

        # Return conformers to original type
        self._conformers = [unit.Quantity(conf, unit.angstrom) for conf in conformers]

    def apply_elf_conformer_selection(
        self,
        percentage: float = 2.0,
        limit: int = 10,
        toolkit_registry: Optional[
            Union[ToolkitRegistry, ToolkitWrapper]
        ] = GLOBAL_TOOLKIT_REGISTRY,
        **kwargs,
    ):
        """Select a set of diverse conformers from the molecule's conformers with ELF.

        Applies the `Electrostatically Least-interacting Functional groups method
        <https://docs.eyesopen.com/toolkits/python/quacpactk/molchargetheory.html#elf-conformer-selection>`_
        to select a set of diverse conformers which have minimal
        electrostatically strongly interacting functional groups from the
        molecule's conformers.

        Parameters
        ----------
        toolkit_registry
            The underlying toolkit to use to select the ELF conformers.
        percentage
            The percentage of conformers with the lowest electrostatic
            interaction energies to greedily select from.
        limit
            The maximum number of conformers to select.

        Notes
        -----
        * The input molecule should have a large set of conformers already
          generated to select the ELF conformers from.
        * The selected conformers will be retained in the `conformers` list
          while unselected conformers will be discarded.

        See Also
        --------
        openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.apply_elf_conformer_selection
        openff.toolkit.utils.toolkits.RDKitToolkitWrapper.apply_elf_conformer_selection
        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            toolkit_registry.call(
                "apply_elf_conformer_selection",
                molecule=self,
                percentage=percentage,
                limit=limit,
                **kwargs,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            toolkit.apply_elf_conformer_selection(  # type: ignore[attr-defined]
                molecule=self, percentage=percentage, limit=limit, **kwargs
            )
        else:
            raise InvalidToolkitRegistryError(
                f"Invalid toolkit_registry passed to apply_elf_conformer_selection."
                f"Expected ToolkitRegistry or ToolkitWrapper. Got "
                f"{type(toolkit_registry)}"
            )

    def compute_partial_charges_am1bcc(
        self,
        use_conformers=None,
        strict_n_conformers=False,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        .. deprecated:: 0.11.0

            This method was deprecated in v0.11.0 and will soon be removed.
            Use :py:meth:`assign_partial_charges(partial_charge_method='am1bcc')
            <Molecule.assign_partial_charges>` instead.

        Calculate partial atomic charges for this molecule using AM1-BCC run by an underlying toolkit
        and assign them to this molecule's ``partial_charges`` attribute.

        Parameters
        ----------
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is provided for the given charge method.
            If this is False and an invalid number of conformers is found, a warning will be raised.
        use_conformers : iterable of openmm.unit.Quantity-wrapped numpy arrays, each with shape (n_atoms, 3)
            and dimension of distance. Optional, default=None Coordinates to use for partial charge calculation.
            If None, an appropriate number of conformers for the given charge method will be generated.
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry
        or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for the calculation

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.generate_conformers()
        >>> molecule.compute_partial_charges_am1bcc()

        Raises
        ------
        InvalidToolkitRegistryError
            If an invalid object is passed as the toolkit_registry parameter

        """
        # TODO: Remove in version 0.12.0
        warnings.warn(
            "compute_partial_charges_am1bcc is deprecated and will be removed in version 0.12.0. "
            "Use assign_partial_charges(partial_charge_method='am1bcc') instead.",
            UserWarning,
        )
        self.assign_partial_charges(
            partial_charge_method="am1bcc",
            use_conformers=use_conformers,
            strict_n_conformers=strict_n_conformers,
            toolkit_registry=toolkit_registry,
        )

    def assign_partial_charges(
        self,
        partial_charge_method: str,
        strict_n_conformers=False,
        use_conformers=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        normalize_partial_charges=True,
    ):
        """
        Calculate partial atomic charges and store them in the molecule.

        ``assign_partial_charges`` computes charges using the specified toolkit
        and assigns the new values to the ``partial_charges`` attribute.
        Supported charge methods vary from toolkit to toolkit, but some
        supported methods are:

        - ``"am1bcc"``
        - ``"am1bccelf10"`` (requires OpenEye Toolkits)
        - ``"am1-mulliken"``
        - ``"mmff94"``
        - ``"gasteiger"``

        For more supported charge methods and details, see the corresponding
        methods in each toolkit wrapper:

        - :meth:`OpenEyeToolkitWrapper.assign_partial_charges \
          <openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.assign_partial_charges>`
        - :meth:`RDKitToolkitWrapper.assign_partial_charges \
          <openff.toolkit.utils.toolkits.RDKitToolkitWrapper.assign_partial_charges>`
        - :meth:`AmberToolsToolkitWrapper.assign_partial_charges \
          <openff.toolkit.utils.toolkits.AmberToolsToolkitWrapper.assign_partial_charges>`
        - :meth:`BuiltInToolkitWrapper.assign_partial_charges \
          <openff.toolkit.utils.toolkits.BuiltInToolkitWrapper.assign_partial_charges>`

        Parameters
        ----------
        partial_charge_method : string
            The partial charge calculation method to use for partial charge
            calculation.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is
            provided for the given charge method. If this is False and an
            invalid number of conformers is found, a warning will be raised.
        use_conformers : iterable of openmm.unit.Quantity-wrapped numpy arrays, each with shape (n_atoms, 3) and
            dimension of distance. Optional, default=None
            Coordinates to use for partial charge calculation. If None, an
            appropriate number of conformers will be generated.
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or
            openff.toolkit.utils.toolkits.ToolkitWrapper,
            optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for the
            calculation.
        normalize_partial_charges : bool, default=True
            Whether to offset partial charges so that they sum to the total
            formal charge of the molecule. This is used to prevent accumulation
            of rounding errors when the partial charge assignment method returns
            values at limited precision.

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.assign_partial_charges('am1-mulliken')

        Raises
        ------
        InvalidToolkitRegistryError
            If an invalid object is passed as the toolkit_registry parameter

        See Also
        --------
        openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.assign_partial_charges
        openff.toolkit.utils.toolkits.RDKitToolkitWrapper.assign_partial_charges
        openff.toolkit.utils.toolkits.AmberToolsToolkitWrapper.assign_partial_charges
        openff.toolkit.utils.toolkits.BuiltInToolkitWrapper.assign_partial_charges
        """

        # Raise a warning when users try to apply these charge methods to "large" molecules
        WARN_LARGE_MOLECULES: Set[str] = {
            "am1bcc",
            "am1bccelf10",
            "am1-mulliken",
            "am1bccnosymspt",
            "am1elf10",
        }

        if partial_charge_method in WARN_LARGE_MOLECULES:
            if self.n_atoms > 150:
                warnings.warn(
                    f"Warning! Partial charge method '{partial_charge_method}' is not designed "
                    "for use on large (i.e. > 150 atoms) molecules and may crash or take hours to "
                    f"run on this molecule (found {self.n_atoms} atoms). For more, see "
                    "https://docs.openforcefield.org/projects/toolkit/en/stable/faq.html"
                    "#parameterizing-my-system-which-contains-a-large-molecule-is-taking-forever-whats-wrong",
                )

        if isinstance(toolkit_registry, ToolkitRegistry):
            # We may need to try several toolkitwrappers to find one
            # that supports the desired partial charge method, so we
            # tell the ToolkitRegistry to continue trying ToolkitWrappers
            # if one raises an error (raise_exception_types=[])
            toolkit_registry.call(
                "assign_partial_charges",
                molecule=self,
                partial_charge_method=partial_charge_method,
                use_conformers=use_conformers,
                strict_n_conformers=strict_n_conformers,
                normalize_partial_charges=normalize_partial_charges,
                raise_exception_types=[],
                _cls=self.__class__,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit_wrapper: ToolkitWrapper = toolkit_registry
            toolkit_wrapper.assign_partial_charges(  # type: ignore[attr-defined]
                self,
                partial_charge_method=partial_charge_method,
                use_conformers=use_conformers,
                strict_n_conformers=strict_n_conformers,
                normalize_partial_charges=normalize_partial_charges,
                _cls=self.__class__,
            )
        else:
            raise InvalidToolkitRegistryError(
                f"Invalid toolkit_registry passed to assign_partial_charges."
                f"Expected ToolkitRegistry or ToolkitWrapper. Got  {type(toolkit_registry)}"
            )

    def _normalize_partial_charges(self):
        """
        Add offsets to each partial charge to ensure that they sum to the formal charge of the molecule,
        to the limit of a python float's precision. Modifies the partial charges in-place.
        """
        expected_charge = self.total_charge

        current_charge = 0.0 * unit.elementary_charge
        for pc in self.partial_charges:
            current_charge += pc

        charge_offset = (expected_charge - current_charge) / self.n_atoms

        self.partial_charges += charge_offset

    def assign_fractional_bond_orders(
        self,
        bond_order_model=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        use_conformers=None,
    ):
        """
        Update and store list of bond orders this molecule.

        Bond orders are stored on each bond, in the
        ``bond.fractional_bond_order`` attribute.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or
            openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion
        bond_order_model : string, optional. Default=None
            The bond order model to use for fractional bond order calculation. If ``None``, ``"am1-wiberg"`` is used.
        use_conformers : iterable of openmm.unit.Quantity(np.array) with shape (n_atoms, 3) and dimension of distance,
            optional, default=None
            The conformers to use for fractional bond order calculation. If ``None``, an appropriate number
            of conformers will be generated by an available ``ToolkitWrapper``.

        Examples
        --------

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.assign_fractional_bond_orders()

        Raises
        ------
        InvalidToolkitRegistryError
            If an invalid object is passed as the toolkit_registry parameter

        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call(
                "assign_fractional_bond_orders",
                self,
                bond_order_model=bond_order_model,
                use_conformers=use_conformers,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            return toolkit.assign_fractional_bond_orders(
                self, bond_order_model=bond_order_model, use_conformers=use_conformers
            )
        else:
            raise InvalidToolkitRegistryError(
                f"Invalid toolkit_registry passed to assign_fractional_bond_orders. "
                f"Expected ToolkitRegistry or ToolkitWrapper. Got {type(toolkit_registry)}."
            )

    def _invalidate_cached_properties(self):
        """
        Indicate that the chemical entity has been altered.
        """
        # if hasattr(self, '_cached_properties'):
        #    delattr(self, '_cached_properties')
        self._conformers = None
        self._partial_charges = None
        self._propers = None
        self._impropers = None

        self._cached_smiles = None
        # TODO: Clear fractional bond orders
        self._ordered_connection_table_hash = None
        for atom in self.atoms:
            if "molecule_atom_index" in atom.__dict__:
                del atom.__dict__["molecule_atom_index"]

    def to_networkx(self):
        """Generate a NetworkX undirected graph from the molecule.

        Nodes are Atoms labeled with atom indices and atomic elements (via the ``element`` node atrribute).
        Edges denote chemical bonds between Atoms.

        .. todo ::

           * Do we need a ``from_networkx()`` method? If so, what would the Graph be required to provide?
           * Should edges be labeled with discrete bond types in some aromaticity model?
           * Should edges be labeled with fractional bond order if a method is specified?
           * Should we add other per-atom and per-bond properties (e.g. partial charges) if present?
           * Can this encode bond/atom chirality?


        Returns
        -------
        graph : networkx.Graph
            The resulting graph, with nodes (atoms) labeled with atom indices, elements, stereochemistry and
            aromaticity flags and bonds with two atom indices, bond order, stereochemistry, and aromaticity flags

        Examples
        --------
        Retrieve the bond graph for imatinib (OpenEye toolkit required)

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> nxgraph = molecule.to_networkx()

        """
        import networkx as nx

        G = nx.Graph()
        for atom in self.atoms:
            G.add_node(
                atom.molecule_atom_index,
                atomic_number=atom.atomic_number,
                is_aromatic=atom.is_aromatic,
                stereochemistry=atom.stereochemistry,
                formal_charge=atom.formal_charge,
            )
            # G.add_node(atom.molecule_atom_index, attr_dict={'atomic_number': atom.atomic_number})
        for bond in self.bonds:
            G.add_edge(
                bond.atom1_index,
                bond.atom2_index,
                bond_order=bond.bond_order,
                is_aromatic=bond.is_aromatic,
                stereochemistry=bond.stereochemistry,
            )
            # G.add_edge(bond.atom1_index, bond.atom2_index, attr_dict={'order':bond.bond_order})

        return G

    def find_rotatable_bonds(
        self, ignore_functional_groups=None, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY
    ):
        """
        Find all bonds classed as rotatable ignoring any matched to the ``ignore_functional_groups`` list.

        Parameters
        ----------
        ignore_functional_groups: optional, List[str], default=None,
            A list of bond SMARTS patterns to be ignored when finding rotatable bonds.

        toolkit_registry: openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapperl, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMARTS matching

        Returns
        -------
        bonds: List[openff.toolkit.topology.molecule.Bond]
            The list of openff.toolkit.topology.molecule.Bond instances which are rotatable.
        """

        # general rotatable bond smarts taken from RDKit
        # https://github.com/rdkit/rdkit/blob/1bf6ef3d65f5c7b06b56862b3fb9116a3839b229/rdkit/Chem/Lipinski.py#L47%3E
        rotatable_bond_smarts = "[!$(*#*)&!D1:1]-&!@[!$(*#*)&!D1:2]"

        # get all of the general matches
        general_matches = self.chemical_environment_matches(
            query=rotatable_bond_smarts, toolkit_registry=toolkit_registry
        )

        # this will give all forwards and backwards matches, so condense them down with this function
        def condense_matches(matches):
            condensed_matches = set()
            for m in matches:
                condensed_matches.add(tuple(sorted(m)))
            return condensed_matches

        general_bonds = condense_matches(general_matches)

        # now refine the list using the ignore groups
        if ignore_functional_groups is not None:
            matches_to_ignore = set()

            # make ignore_functional_groups an iterable object
            if isinstance(ignore_functional_groups, str):
                ignore_functional_groups = [ignore_functional_groups]
            else:
                try:
                    iter(ignore_functional_groups)
                except TypeError:
                    ignore_functional_groups = [ignore_functional_groups]

            # find the functional groups to remove
            for functional_group in ignore_functional_groups:
                # note I run the searches through this function so they have to be SMIRKS?
                ignore_matches = self.chemical_environment_matches(
                    query=functional_group, toolkit_registry=toolkit_registry
                )
                ignore_matches = condense_matches(ignore_matches)
                # add the new matches to the matches to ignore
                matches_to_ignore.update(ignore_matches)

            # now remove all the matches
            for match in matches_to_ignore:
                try:
                    general_bonds.remove(match)
                # if the key is not in the list, the ignore pattern was not valid
                except KeyError:
                    continue

        # gather a list of bond instances to return
        rotatable_bonds = [self.get_bond_between(*bond) for bond in general_bonds]
        return rotatable_bonds

    def _add_atom(
        self,
        atomic_number,
        formal_charge,
        is_aromatic,
        stereochemistry=None,
        name=None,
        metadata=None,
        invalidate_cache: bool = True,
    ):
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
            Either 'R' or 'S' for specified stereochemistry, or None if stereochemistry is irrelevant
        name : str, optional, default=None
            An optional name for the atom
        metadata : dict[str: (int, str)], default=None
            An optional dictionary where keys are strings and values are strings or ints. This is intended
            to record atom-level information used to inform hierarchy definition and iteration, such as
            grouping atom by residue and chain.
        invalidate_cache : bool, default=True
            Whether or not to invalidate the cache of the molecule upon the addition of this atom. This should
            be left to its default value (`True`) for safety.

        Returns
        -------
        index : int
            The index of the atom in the molecule

        Examples
        --------

        Define a methane molecule

        >>> molecule = Molecule()
        >>> molecule.name = 'methane'
        >>> C = molecule.add_atom(6, 0, False)
        >>> H1 = molecule.add_atom(1, 0, False)
        >>> H2 = molecule.add_atom(1, 0, False)
        >>> H3 = molecule.add_atom(1, 0, False)
        >>> H4 = molecule.add_atom(1, 0, False)
        >>> bond_idx = molecule.add_bond(C, H1, False, 1)
        >>> bond_idx = molecule.add_bond(C, H2, False, 1)
        >>> bond_idx = molecule.add_bond(C, H3, False, 1)
        >>> bond_idx = molecule.add_bond(C, H4, False, 1)

        """
        # Create an atom
        atom = Atom(
            atomic_number,
            formal_charge,
            is_aromatic,
            stereochemistry=stereochemistry,
            name=name,
            metadata=metadata,
            molecule=self,
        )
        self._atoms.append(atom)
        if invalidate_cache:
            self._invalidate_cached_properties()

        # Since we just appended it, we can just return the length - 1
        return len(self._atoms) - 1

    def _add_bond(
        self,
        atom1,
        atom2,
        bond_order,
        is_aromatic,
        stereochemistry=None,
        fractional_bond_order=None,
        invalidate_cache: bool = True,
    ):
        """
        Add a bond between two specified atom indices

        Parameters
        ----------
        atom1 : int or openff.toolkit.topology.molecule.Atom
            Index of first atom or first atom
        atom2_index : int or openff.toolkit.topology.molecule.Atom
            Index of second atom or second atom
        bond_order : int
            Integral bond order of Kekulized form
        is_aromatic : bool
            True if this bond is aromatic, False otherwise
        stereochemistry : str, optional, default=None
            Either 'E' or 'Z' for specified stereochemistry, or None if stereochemistry is irrelevant
        fractional_bond_order : float, optional, default=None
            The fractional (eg. Wiberg) bond order
        invalidate_cache : bool, default=True
            Whether or not to invalidate the cache of the molecule upon the addition of this atom. This should
            be left to its default value (`True`) for safety.

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
            raise Exception(
                "Invalid inputs to molecule._add_bond. Expected ints or Atoms. "
                f"Received {atom1} (type {type(atom1)}) and {atom2} (type {type(atom2)}) "
            )
        # TODO: Check to make sure bond does not already exist
        if atom1_atom.is_bonded_to(atom2_atom):
            raise Exception(
                f"Bond already exists between {atom1_atom} and {atom2_atom})"
            )
        bond = Bond(
            atom1_atom,
            atom2_atom,
            bond_order,
            is_aromatic,
            stereochemistry=stereochemistry,
            fractional_bond_order=fractional_bond_order,
        )
        self._bonds.append(bond)
        if invalidate_cache:
            self._invalidate_cached_properties()

        # Since we just appended it, we can just return the length - 1
        return len(self._bonds) - 1

    def _add_conformer(self, coordinates):
        """
        Add a conformation of the molecule

        Parameters
        ----------
        coordinates: openmm.unit.Quantity(np.array) with shape (n_atoms, 3) and dimension of distance
            Coordinates of the new conformer, with the first dimension of the array corresponding to the atom index in
            the molecule's indexing system.

        Returns
        -------
        index: int
            The index of this conformer
        """
        if coordinates.shape != (self.n_atoms, 3):
            raise InvalidConformerError(
                "molecule.add_conformer given input of the wrong shape: "
                f"Given {coordinates.shape}, expected {(self.n_atoms, 3)}"
            )

        if isinstance(coordinates, unit.Quantity):
            if not coordinates.units.is_compatible_with(unit.angstrom):
                raise IncompatibleUnitError(
                    "Coordinates passed to Molecule._add_conformer with incompatible units. "
                    "Ensure that units are dimension of length."
                )

        elif hasattr(coordinates, "unit"):
            from openff.units.openmm import from_openmm
            from openmm import unit as openmm_unit

            if not isinstance(coordinates, openmm_unit.Quantity):
                raise IncompatibleUnitError(
                    "Unsupported type passed to Molecule._add_conformer setter. "
                    "Found object of type {type(other)}."
                )

            if not coordinates.unit.is_compatible(openmm_unit.meter):
                raise IncompatibleUnitError(
                    "Coordinates passed to Molecule._add_conformer with units of incompatible dimensionality. "
                    f"Adding conformers with OpenMM-style units is supported, by found units of {coordinates.unit}. "
                    "Ensure that units are dimension of length."
                )

            coordinates = from_openmm(coordinates)

        else:
            raise IncompatibleUnitError(
                "Unknown object passed to Molecule._add_conformer. Expected types include "
                f"openmm.unit.Quantity and openff.units.unit.Quantity, found type {type(coordinates)}."
            )

        tmp_conf = unit.Quantity(
            np.zeros(shape=(self.n_atoms, 3), dtype=float), unit.angstrom
        )
        try:
            tmp_conf[:] = coordinates
        except AttributeError as e:
            print(e)

        if self._conformers is None:
            # TODO should we checking that the exact same conformer is not in the list already?
            self._conformers = []
        self._conformers.append(tmp_conf)
        return len(self._conformers)

    @property
    def partial_charges(self):
        """
        Returns the partial charges (if present) on the molecule.

        Returns
        -------
        partial_charges : a openmm.unit.Quantity - wrapped numpy array [1 x n_atoms] or None
            The partial charges on the molecule's atoms. Returns None if no charges have been specified.
        """
        return self._partial_charges

    @partial_charges.setter
    def partial_charges(self, charges):
        """
        Set the atomic partial charges for this molecule.

        Parameters
        ----------
        charges : None or a openmm.unit.Quantity - wrapped numpy array [1 x n_atoms]
            The partial charges to assign to the molecule. If not None, must be in units compatible with
            openmm.unit.elementary_charge

        """
        if charges is None:
            self._partial_charges = None
        elif charges.shape == (self.n_atoms,):
            if isinstance(charges, unit.Quantity):
                if charges.units in unit.elementary_charge.compatible_units():
                    self._partial_charges = charges
            if hasattr(charges, "unit"):
                from openmm import unit as openmm_unit

                if not isinstance(charges, openmm_unit.Quantity):
                    raise IncompatibleUnitError(
                        "Unsupported type passed to partial_charges setter. "
                        "Found object of type {type(charges)}."
                    )

                elif isinstance(charges, openmm_unit.Quantity):
                    from openff.units.openmm import from_openmm

                    converted = from_openmm(charges)
                    if converted.units in unit.elementary_charge.compatible_units():
                        self._partial_charges = converted

    @property
    def n_particles(self) -> int:
        """DEPRECATED: Use Molecule.n_atoms instead."""
        _molecule_deprecation("n_particles", "n_atoms")
        return self.n_atoms

    @property
    def n_atoms(self) -> int:
        """
        The number of Atom objects.
        """
        return len(self._atoms)

    @property
    def n_bonds(self):
        """
        The number of Bond objects in the molecule.
        """
        return sum([1 for bond in self.bonds])

    @property
    def n_angles(self) -> int:
        """Number of angles in the molecule."""
        self._construct_angles()
        return len(self._angles)

    @property
    def n_propers(self) -> int:
        """Number of proper torsions in the molecule."""
        self._construct_torsions()
        assert (
            self._propers is not None
        ), "_construct_torsions always sets _propers to a set"
        return len(self._propers)

    @property
    def n_impropers(self) -> int:
        """Number of possible improper torsions in the molecule."""
        self._construct_torsions()
        assert (
            self._impropers is not None
        ), "_construct_torsions always sets _impropers to a set"
        return len(self._impropers)

    @property
    def particles(self) -> List[Atom]:
        """DEPRECATED: Use Molecule.atoms instead."""
        _molecule_deprecation("particles", "atoms")
        return self.atoms

    def particle(self, index: int) -> Atom:
        """DEPRECATED: Use Molecule.atom instead."""
        _molecule_deprecation("particle", "atom")
        return self.atom(index)

    def particle_index(self, particle: Atom) -> int:
        """DEPRECATED: Use Molecule.atom_index instead."""
        _molecule_deprecation("particle_index", "atom_index")
        return self.atom_index(particle)

    @property
    def atoms(self):
        """
        Iterate over all Atom objects in the molecule.
        """
        return self._atoms

    def atom(self, index: int) -> Atom:
        """
        Get the atom with the specified index.

        Parameters
        ----------
        index : int

        Returns
        -------
        atom : openff.toolkit.topology.Atom
        """
        return self._atoms[index]

    def atom_index(self, atom: Atom) -> int:
        """
        Returns the index of the given atom in this molecule

        .. TODO: document behaviour when atom is not present in self

        Parameters
        ----------
        atom : openff.toolkit.topology.Atom

        Returns
        -------
        index : int
            The index of the given atom in this molecule
        """
        return atom.molecule_atom_index

    @property
    def conformers(self):
        """
        Returns the list of conformers for this molecule.

        Conformers are presented as a list of ``Quantity``-wrapped NumPy
        arrays, of shape (3 x n_atoms) and with dimensions of [Distance]. The
        return value is the actual list of conformers, and changes to the
        contents affect the original ``FrozenMolecule``.
        """
        return self._conformers

    @property
    def n_conformers(self) -> int:
        """
        The number of conformers for this molecule.
        """
        if self._conformers is None:
            return 0
        return len(self._conformers)

    @property
    def bonds(self) -> List[Bond]:
        """
        Iterate over all Bond objects in the molecule.
        """
        return self._bonds

    def bond(self, index: int) -> Bond:
        """
        Get the bond with the specified index.

        Parameters
        ----------
        index : int

        Returns
        -------
        bond : openff.toolkit.topology.Bond
        """
        return self._bonds[index]

    @property
    def angles(self) -> Set[Tuple[Atom, Atom, Atom]]:
        """
        Get an iterator over all i-j-k angles.
        """
        self._construct_angles()
        return self._angles

    @property
    def torsions(self) -> Set[Tuple[Atom, Atom, Atom, Atom]]:
        """
        Get an iterator over all i-j-k-l torsions.
        Note that i-j-k-i torsions (cycles) are excluded.

        Returns
        -------
        torsions : iterable of 4-Atom tuples
        """
        self._construct_torsions()
        assert (
            self._torsions is not None
        ), "_construct_torsions always sets _torsions to a set"
        return self._torsions

    @property
    def propers(self) -> Set[Tuple[Atom, Atom, Atom, Atom]]:
        """
        Iterate over all proper torsions in the molecule

        .. todo::

           * Do we need to return a ``Torsion`` object that collects information about fractional bond orders?
        """
        self._construct_torsions()
        assert (
            self._propers is not None
        ), "_construct_torsions always sets _propers to a set"
        return self._propers

    @property
    def impropers(self) -> Set[Tuple[Atom, Atom, Atom, Atom]]:
        """
        Iterate over all improper torsions in the molecule.

        .. todo ::
           * Do we need to return a ``Torsion`` object that collects information about fractional bond orders?

        Returns
        -------
        impropers : set of tuple
            An iterator of tuples, each containing the atoms making
            up a possible improper torsion.

        See Also
        --------
        smirnoff_impropers, amber_impropers
        """
        self._construct_torsions()
        assert (
            self._impropers is not None
        ), "_construct_torsions always sets _impropers to a set"
        return self._impropers

    @property
    def smirnoff_impropers(self) -> Set[Tuple[Atom, Atom, Atom, Atom]]:
        """
        Iterate over all impropers with trivalent centers, reporting the central atom second.

        The central atom is reported second in each torsion. This method reports
        an improper for each trivalent atom in the molecule, whether or not any
        given force field would assign it improper torsion parameters.

        Also note that this will return 6 possible atom orderings around each improper
        center. In current SMIRNOFF parameterization, three of these six
        orderings will be used for the actual assignment of the improper term
        and measurement of the angles. These three orderings capture the three unique
        angles that could be calculated around the improper center, therefore the sum
        of these three terms will always return a consistent energy.

        The exact three orderings that will be applied during parameterization can not be
        determined in this method, since it requires sorting the atom indices, and
        those indices may change when this molecule is added to a Topology.

        For more details on the use of three-fold ('trefoil') impropers, see
        https://openforcefield.github.io/standards/standards/smirnoff/#impropertorsions

        Returns
        -------
        impropers : set of tuple
            An iterator of tuples, each containing the indices of atoms making
            up a possible improper torsion. The central atom is listed second
            in each tuple.

        See Also
        --------
        impropers, amber_impropers

        """
        # TODO: Replace with non-cheminformatics-toolkit method
        #       (ie. just looping over all atoms and finding ones that have 3 bonds?)

        smirnoff_improper_smarts = "[*:1]~[X3:2](~[*:3])~[*:4]"
        improper_idxs = self.chemical_environment_matches(smirnoff_improper_smarts)
        smirnoff_impropers = {
            (self.atom(imp[0]), self.atom(imp[1]), self.atom(imp[2]), self.atom(imp[3]))
            for imp in improper_idxs
        }
        return smirnoff_impropers

    @property
    def amber_impropers(self) -> Set[Tuple[Atom, Atom, Atom, Atom]]:
        """
        Iterate over all impropers with trivalent centers, reporting the central atom first.

        The central atom is reported first in each torsion. This method reports
        an improper for each trivalent atom in the molecule, whether or not any
        given force field would assign it improper torsion parameters.

        Also note that this will return 6 possible atom orderings around each
        improper center. In current AMBER parameterization, one of these six
        orderings will be used for the actual assignment of the improper term
        and measurement of the angle. This method does not encode the logic to
        determine which of the six orderings AMBER would use.

        Returns
        -------
        impropers : set of tuple
            An iterator of tuples, each containing the indices of atoms making
            up a possible improper torsion. The central atom is listed first in
            each tuple.

        See Also
        --------
        impropers, smirnoff_impropers

        """
        # TODO: Replace with non-cheminformatics-toolkit method
        #       (ie. just looping over all atoms and finding ones that have 3 bonds?)
        amber_improper_smarts = "[X3:1](~[*:2])(~[*:3])~[*:4]"
        improper_idxs = self.chemical_environment_matches(amber_improper_smarts)
        amber_impropers = {
            (self.atom(imp[0]), self.atom(imp[1]), self.atom(imp[2]), self.atom(imp[3]))
            for imp in improper_idxs
        }
        return amber_impropers

    def nth_degree_neighbors(self, n_degrees):
        """
        Return canonicalized pairs of atoms whose shortest separation is `exactly` n bonds.
        Only pairs with increasing atom indices are returned.

        Parameters
        ----------
        n: int
            The number of bonds separating atoms in each pair

        Returns
        -------
        neighbors: iterator of tuple of Atom
            Tuples (len 2) of atom that are separated by ``n`` bonds.

        Notes
        -----

        The criteria used here relies on minimum distances; when there are multiple valid
        paths between atoms, such as atoms in rings, the shortest path is considered.
        For example, two atoms in "meta" positions with respect to each other in a benzene
        are separated by two paths, one length 2 bonds and the other length 4 bonds. This
        function would consider them to be 2 apart and would not include them if ``n=4`` was
        passed.

        """
        if n_degrees <= 0:
            raise ValueError(
                "Cannot consider neighbors separated by 0 or fewer atoms. Asked to consider "
                f"path lengths of {n_degrees}."
            )
        else:
            return _nth_degree_neighbors_from_graphlike(
                graphlike=self, n_degrees=n_degrees
            )

    @property
    def total_charge(self):
        """
        Return the total charge on the molecule
        """
        charge_sum = 0.0 * unit.elementary_charge
        for atom in self.atoms:
            charge_sum += atom.formal_charge
        return charge_sum

    @property
    def name(self) -> str:
        """
        The name (or title) of the molecule
        """
        return self._name

    @name.setter
    def name(self, other):
        """
        Set the name of this molecule
        """
        if other is None:
            self._name = ""
        elif type(other) is str:
            self._name = other
        else:
            raise Exception("Molecule name must be a string")

    @property
    def properties(self) -> Dict[str, Any]:
        """
        The properties dictionary of the molecule
        """
        return self._properties

    @property
    def hill_formula(self) -> str:
        """
        Get the Hill formula of the molecule
        """
        return self.to_hill_formula()

    def to_hill_formula(self) -> str:
        """
        Generate the Hill formula of this molecule.

        Returns
        ----------
        formula : the Hill formula of the molecule

        Raises
        -----------
        NotImplementedError : if the molecule is not of one of the specified types.
        """
        atom_nums = [atom.atomic_number for atom in self.atoms]

        return _atom_nums_to_hill_formula(atom_nums)

    @staticmethod
    def _object_to_hill_formula(obj: Union["Molecule", "nx.Graph"]) -> str:
        """Take a Molecule or NetworkX graph and generate its Hill formula.
        This provides a backdoor to the old functionality of Molecule.to_hill_formula, which
        was a static method that duck-typed inputs of Molecule or graph objects."""
        import networkx as nx

        if isinstance(obj, FrozenMolecule):
            return obj.to_hill_formula()
        elif isinstance(obj, nx.Graph):
            return _networkx_graph_to_hill_formula(obj)
        else:
            raise RuntimeError(
                f"Unsupport object of type {type(obj)} passed to "
                "Molecule._object_to_hill_formula"
            )

    def chemical_environment_matches(
        self,
        query,
        unique=False,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """Find matches in the molecule for a SMARTS string or ``ChemicalEnvironment`` query

        Parameters
        ----------
        query : str or ChemicalEnvironment
            SMARTS string (with one or more tagged atoms) or ``ChemicalEnvironment`` query.
            Query will internally be resolved to SMIRKS using ``query.asSMIRKS()`` if it has an ``.asSMIRKS`` method.
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for chemical environment matches

        Returns
        -------
        matches : list of atom index tuples
            A list of tuples, containing the indices of the matching atoms.

        Examples
        --------
        Retrieve all the carbon-carbon bond matches in a molecule

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> matches = molecule.chemical_environment_matches('[#6X3:1]~[#6X3:2]')

        .. todo ::

           * Do we want to generalize ``query`` to allow other kinds of queries,
             such as mdtraj DSL, pymol selections, atom index slices, etc? We
             could call it ``topology.matches(query)`` instead of
             ``chemical_environment_matches``

        """
        # Resolve to SMIRKS if needed
        # TODO: Update this to use updated ChemicalEnvironment API
        if hasattr(query, "smirks"):
            smirks = query.smirks
        elif type(query) == str:
            smirks = query
        else:
            raise ValueError("'query' must be either a string or a ChemicalEnvironment")

        # Use specified cheminformatics toolkit to determine matches with specified aromaticity model
        # TODO: Simplify this by requiring a toolkit registry for the molecule?
        # TODO: Do we have to pass along an aromaticity model?
        if isinstance(toolkit_registry, ToolkitRegistry):
            matches = toolkit_registry.call(
                "find_smarts_matches",
                self,
                smirks,
                unique=unique,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            matches = toolkit_registry.find_smarts_matches(
                self,
                smirks,
                unique=unique,
            )
        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        return matches

    @classmethod
    def from_iupac(
        cls,
        iupac_name,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
        **kwargs,
    ):
        """Generate a molecule from IUPAC or common name

        .. note :: This method requires the OpenEye toolkit to be installed.

        Parameters
        ----------
        iupac_name : str
            IUPAC name of molecule to be generated
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for chemical environment matches
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if molecule contains undefined stereochemistry.

        Returns
        -------
        molecule : Molecule
            The resulting molecule with position

        Examples
        --------

        Create a molecule from an IUPAC name

        >>> molecule = Molecule.from_iupac('4-[(4-methylpiperazin-1-yl)methyl]-N-(4-methyl-3-{[4-(pyridin-3-yl)pyrimidin-2-yl]amino}phenyl)benzamide')  # noqa

        Create a molecule from a common name

        >>> molecule = Molecule.from_iupac('imatinib')

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            molecule = toolkit_registry.call(
                "from_iupac",
                iupac_name,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
                **kwargs,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            molecule = toolkit.from_iupac(
                iupac_name,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
                **kwargs,
            )
        else:
            raise Exception(
                "Invalid toolkit_registry passed to from_iupac. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}."
            )

        return molecule

    def to_iupac(self, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """Generate IUPAC name from Molecule

        Returns
        -------
        iupac_name : str
            IUPAC name of the molecule

        .. note :: This method requires the OpenEye toolkit to be installed.

        Examples
        --------

        >>> from openff.toolkit.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> iupac_name = molecule.to_iupac()

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            to_iupac_method = toolkit_registry.resolve("to_iupac")
        elif isinstance(toolkit_registry, ToolkitWrapper):
            to_iupac_method = toolkit_registry.to_iupac
        else:
            raise Exception(
                "Invalid toolkit_registry passed to to_iupac. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        # TODO: Can `to_iupac` fail if given a well-behaved OFFMol/OEMol?
        result = to_iupac_method(self)
        return result

    @classmethod
    def from_topology(cls, topology):
        """Return a Molecule representation of an OpenFF Topology containing a single Molecule object.

        Parameters
        ----------
        topology : openff.toolkit.topology.Topology
            The :class:`Topology` object containing a single :class:`Molecule` object.
            Note that OpenMM and MDTraj ``Topology`` objects are not supported.

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            The Molecule object in the topology

        Raises
        ------
        ValueError
            If the topology does not contain exactly one molecule.

        Examples
        --------

        Create a molecule from a Topology object that contains exactly one molecule

        >>> molecule = Molecule.from_topology(topology)  # doctest: +SKIP

        """
        # TODO: Ensure we are dealing with an OpenFF Topology object
        if topology.n_molecules != 1:
            raise ValueError("Topology must contain exactly one molecule")
        molecule = [i for i in topology.molecules][0]
        return cls(molecule)

    def to_topology(self):
        """
        Return an OpenFF Topology representation containing one copy of this molecule

        Returns
        -------
        topology : openff.toolkit.topology.Topology
            A Topology representation of this molecule

        Examples
        --------

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> topology = molecule.to_topology()

        """
        from openff.toolkit.topology import Topology

        return Topology.from_molecules(self)

    @classmethod
    def from_file(
        cls,
        file_path,
        file_format=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
    ):
        """
        Create one or more molecules from a file

        .. todo::

           * Extend this to also include some form of .offmol Open Force Field Molecule format?
           * Generalize this to also include file-like objects?

        Parameters
        ----------
        file_path : str or file-like object
            The path to the file or file-like object to stream one or more molecules from.
        file_format : str, optional, default=None
            Format specifier, usually file suffix (eg. 'MOL2', 'SMI')
            Note that not all toolkits support all formats. Check ToolkitWrapper.toolkit_file_read_formats for your
            loaded toolkits for details.
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or
            openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for file loading. If a Toolkit is passed, only
            the highest-precedence toolkit is used
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecules : Molecule or list of Molecules
            If there is a single molecule in the file, a Molecule is returned;
            otherwise, a list of Molecule objects is returned.

        Examples
        --------
        >>> from openff.toolkit.tests.utils import get_monomer_mol2_file_path
        >>> mol2_file_path = get_monomer_mol2_file_path('cyclohexane')
        >>> molecule = Molecule.from_file(mol2_file_path)

        """

        if file_format is None:
            if not isinstance(file_path, str):
                raise Exception(
                    "If providing a file-like object for reading molecules, the format must be specified"
                )
            # Assume that files ending in ".gz" should use their second-to-last suffix for compatibility check
            # TODO: Will all cheminformatics packages be OK with gzipped files?
            if file_path[-3:] == ".gz":
                file_format = file_path.split(".")[-2]
            else:
                file_format = file_path.split(".")[-1]
        file_format = file_format.upper()

        if file_format == "XYZ":
            raise UnsupportedFileTypeError(
                "Parsing `.xyz` files is not currently supported because they lack sufficient "
                "chemical information to be used with SMIRNOFF force fields. For more information, "
                "see https://open-forcefield-toolkit.readthedocs.io/en/latest/faq.html or to provide "
                "feedback please visit https://github.com/openforcefield/openff-toolkit/issues/1145."
            )

        # Determine which toolkit to use (highest priority that's compatible with input type)
        if isinstance(toolkit_registry, ToolkitRegistry):
            # TODO: Encapsulate this logic into ToolkitRegistry.call()?
            toolkit = None
            supported_read_formats = {}
            for query_toolkit in toolkit_registry.registered_toolkits:
                if file_format in query_toolkit.toolkit_file_read_formats:
                    toolkit = query_toolkit
                    break
                supported_read_formats[
                    query_toolkit.toolkit_name
                ] = query_toolkit.toolkit_file_read_formats
            if toolkit is None:
                msg = (
                    f"No toolkits in registry can read file {file_path} (format {file_format}). Supported "
                    f"formats in the provided ToolkitRegistry are {supported_read_formats}. "
                )
                # Per issue #407, not allowing RDKit to read mol2 has confused a lot of people. Here we add text
                # to the error message that will hopefully reduce this confusion.
                if file_format == "MOL2" and RDKitToolkitWrapper.is_available():
                    msg += (
                        "RDKit does not fully support input of molecules from mol2 format unless they "
                        "have Corina atom types, and this is not common in the simulation community. For this "
                        "reason, the Open Force Field Toolkit does not use "
                        "RDKit to read .mol2. Consider reading from SDF instead. If you would like to attempt "
                        "to use RDKit to read mol2 anyway, you can load the molecule of interest into an RDKit "
                        "molecule and use openff.toolkit.topology.Molecule.from_rdkit, but we do not recommend this."
                    )
                elif file_format == "PDB" and RDKitToolkitWrapper.is_available():
                    msg += (
                        "RDKit can not safely read PDBs on their own. Information about bond order and aromaticity "
                        "is likely to be lost. PDBs can be used along with a valid smiles string with RDKit using "
                        "the constructor Molecule.from_pdb_and_smiles(file_path, smiles)"
                    )
                raise NotImplementedError(msg)

        elif isinstance(toolkit_registry, ToolkitWrapper):
            # TODO: Encapsulate this logic in ToolkitWrapper?
            toolkit = toolkit_registry
            if file_format not in toolkit.toolkit_file_read_formats:
                msg = (
                    f"Toolkit {toolkit.toolkit_name} can not read file {file_path} (format {file_format}). Supported "
                    f"formats for this toolkit are {toolkit.toolkit_file_read_formats}."
                )
                if toolkit.toolkit_name == "The RDKit" and file_format == "PDB":
                    msg += (
                        "RDKit can however read PDBs with a valid smiles string using the "
                        "Molecule.from_pdb_and_smiles(file_path, smiles) constructor"
                    )
                raise NotImplementedError(msg)
        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        mols = list()

        if isinstance(file_path, str):
            mols = toolkit.from_file(
                file_path,
                file_format=file_format,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
            )
        elif hasattr(file_path, "read"):
            file_obj = file_path
            mols = toolkit.from_file_obj(
                file_obj,
                file_format=file_format,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
            )

        if len(mols) == 0:
            raise Exception(f"Unable to read molecule from file: {file_path}")
        elif len(mols) == 1:
            return mols[0]

        return mols

    @classmethod
    def from_pdb(cls, file_path, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        .. deprecated:: 0.11.0
            ``from_pdb`` is deprecated and will soon be removed. Use
            :py:meth:`from_polymer_pdb` instead.
        """
        warnings.warn(
            "Molecule.from_pdb will soon be deprecated in favor of the more explicit "
            "Molecule.from_polymer_pdb",
            UserWarning,
        )
        return cls.from_polymer_pdb(file_path, toolkit_registry=toolkit_registry)

    @classmethod
    @requires_package("openmm")
    def from_polymer_pdb(
        cls, file_path: Union[str, TextIO], toolkit_registry=GLOBAL_TOOLKIT_REGISTRY
    ):
        """
        Loads a polymer from a PDB file.

        Currently only supports proteins with canonical amino acids that are
        either uncapped or capped by ACE/NME groups, but may later be extended
        to handle other common polymers, or accept user-defined polymer
        templates.

        Metadata such as residues, chains, and atom names are recorded in the
        ``Atom.properties`` attribute, which is a dictionary mapping from
        strings like "residue" to the appropriate value. ``from_polymer_pdb``
        returns a molecule that can be iterated over with the ``.residues`` and
        ``.chains`` attributes, as well as the usual ``.atoms``.

        This method proceeds in the following order:

        * Loads the polymer substructure template file (distributed with the
          OpenFF Toolkit)
        * Loads the PDB into an OpenMM :class:`openmm.app.PDBFile` object
        * Turns OpenMM topology into a temporarily invalid rdkit Molecule
        * Adds chemical information to the molecule:
            * For each substructure loaded from the substructure template file
                * Uses rdkit to find matches between the substructure and the
                  molecule
                * For any matches, assigns the atom formal charge and bond order
                  info from the substructure to the rdkit molecule, then marks
                  the atoms and bonds as having been assigned so they can not be
                  overwritten by subsequent isomorphisms
        * Take coordinates from the OpenMM Topology and add them as a conformer
        * Convert the rdkit Molecule to OpenFF

        Parameters
        ----------
        file_path : str or file object
            PDB information to be passed to OpenMM PDBFile object for loading
        toolkit_registry = ToolkitWrapper or ToolkitRegistry. Default = None
            Either a ToolkitRegistry, ToolkitWrapper

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
        """
        import openmm.unit as openmm_unit
        from openmm.app import PDBFile

        pdb = PDBFile(file_path)

        substructure_file_path = get_data_file_path(
            "proteins/aa_residues_substructures_explicit_bond_orders_with_caps.json"
        )

        with open(substructure_file_path, "r") as subfile:
            substructure_dictionary = json.load(subfile)

        offmol = toolkit_registry.call(
            "_polymer_openmm_topology_to_offmol", pdb.topology, substructure_dictionary
        )

        coords = unit.Quantity(
            np.array(
                [
                    [*vec3.value_in_unit(openmm_unit.angstrom)]
                    for vec3 in pdb.getPositions()
                ]
            ),
            unit.angstrom,
        )
        offmol.add_conformer(coords)
        offmol = toolkit_registry.call("_assign_aromaticity_and_stereo_from_3d", offmol)
        for i, atom in enumerate(pdb.topology.atoms()):
            offmol.atoms[i].name = atom.name
            offmol.atoms[i].metadata["residue_name"] = atom.residue.name
            offmol.atoms[i].metadata["residue_number"] = atom.residue.id
            offmol.atoms[i].metadata["chain_id"] = atom.residue.chain.id
        offmol.add_default_hierarchy_schemes()

        return offmol

    def _to_xyz_file(self, file_path):
        """
        Write the current molecule and its conformers to a multiframe xyz file, if the molecule
        has no current coordinates all atoms will be set to 0,0,0 in keeping with the behaviour of the
        backend toolkits.

        Information on the type of XYZ file written can be found here <http://openbabel.org/wiki/XYZ_(format)>.

        Parameters
        ----------
        file_path : str or file-like object
            A file-like object or the path to the file to be written.
        """

        # If we do not have a conformer make one with all zeros
        if self.n_conformers == 0:
            conformers = [
                unit.Quantity(np.zeros((self.n_atoms, 3), dtype=float), unit.angstrom)
            ]

        else:
            conformers = self._conformers

        if len(conformers) == 1:
            end = ""
            title = (
                lambda frame: f'{self.name if self.name != "" else self.hill_formula}{frame}\n'
            )
        else:
            end = 1
            title = (
                lambda frame: f'{self.name if self.name != "" else self.hill_formula} Frame {frame}\n'
            )

        # check if we have a file path or an open file object
        if isinstance(file_path, str):
            xyz_data = open(file_path, "w")
        else:
            xyz_data = file_path

        # add the data to the xyz_data list
        for i, geometry in enumerate(conformers, 1):
            xyz_data.write(f"{self.n_atoms}\n" + title(end))
            for j, atom_coords in enumerate(geometry.m_as(unit.angstrom)):
                x, y, z = atom_coords
                xyz_data.write(
                    f"{SYMBOLS[self.atoms[j].atomic_number]}       {x: .10f}   {y: .10f}   {z: .10f}\n"
                )

            # now we up the frame count
            end = i + 1

        # now close the file
        xyz_data.close()

    def to_file(self, file_path, file_format, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """Write the current molecule to a file or file-like object

        Parameters
        ----------
        file_path : str or file-like object
            A file-like object or the path to the file to be written.
        file_format : str
            Format specifier, one of ['MOL2', 'MOL2H', 'SDF', 'PDB', 'SMI', 'CAN', 'TDT']
            Note that not all toolkits support all formats
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for file writing. If a Toolkit is passed,
            only the highest-precedence toolkit is used

        Raises
        ------
        ValueError
            If the requested file_format is not supported by one of the installed cheminformatics toolkits

        Examples
        --------

        >>> molecule = Molecule.from_iupac('imatinib')
        >>> molecule.to_file('imatinib.mol2', file_format='mol2')  # doctest: +SKIP
        >>> molecule.to_file('imatinib.sdf', file_format='sdf')  # doctest: +SKIP
        >>> molecule.to_file('imatinib.pdb', file_format='pdb')  # doctest: +SKIP

        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            pass
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            toolkit_registry = ToolkitRegistry(toolkit_precedence=[])
            toolkit_registry.add_toolkit(toolkit)
        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        file_format = file_format.upper()
        # check if xyz, use the toolkit independent method.
        if file_format == "XYZ":
            return self._to_xyz_file(file_path=file_path)

        # Take the first toolkit that can write the desired output format
        toolkit = None
        for query_toolkit in toolkit_registry.registered_toolkits:
            if file_format in query_toolkit.toolkit_file_write_formats:
                toolkit = query_toolkit
                break

        # Raise an exception if no toolkit was found to provide the requested file_format
        if toolkit is None:
            supported_formats = {}
            for toolkit in toolkit_registry.registered_toolkits:
                supported_formats[
                    toolkit.toolkit_name
                ] = toolkit.toolkit_file_write_formats
            raise ValueError(
                f"The requested file format ({file_format}) is not available from any of the installed toolkits "
                f"(supported formats: {supported_formats})"
            )

        # Write file
        if type(file_path) == str:
            # Open file for writing
            toolkit.to_file(self, file_path, file_format)
        else:
            toolkit.to_file_obj(self, file_path, file_format)

    def enumerate_tautomers(
        self, max_states=20, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY
    ):
        """
        Enumerate the possible tautomers of the current molecule

        Parameters
        ----------
        max_states: int optional, default=20
            The maximum amount of molecules that should be returned

        toolkit_registry: openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use to enumerate the tautomers.

        Returns
        -------
        molecules: List[openff.toolkit.topology.Molecule]
            A list of openff.toolkit.topology.Molecule instances not including the input molecule.
        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            molecules = toolkit_registry.call(
                "enumerate_tautomers", molecule=self, max_states=max_states
            )

        elif isinstance(toolkit_registry, ToolkitWrapper):
            molecules = toolkit_registry.enumerate_tautomers(
                self, max_states=max_states
            )

        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        return molecules

    def enumerate_stereoisomers(
        self,
        undefined_only=False,
        max_isomers=20,
        rationalise=True,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Enumerate the stereocenters and bonds of the current molecule.

        Parameters
        ----------
        undefined_only: bool optional, default=False
            If we should enumerate all stereocenters and bonds or only those with undefined stereochemistry

        max_isomers: int optional, default=20
            The maximum amount of molecules that should be returned

        rationalise: bool optional, default=True
            If we should try to build and rationalise the molecule to ensure it can exist

        toolkit_registry: openff.toolkit.utils.toolkits.ToolkitRegistry or
            lopenff.toolkit.utils.toolkits.ToolkitWrapper, default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use to enumerate the stereoisomers.

        Returns
        --------
        molecules: List[openff.toolkit.topology.Molecule]
            A list of :class:`Molecule` instances not including the input molecule.

        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            molecules = toolkit_registry.call(
                "enumerate_stereoisomers",
                molecule=self,
                undefined_only=undefined_only,
                max_isomers=max_isomers,
                rationalise=rationalise,
            )

        elif isinstance(toolkit_registry, ToolkitWrapper):
            molecules = toolkit_registry.enumerate_stereoisomers(
                self,
                undefined_only=undefined_only,
                max_isomers=max_isomers,
                rationalise=rationalise,
            )

        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        return molecules

    @OpenEyeToolkitWrapper.requires_toolkit()
    def enumerate_protomers(self, max_states=10):
        """
        Enumerate the formal charges of a molecule to generate different protomoers.

        Parameters
        ----------
        max_states: int optional, default=10,
            The maximum number of protomer states to be returned.

        Returns
        -------
        molecules: List[openff.toolkit.topology.Molecule],
            A list of the protomers of the input molecules not including the input.
        """

        toolkit = OpenEyeToolkitWrapper()
        molecules = toolkit.enumerate_protomers(molecule=self, max_states=max_states)

        return molecules

    @classmethod
    @RDKitToolkitWrapper.requires_toolkit()
    def from_rdkit(
        cls, rdmol, allow_undefined_stereo=False, hydrogens_are_explicit=False
    ):
        """
        Create a Molecule from an RDKit molecule.

        Requires the RDKit to be installed.

        Parameters
        ----------
        rdmol : rkit.RDMol
            An RDKit molecule
        allow_undefined_stereo : bool, default=False
            If ``False``, raises an exception if ``rdmol`` contains undefined stereochemistry.
        hydrogens_are_explicit : bool, default=False
            If ``False``, RDKit will perform hydrogen addition using ``Chem.AddHs``

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An OpenFF molecule

        Examples
        --------

        Create a molecule from an RDKit molecule

        >>> from rdkit import Chem
        >>> from openff.toolkit.tests.utils import get_data_file_path
        >>> rdmol = Chem.MolFromMolFile(get_data_file_path('systems/monomers/ethanol.sdf'))
        >>> molecule = Molecule.from_rdkit(rdmol)

        """
        toolkit = RDKitToolkitWrapper()
        molecule = toolkit.from_rdkit(
            rdmol,
            allow_undefined_stereo=allow_undefined_stereo,
            hydrogens_are_explicit=hydrogens_are_explicit,
            _cls=cls,
        )
        return molecule

    def to_rdkit(
        self,
        aromaticity_model=DEFAULT_AROMATICITY_MODEL,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Create an RDKit molecule

        Requires the RDKit to be installed.

        Parameters
        ----------
        aromaticity_model : str, optional, default=DEFAULT_AROMATICITY_MODEL
            The aromaticity model to use

        Returns
        -------
        rdmol : rdkit.RDMol
            An RDKit molecule

        Examples
        --------

        Convert a molecule to RDKit

        >>> from openff.toolkit.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> rdmol = molecule.to_rdkit()

        """
        # toolkit = RDKitToolkitWrapper()
        if isinstance(toolkit_registry, ToolkitWrapper):
            return toolkit_registry.to_rdkit(self, aromaticity_model=aromaticity_model)
        else:
            return toolkit_registry.call(
                "to_rdkit", self, aromaticity_model=aromaticity_model
            )

    @classmethod
    @OpenEyeToolkitWrapper.requires_toolkit()
    def from_openeye(cls, oemol, allow_undefined_stereo=False):
        """
        Create a ``Molecule`` from an OpenEye molecule.

        Requires the OpenEye toolkit to be installed.

        Parameters
        ----------
        oemol : openeye.oechem.OEMol
            An OpenEye molecule
        allow_undefined_stereo : bool, default=False
            If ``False``, raises an exception if oemol contains undefined stereochemistry.

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An OpenFF molecule

        Examples
        --------

        Create a ``Molecule`` from an OpenEye OEMol

        >>> from openeye import oechem
        >>> from openff.toolkit.tests.utils import get_data_file_path
        >>> ifs = oechem.oemolistream(get_data_file_path('systems/monomers/ethanol.mol2'))
        >>> oemols = list(ifs.GetOEGraphMols())
        >>> molecule = Molecule.from_openeye(oemols[0])

        """
        toolkit = OpenEyeToolkitWrapper()
        molecule = toolkit.from_openeye(
            oemol, allow_undefined_stereo=allow_undefined_stereo, _cls=cls
        )
        return molecule

    @requires_package("qcelemental")
    def to_qcschema(self, multiplicity=1, conformer=0, extras=None):
        """
        Create a QCElemental Molecule.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        multiplicity : int, default=1,
            The multiplicity of the molecule;
            sets ``molecular_multiplicity`` field for QCElemental Molecule.

        conformer : int, default=0,
            The index of the conformer to use for the QCElemental Molecule geometry.

        extras : dict, default=None
            A dictionary that should be included in the ``extras`` field on the QCElemental Molecule.
            This can be used to include extra information, such as a smiles representation.

        Returns
        ---------
        qcelemental.models.Molecule
            A validated QCElemental Molecule.

        Examples
        --------

        Create a QCElemental Molecule:

        >>> import qcelemental as qcel
        >>> mol = Molecule.from_smiles('CC')
        >>> mol.generate_conformers(n_conformers=1)
        >>> qcemol = mol.to_qcschema()

        Raises
        --------
        MissingOptionalDependencyError
            If qcelemental is not installed, the qcschema can not be validated.
        InvalidConformerError
            No conformer found at the given index.

        """
        import qcelemental as qcel

        # get/ check the geometry
        try:
            geometry = self.conformers[conformer].m_as(unit.bohr)
        except (IndexError, TypeError):
            raise InvalidConformerError(
                "The molecule must have a conformation to produce a valid qcschema; "
                f"no conformer was found at index {conformer}."
            )

        # Gather the required qcschema data
        charge = self.total_charge.m_as(unit.elementary_charge)
        connectivity = [
            (bond.atom1_index, bond.atom2_index, bond.bond_order) for bond in self.bonds
        ]
        symbols = [SYMBOLS[atom.atomic_number] for atom in self.atoms]
        if extras is not None:
            extras[
                "canonical_isomeric_explicit_hydrogen_mapped_smiles"
            ] = self.to_smiles(mapped=True)
        else:
            extras = {
                "canonical_isomeric_explicit_hydrogen_mapped_smiles": self.to_smiles(
                    mapped=True
                )
            }

        schema_dict = {
            "symbols": symbols,
            "geometry": geometry,
            # If we have no bonds we must supply None
            "connectivity": connectivity if connectivity else None,
            "molecular_charge": charge,
            "molecular_multiplicity": multiplicity,
            "extras": extras,
        }

        return qcel.models.Molecule.from_data(schema_dict, validate=True)

    @classmethod
    def from_mapped_smiles(
        cls,
        mapped_smiles,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
    ):
        """
        Create an :class:`Molecule` from a mapped SMILES made with cmiles.
        The molecule will be in the order of the indexing in the mapped smiles string.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        mapped_smiles: str
            A CMILES-style mapped smiles string with explicit hydrogens.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapper, optional
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion

        allow_undefined_stereo : bool, default=False
            If false, raises an exception if oemol contains undefined stereochemistry.

        Returns
        ----------
        offmol : openff.toolkit.topology.molecule.Molecule
            An OpenFF molecule instance.

        Raises
        --------
        SmilesParsingError
            If the given SMILES had no indexing picked up by the toolkits.
        """

        # create the molecule from the smiles and check we have the right number of indexes
        # in the mapped SMILES
        offmol = cls.from_smiles(
            mapped_smiles,
            hydrogens_are_explicit=True,
            toolkit_registry=toolkit_registry,
            allow_undefined_stereo=allow_undefined_stereo,
        )

        # check we found some mapping and remove it as we do not want to expose atom maps
        try:
            mapping = offmol._properties.pop("atom_map")
        except KeyError:
            raise SmilesParsingError(
                "The given SMILES has no indexing, please generate a valid explicit hydrogen "
                "mapped SMILES using cmiles."
            )

        if len(mapping) != offmol.n_atoms:
            raise SmilesParsingError(
                "The mapped smiles does not contain enough indexes to remap the molecule."
            )

        # remap the molecule using the atom map found in the smiles
        # the order is mapping = Dict[current_index: new_index]
        # first renumber the mapping dict indexed from 0, currently from 1 as 0 indicates no mapping in toolkits
        adjusted_mapping = dict((current, new - 1) for current, new in mapping.items())

        return offmol.remap(adjusted_mapping, current_to_new=True)

    @classmethod
    @requires_package("qcelemental")
    def from_qcschema(
        cls,
        qca_record,
        client=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo=False,
    ):
        """
        Create a Molecule from a QCArchive molecule record or dataset entry
        based on attached cmiles information.

        For a molecule record, a conformer will be set from its geometry.

        For a dataset entry, if a corresponding client instance is provided,
        the starting geometry for that entry will be used as a conformer.

        A QCElemental Molecule produced from ``Molecule.to_qcschema`` can be round-tripped
        through this method to produce a new, valid Molecule.

        Parameters
        ----------
        qca_record : dict
            A QCArchive molecule record or dataset entry.

        client : optional, default=None,
            A qcportal.FractalClient instance to use for fetching an initial geometry.
            Only used if ``qca_record`` is a dataset entry.

        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or
            openff.toolkit.utils.toolkits.ToolkitWrapper, optional
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMILES-to-molecule conversion

        allow_undefined_stereo : bool, default=False
            If false, raises an exception if qca_record contains undefined stereochemistry.

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An OpenFF molecule instance.

        Examples
        --------
        Get Molecule from a QCArchive molecule record:

        >>> from qcportal import FractalClient
        >>> client = FractalClient()
        >>> offmol = Molecule.from_qcschema(client.query_molecules(molecular_formula="C16H20N3O5")[0])

        Get Molecule from a QCArchive optimization entry:

        >>> from qcportal import FractalClient
        >>> client = FractalClient()
        >>> optds = client.get_collection("OptimizationDataset",
                                          "SMIRNOFF Coverage Set 1")
        >>> offmol = Molecule.from_qcschema(optds.get_entry('coc(o)oc-0'))

        Same as above, but with conformer(s) from initial molecule(s) by providing client to database:

        >>> offmol = Molecule.from_qcschema(optds.get_entry('coc(o)oc-0'), client=client)

        Raises
        -------
        AttributeError
            - If the record dict can not be made from ``qca_record``.
            - If a ``client`` is passed and it could not retrieve the initial molecule.

        KeyError
            If the dict does not contain the ``canonical_isomeric_explicit_hydrogen_mapped_smiles``.

        InvalidConformerError
            Silent error, if the conformer could not be attached.
        """

        # We can accept the Dataset entry record or the dict with JSON encoding
        # lets get it all in the dict rep
        if not isinstance(qca_record, dict):
            try:
                qca_record = qca_record.dict(encoding="json")
            except AttributeError:
                raise AttributeError(
                    "The object passed could not be converted to a dict with json encoding"
                )

        # identify if this is a dataset entry
        if "attributes" in qca_record:
            mapped_smiles = qca_record["attributes"][
                "canonical_isomeric_explicit_hydrogen_mapped_smiles"
            ]
            if client is not None:
                # try and find the initial molecule conformations and attach them
                # collect the input molecules
                try:
                    input_mols = client.query_molecules(
                        id=qca_record["initial_molecules"]
                    )
                except KeyError:
                    # this must be an optimisation record
                    input_mols = client.query_molecules(
                        id=qca_record["initial_molecule"]
                    )
                except AttributeError:
                    raise AttributeError(
                        "The provided client can not query molecules, make sure it is an instance of"
                        "qcportal.client.FractalClient() with the correct address."
                    )
            else:
                input_mols = []

        # identify if this is a molecule record
        elif "extras" in qca_record:
            mapped_smiles = qca_record["extras"][
                "canonical_isomeric_explicit_hydrogen_mapped_smiles"
            ]
            input_mols = [qca_record]
        else:
            raise KeyError(
                "The record must contain the hydrogen mapped smiles to be safely made from the archive. "
                "It is not present in either 'attributes' or 'extras' on the provided `qca_record`"
            )

        # make a new molecule that has been reordered to match the cmiles mapping
        offmol = cls.from_mapped_smiles(
            mapped_smiles,
            toolkit_registry=toolkit_registry,
            allow_undefined_stereo=allow_undefined_stereo,
        )

        # now for each molecule convert and attach the input geometry
        initial_ids = {}
        for molecule in input_mols:
            if not isinstance(molecule, dict):
                mol = molecule.dict(encoding="json")
            else:
                mol = molecule

            geometry = unit.Quantity(
                np.array(mol["geometry"], float).reshape(-1, 3), unit.bohr
            )
            try:
                offmol._add_conformer(geometry.to(unit.angstrom))
                # in case this molecule didn't come from a server at all
                if "id" in mol:
                    initial_ids[mol["id"]] = offmol.n_conformers - 1
            except InvalidConformerError:
                print(
                    "Invalid conformer for this molecule, the geometry could not be attached."
                )

        # attach a dict that has the initial molecule ids and the number of the conformer it is stored in
        # if it's empty, don't bother
        if initial_ids:
            offmol._properties["initial_molecules"] = initial_ids

        return offmol

    @classmethod
    @RDKitToolkitWrapper.requires_toolkit()
    def from_pdb_and_smiles(cls, file_path, smiles, allow_undefined_stereo=False):
        """
        Create a Molecule from a pdb file and a SMILES string using RDKit.

        Requires RDKit to be installed.

        .. warning :: This API is experimental and subject to change.

        The molecule is created and sanitised based on the SMILES string, we then find a mapping
        between this molecule and one from the PDB based only on atomic number and connections.
        The SMILES molecule is then reindexed to match the PDB, the conformer is attached, and the
        molecule returned.

        Note that any stereochemistry in the molecule is set by the SMILES, and not the coordinates
        of the PDB.

        Parameters
        ----------
        file_path: str
            PDB file path
        smiles : str
            a valid smiles string for the pdb, used for stereochemistry, formal charges, and bond order
        allow_undefined_stereo : bool, default=False
            If false, raises an exception if SMILES contains undefined stereochemistry.

        Returns
        --------
        molecule : openff.toolkit.Molecule
            An OFFMol instance with ordering the same as used in the PDB file.

        Raises
        ------
        InvalidConformerError
            If the SMILES and PDB molecules are not isomorphic.
        """

        toolkit = RDKitToolkitWrapper()
        return toolkit.from_pdb_and_smiles(
            file_path, smiles, allow_undefined_stereo, _cls=cls
        )

    def canonical_order_atoms(self, toolkit_registry=GLOBAL_TOOLKIT_REGISTRY):
        """
        Produce a copy of the molecule with the atoms reordered canonically.

        Each toolkit defines its own canonical ordering of atoms. The canonical
        order may change from toolkit version to toolkit version or between
        toolkits.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        toolkit_registry : openff.toolkit.utils.toolkits.ToolkitRegistry or
            openff.toolkit.utils.toolkits.ToolkitWrapper, optional
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for
            SMILES-to-molecule conversion

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
            An new OpenFF style molecule with atoms in the canonical order.
        """

        if isinstance(toolkit_registry, ToolkitRegistry):
            return toolkit_registry.call("canonical_order_atoms", self)
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            return toolkit.canonical_order_atoms(self)
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to from_smiles. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}."
            )

    def remap(self, mapping_dict, current_to_new=True):
        """
        Remap all of the indexes in the molecule to match the given mapping dict

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        mapping_dict : dict,
            A dictionary of the mapping between indexes, this should start from 0.
        current_to_new : bool, default=True
            If this is ``True``, then ``mapping_dict`` is of the form ``{current_index: new_index}``;
            otherwise, it is of the form ``{new_index: current_index}``

        Returns
        -------
        new_molecule :  openff.toolkit.topology.molecule.Molecule
            An openff.toolkit.Molecule instance with all attributes transferred, in the PDB order.
        """

        # make sure the size of the mapping matches the current molecule
        if len(mapping_dict) != self.n_atoms:
            raise ValueError(
                f"The number of mapping indices({len(mapping_dict)}) does not match the number of"
                f"atoms in this molecule({self.n_atoms})"
            )

        # make two mapping dicts we need new to old for atoms
        # and old to new for bonds
        if current_to_new:
            cur_to_new = mapping_dict
            new_to_cur = dict(zip(mapping_dict.values(), mapping_dict.keys()))
        else:
            new_to_cur = mapping_dict
            cur_to_new = dict(zip(mapping_dict.values(), mapping_dict.keys()))

        new_molecule = self.__class__()
        new_molecule.name = self.name

        try:
            # add the atoms list
            for i in range(self.n_atoms):
                # get the old atom info
                old_atom = self._atoms[new_to_cur[i]]
                new_molecule._add_atom(**old_atom.to_dict())
        # this is the first time we access the mapping; catch an index error here corresponding to mapping that starts
        # from 0 or higher
        except (KeyError, IndexError):
            raise IndexError(
                f"The mapping supplied is missing a relation corresponding to atom({i})"
            )

        # add the bonds but with atom indexes in a sorted ascending order
        for bond in self._bonds:
            atoms = sorted([cur_to_new[bond.atom1_index], cur_to_new[bond.atom2_index]])
            bond_dict = bond.to_dict()
            bond_dict["atom1"] = atoms[0]
            bond_dict["atom2"] = atoms[1]
            new_molecule._add_bond(**bond_dict)

        # we can now resort the bonds
        sorted_bonds = sorted(
            new_molecule.bonds, key=operator.attrgetter("atom1_index", "atom2_index")
        )
        new_molecule._bonds = sorted_bonds

        # remap the charges
        if self.partial_charges is not None:
            new_charges = np.zeros(self.n_atoms)
            for i in range(self.n_atoms):
                new_charges[i] = self.partial_charges[new_to_cur[i]].m_as(
                    unit.elementary_charge
                )
            new_molecule.partial_charges = new_charges * unit.elementary_charge

        # remap the conformers there can be more than one
        if self.conformers is not None:
            for conformer in self.conformers:
                new_conformer = np.zeros((self.n_atoms, 3))
                for i in range(self.n_atoms):
                    new_conformer[i] = conformer[new_to_cur[i]].m_as(unit.angstrom)
                new_molecule._add_conformer(new_conformer * unit.angstrom)

        # move any properties across
        new_molecule._properties = self._properties

        return new_molecule

    def to_openeye(
        self,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        aromaticity_model=DEFAULT_AROMATICITY_MODEL,
    ):
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
        # toolkit = OpenEyeToolkitWrapper()
        if isinstance(toolkit_registry, ToolkitWrapper):
            return toolkit_registry.to_openeye(
                self, aromaticity_model=aromaticity_model
            )
        else:
            return toolkit_registry.call(
                "to_openeye", self, aromaticity_model=aromaticity_model
            )

    def _construct_angles(self):
        """
        Get an iterator over all i-j-k angles.
        """
        # TODO: Build Angle objects instead of tuple of atoms.
        if not hasattr(self, "_angles"):
            self._construct_bonded_atoms_list()
            self._angles = set()
            for atom1 in self._atoms:
                for atom2 in self._bondedAtoms[atom1]:
                    for atom3 in self._bondedAtoms[atom2]:
                        if atom1 == atom3:
                            continue
                        # TODO: Encapsulate this logic into an Angle class.
                        if atom1.molecule_atom_index < atom3.molecule_atom_index:
                            self._angles.add((atom1, atom2, atom3))
                        else:
                            self._angles.add((atom3, atom2, atom1))

    def _construct_torsions(self):
        """
        Construct sets containing the atoms improper and proper torsions
        """
        # TODO: Build Proper/ImproperTorsion objects instead of tuple of atoms.
        if not hasattr(self, "_torsions"):
            self._construct_bonded_atoms_list()

            self._propers = set()
            self._impropers = set()
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

                            if atom1.molecule_atom_index < atom4.molecule_atom_index:
                                torsion = (atom1, atom2, atom3, atom4)
                            else:
                                torsion = (atom4, atom3, atom2, atom1)

                            self._propers.add(torsion)

                        for atom3i in self._bondedAtoms[atom2]:
                            if atom3i == atom3:
                                continue
                            if atom3i == atom1:
                                continue

                            improper = (atom1, atom2, atom3, atom3i)
                            self._impropers.add(improper)

            self._torsions = self._propers | self._impropers

    def _construct_bonded_atoms_list(self):
        """
        Construct list of all atoms each atom is bonded to.

        """
        # TODO: Add this to cached_properties
        if not hasattr(self, "_bondedAtoms"):
            # self._atoms = [ atom for atom in self.atoms() ]
            self._bondedAtoms = dict()
            for atom in self._atoms:
                self._bondedAtoms[atom] = set()
            for bond in self._bonds:
                atom1 = self.atoms[bond.atom1_index]
                atom2 = self.atoms[bond.atom2_index]
                self._bondedAtoms[atom1].add(atom2)
                self._bondedAtoms[atom2].add(atom1)

    def _is_bonded(self, atom_index_1, atom_index_2):
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


        """
        self._construct_bonded_atoms_list()
        atom1 = self._atoms[atom_index_1]
        atom2 = self._atoms[atom_index_2]
        return atom2 in self._bondedAtoms[atom1]

    def get_bond_between(self, i, j):
        """Returns the bond between two atoms

        Parameters
        ----------
        i, j : int or Atom
            Atoms or atom indices to check

        Returns
        -------
        bond : Bond
            The bond between i and j.

        """
        if isinstance(i, int) and isinstance(j, int):
            atom_i = self._atoms[i]
            atom_j = self._atoms[j]
        elif isinstance(i, Atom) and isinstance(j, Atom):
            atom_i = i
            atom_j = j
        else:
            raise TypeError(
                "Invalid input passed to get_bond_between(). Expected ints or Atoms, "
                f"got {j} and {j}."
            )

        for bond in atom_i.bonds:

            for atom in bond.atoms:

                if atom == atom_i:
                    continue

                if atom == atom_j:
                    return bond

        from openff.toolkit.topology import NotBondedError

        raise NotBondedError(f"No bond between atom {i} and {j}")


class Molecule(FrozenMolecule):
    """
    Mutable chemical representation of a molecule, such as a small molecule or biopolymer.

    .. todo :: What other API calls would be useful for supporting biopolymers as small molecules? Perhaps iterating
        over chains and residues?

    Examples
    --------

    Create a molecule from an sdf file

    >>> from openff.toolkit.utils import get_data_file_path
    >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
    >>> molecule = Molecule(sdf_filepath)

    Convert to OpenEye OEMol object

    >>> oemol = molecule.to_openeye()

    Create a molecule from an OpenEye molecule

    >>> molecule = Molecule.from_openeye(oemol)

    Convert to RDKit Mol object

    >>> rdmol = molecule.to_rdkit()

    Create a molecule from an RDKit molecule

    >>> molecule = Molecule.from_rdkit(rdmol)

    Create a molecule from IUPAC name (requires the OpenEye toolkit)

    >>> molecule = Molecule.from_iupac('imatinib')

    Create a molecule from SMILES

    >>> molecule = Molecule.from_smiles('Cc1ccccc1')

    .. warning :: This API is experimental and subject to change.

    """

    def __init__(self, *args, **kwargs):
        """
        Create a new Molecule object

        Parameters
        ----------
        other : optional, default=None
            If specified, attempt to construct a copy of the molecule from the
            specified object. This can be any one of the following:

            * a :class:`Molecule` object
            * a file that can be used to construct a :class:`Molecule` object
            * an ``openeye.oechem.OEMol``
            * an ``rdkit.Chem.rdchem.Mol``
            * a serialized :class:`Molecule` object

        Examples
        --------

        Create an empty molecule:

        >>> empty_molecule = Molecule()

        Create a molecule from a file that can be used to construct a molecule,
        using either a filename or file-like object:

        >>> from openff.toolkit.utils import get_data_file_path
        >>> sdf_filepath = get_data_file_path('molecules/ethanol.sdf')
        >>> molecule = Molecule(sdf_filepath)
        >>> molecule = Molecule(open(sdf_filepath, 'r'), file_format='sdf')

        >>> import gzip
        >>> mol2_gz_filepath = get_data_file_path('molecules/toluene.mol2.gz')
        >>> molecule = Molecule(gzip.GzipFile(mol2_gz_filepath, 'r'), file_format='mol2')

        Create a molecule from another molecule:

        >>> molecule_copy = Molecule(molecule)

        Convert to OpenEye OEMol object

        >>> oemol = molecule.to_openeye()

        Create a molecule from an OpenEye molecule:

        >>> molecule = Molecule(oemol)

        Convert to RDKit Mol object

        >>> rdmol = molecule.to_rdkit()

        Create a molecule from an RDKit molecule:

        >>> molecule = Molecule(rdmol)

        Convert the molecule into a dictionary and back again:

        >>> serialized_molecule = molecule.to_dict()
        >>> molecule_copy = Molecule(serialized_molecule)

        .. todo ::

           * If a filename or file-like object is specified but the file
             contains more than one molecule, what is the proper behavior?
             Read just the first molecule, or raise an exception if more
             than one molecule is found?

           * Should we also support SMILES strings or IUPAC names for
             ``other``?

        """
        super(Molecule, self).__init__(*args, **kwargs)

    # TODO: Change this to add_atom(Atom) to improve encapsulation and extensibility?
    def add_atom(
        self,
        atomic_number,
        formal_charge,
        is_aromatic,
        stereochemistry=None,
        name=None,
        metadata=None,
    ):
        """
        Add an atom to the molecule.

        Parameters
        ----------
        atomic_number : int
            Atomic number of the atom
        formal_charge : int
            Formal charge of the atom
        is_aromatic : bool
            If ``True``, atom is aromatic; if ``False``, not aromatic
        stereochemistry : str, optional, default=None
            Either ``'R'`` or ``'S'`` for specified stereochemistry, or ``None`` if stereochemistry is irrelevant
        name : str, optional
            An optional name for the atom
        metadata : dict[str: (int, str)], default=None
            An optional dictionary where keys are strings and values are strings or ints. This is intended
            to record atom-level information used to inform hierarchy definition and iteration, such as
            grouping atom by residue and chain.

        Returns
        -------
        index : int
            The index of the atom in the molecule

        Examples
        --------

        Define a methane molecule

        >>> molecule = Molecule()
        >>> molecule.name = 'methane'
        >>> C = molecule.add_atom(6, 0, False)
        >>> H1 = molecule.add_atom(1, 0, False)
        >>> H2 = molecule.add_atom(1, 0, False)
        >>> H3 = molecule.add_atom(1, 0, False)
        >>> H4 = molecule.add_atom(1, 0, False)
        >>> bond_idx = molecule.add_bond(C, H1, False, 1)
        >>> bond_idx = molecule.add_bond(C, H2, False, 1)
        >>> bond_idx = molecule.add_bond(C, H3, False, 1)
        >>> bond_idx = molecule.add_bond(C, H4, False, 1)

        """
        atom_index = self._add_atom(
            atomic_number,
            formal_charge,
            is_aromatic,
            stereochemistry=stereochemistry,
            name=name,
            metadata=metadata,
        )
        return atom_index

    def add_bond(
        self,
        atom1,
        atom2,
        bond_order,
        is_aromatic,
        stereochemistry=None,
        fractional_bond_order=None,
    ):
        """
        Add a bond between two specified atom indices


        Parameters
        ----------
        atom1 : int or openff.toolkit.topology.molecule.Atom
            Index of first atom
        atom2 : int or openff.toolkit.topology.molecule.Atom
            Index of second atom
        bond_order : int
            Integral bond order of Kekulized form
        is_aromatic : bool
            True if this bond is aromatic, False otherwise
        stereochemistry : str, optional, default=None
            Either ``'E'`` or ``'Z'`` for specified stereochemistry, or ``None`` if stereochemistry is irrelevant
        fractional_bond_order : float, optional, default=None
            The fractional (eg. Wiberg) bond order

        Returns
        -------
        index: int
            Index of the bond in this molecule

        """
        bond_index = self._add_bond(
            atom1,
            atom2,
            bond_order,
            is_aromatic,
            stereochemistry=stereochemistry,
            fractional_bond_order=fractional_bond_order,
        )
        return bond_index

    def add_conformer(self, coordinates):
        """
        Add a conformation of the molecule

        Parameters
        ----------
        coordinates: unit-wrapped np.array with shape (n_atoms, 3) and dimension of distance
            Coordinates of the new conformer, with the first dimension of the array corresponding to the atom index in
            the molecule's indexing system.

        Returns
        -------
        index: int
            The index of this conformer
        """

        # TODO how can be check that a set of coords and no connections
        #   is a conformation that does not change connectivity?

        return self._add_conformer(coordinates)

    def visualize(
        self,
        backend="rdkit",
        width=None,
        height=None,
        show_all_hydrogens=True,
    ):
        """
        Render a visualization of the molecule in Jupyter

        Parameters
        ----------
        backend : str, optional, default='rdkit'
            The visualization engine to use. Choose from:

            - ``"rdkit"``
            - ``"openeye"``
            - ``"nglview"`` (requires conformers)

        width : int, optional, default=500
            Width of the generated representation (only applicable to
            ``backend=openeye`` or ``backend=rdkit``)
        height : int, optional, default=300
            Width of the generated representation (only applicable to
            ``backend=openeye`` or ``backend=rdkit``)
        show_all_hydrogens : bool, optional, default=True
            Whether to explicitly depict all hydrogen atoms. (only applicable to
            ``backend=openeye`` or ``backend=rdkit``)

        Returns
        -------
        object
            Depending on the backend chosen:

            - rdkit → IPython.display.SVG
            - openeye → IPython.display.Image
            - nglview → nglview.NGLWidget

        """
        from openff.toolkit.utils.toolkits import OPENEYE_AVAILABLE, RDKIT_AVAILABLE

        backend = backend.lower()

        if backend == "nglview":
            try:
                import nglview as nv
            except ImportError:
                raise MissingOptionalDependencyError("nglview")

            if width is not None or height is not None:
                # TODO: More specific exception
                raise ValueError(
                    "The width, height, and show_all_hydrogens arguments do not apply to the nglview backend."
                )
            elif not show_all_hydrogens:
                # TODO: More specific exception
                # TODO: Implement this? Should be able to just strip hydrogens from the PDB
                raise ValueError(
                    "show_all_hydrogens=False is not supported by the nglview backend"
                )

            if self.conformers:
                from openff.toolkit.utils.viz import _OFFTrajectoryNGLView

                trajectory_like = _OFFTrajectoryNGLView(self)
                widget = nv.NGLWidget(trajectory_like)
                return widget
            else:
                # TODO: More specific exception
                raise ValueError(
                    "Visualizing with NGLview requires that the molecule has "
                    "conformers."
                )

        width = 500 if width is None else width
        height = 300 if height is None else height
        show_all_hydrogens = True if show_all_hydrogens is None else show_all_hydrogens

        if backend == "rdkit":
            if RDKIT_AVAILABLE:
                from IPython.display import SVG
                from rdkit.Chem.Draw import (  # type: ignore[import]
                    rdDepictor,
                    rdMolDraw2D,
                )
                from rdkit.Chem.rdmolops import RemoveHs  # type: ignore[import]

                rdmol = self.to_rdkit()

                if not show_all_hydrogens:
                    # updateExplicitCount: Keep a record of the hydrogens we remove.
                    # This is used in visualization to distinguish eg radicals from normal species
                    rdmol = RemoveHs(rdmol, updateExplicitCount=True)

                rdDepictor.SetPreferCoordGen(True)
                rdDepictor.Compute2DCoords(rdmol)
                rdmol = rdMolDraw2D.PrepareMolForDrawing(rdmol)

                drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
                drawer.DrawMolecule(rdmol)
                drawer.FinishDrawing()

                return SVG(drawer.GetDrawingText())
            else:
                warnings.warn(
                    "RDKit was requested as a visualization backend but "
                    "it was not found to be installed. Falling back to "
                    "trying to use OpenEye for visualization."
                )
                backend = "openeye"
        if backend == "openeye":
            if OPENEYE_AVAILABLE:
                from IPython.display import Image
                from openeye import oedepict

                oemol = self.to_openeye()

                opts = oedepict.OE2DMolDisplayOptions(
                    width, height, oedepict.OEScale_AutoScale
                )

                if show_all_hydrogens:
                    opts.SetHydrogenStyle(oedepict.OEHydrogenStyle_ImplicitAll)

                oedepict.OEPrepareDepiction(oemol)
                img = oedepict.OEImage(width, height)
                display = oedepict.OE2DMolDisplay(oemol, opts)
                oedepict.OERenderMolecule(img, display)
                png = oedepict.OEWriteImageToString("png", img)
                return Image(png)

        # TODO: More specific exception
        raise ValueError("Could not find an appropriate backend")

    def perceive_residues(self, substructure_file_path=None, strict_chirality=True):
        """
        Perceive a polymer's residues and permit iterating over them.

        Perceives residues by matching substructures in the current molecule
        with a substructure dictionary file, using SMARTS, and assigns residue
        names and numbers to atom metadata. It then constructs a residue hierarchy
        scheme to allow iterating over residues.


        Parameters
        ----------
        substructure_file_path : str, optional, default=None
            Path to substructure library file in JSON format. Defaults to using
            built-in substructure file.
        strict_chirality: bool, optional, default=True
            Whether to use strict chirality symbols (stereomarks) for
            substructure matchings with SMARTS.
        """
        # Read substructure dictionary file
        if not substructure_file_path:
            substructure_file_path = get_data_file_path(
                "proteins/aa_residues_substructures_with_caps.json"
            )
        with open(substructure_file_path, "r") as subfile:
            substructure_dictionary = json.load(subfile)

        # TODO: Think of a better way to deal with no strict chirality case
        # if ignoring strict chirality, remove/update keys in inner dictionary
        if not strict_chirality:
            # make a copy of substructure dict
            substructure_dictionary_no_chirality = deepcopy(substructure_dictionary)
            # Update inner key (SMARTS) maintaining its value
            for res_name, inner_dict in substructure_dictionary.items():
                for smarts, atom_types in inner_dict.items():
                    smarts_no_chirality = smarts.replace("@", "")  # remove @ in smarts
                    substructure_dictionary_no_chirality[res_name][
                        smarts_no_chirality
                    ] = substructure_dictionary_no_chirality[res_name].pop(
                        smarts
                    )  # update key
            # replace with the new substructure dictionary
            substructure_dictionary = substructure_dictionary_no_chirality

        all_matches = list()
        for residue_name, smarts_dict in substructure_dictionary.items():
            matches = dict()
            for smarts in smarts_dict:
                for match in self.chemical_environment_matches(smarts):
                    matches[match] = smarts
                    all_matches.append(
                        {
                            "atom_idxs": match,
                            "atom_idxs_set": set(match),
                            "smarts": smarts,
                            "residue_name": residue_name,
                            "atom_names": smarts_dict[smarts],
                        }
                    )

        # Remove matches that are subsets of other matches
        # give precedence to the SMARTS defined at the end of the file
        match_idxs_to_delete = set()
        for match_idx in range(len(all_matches) - 1, 0, -1):
            this_match_set = all_matches[match_idx]["atom_idxs_set"]
            this_match_set_size = len(this_match_set)
            for match_before_this_idx in range(match_idx):
                match_before_this_set = all_matches[match_before_this_idx][
                    "atom_idxs_set"
                ]
                match_before_this_set_size = len(match_before_this_set)
                n_overlapping_atoms = len(
                    this_match_set.intersection(match_before_this_set)
                )
                if n_overlapping_atoms > 0:
                    if match_before_this_set_size < this_match_set_size:
                        match_idxs_to_delete.add(match_before_this_idx)
                    else:
                        match_idxs_to_delete.add(match_idx)

        match_idxs_to_delete_list = sorted(list(match_idxs_to_delete), reverse=True)
        for match_idx in match_idxs_to_delete_list:
            all_matches.pop(match_idx)

        all_matches.sort(key=lambda x: min(x["atom_idxs"]))

        # Now the matches have been deduplicated and de-subsetted
        for residue_num, match_dict in enumerate(all_matches):
            for smarts_idx, atom_idx in enumerate(match_dict["atom_idxs"]):
                self.atoms[atom_idx].metadata["residue_name"] = match_dict[
                    "residue_name"
                ]
                self.atoms[atom_idx].metadata["residue_number"] = residue_num + 1
                self.atoms[atom_idx].metadata["atom_name"] = match_dict["atom_names"][
                    smarts_idx
                ]

        # Now add the residue hierarchy scheme
        self._add_residue_hierarchy_scheme()

    def _ipython_display_(self):  # pragma: no cover
        from IPython.display import display

        try:
            return display(self.visualize(backend="nglview"))
        except (ImportError, ValueError):
            pass

        try:
            return display(self.visualize(backend="rdkit"))
        except ValueError:
            pass

        try:
            return display(self.visualize(backend="openeye"))
        except ValueError:
            pass


def _networkx_graph_to_hill_formula(graph: "nx.Graph") -> str:
    """
    Convert a NetworkX graph to a Hill formula.

    Parameters
    ----------
    graph : nx.Graph
        The graph to convert.

    Returns
    -------
    str
        The Hill formula corresponding to the graph.

    """
    import networkx as nx

    if not isinstance(graph, nx.Graph):
        raise Exception("The graph must be a NetworkX graph.")

    atom_nums = list(dict(graph.nodes(data="atomic_number", default=1)).values())
    return _atom_nums_to_hill_formula(atom_nums)


def _atom_nums_to_hill_formula(atom_nums: List[int]) -> str:
    """
    Given a `Counter` object of atom counts by atomic number, generate the corresponding
    Hill formula. See https://en.wikipedia.org/wiki/Chemical_formula#Hill_system"""
    from collections import Counter

    SYMBOLS_ = deepcopy(SYMBOLS)
    SYMBOLS_[0] = "X"

    atom_symbol_counts = Counter(SYMBOLS_[atom_num] for atom_num in atom_nums)

    formula = []
    # Check for C and H first, to make a correct hill formula
    for el in ["C", "H"]:
        if el in atom_symbol_counts:
            count = atom_symbol_counts.pop(el)
            formula.append(el)
            if count > 1:
                formula.append(str(count))

    # now get the rest of the elements in alphabetical ordering
    for el in sorted(atom_symbol_counts.keys()):
        count = atom_symbol_counts.pop(el)
        formula.append(el)
        if count > 1:
            formula.append(str(count))

    return "".join(formula)


def _nth_degree_neighbors_from_graphlike(
    graphlike: Union[Molecule, "_SimpleMolecule"], n_degrees: int
) -> Generator[
    Union[Tuple[Atom, Atom], Tuple["_SimpleAtom", "_SimpleAtom"]], None, None
]:
    """
    Given a graph-like object, return a tuple of the nth degree neighbors of each atom.

    The input `graphlike` object must provide a .to_networkx() method and an
    `atoms` property that can be indexed.

    See Molecule.nth_degree_neighbors for more details.

    Parameters
    ----------
    graphlike : Union[Molecule, _SimpleMolecule]
        The graph-like object to get the neighbors of.
    n: int
        The number of bonds separating atoms in each pair

    Returns
    -------
    neighbors: iterator of tuple of Atom
        Tuples (len 2) of atom that are separated by ``n`` bonds.
    """
    graph = graphlike.to_networkx()

    for node_i in graph.nodes:
        for node_j in graph.nodes:
            if node_i == node_j:
                continue

            path_length = nx.shortest_path_length(graph, node_i, node_j)

            if path_length == n_degrees:
                if node_i > node_j:
                    continue
                yield (graphlike.atoms[node_i], graphlike.atoms[node_j])


class HierarchyScheme:
    """
    Perceives hierarchy elements from the metadata of atoms in a ``Molecule``.

    The Open Force Field Toolkit has no native understanding of hierarchical
    atom organisation schemes common to other biomolecular software, such as
    "residues" or "chains" (see :ref:`userguide_hierarchy`). To facilitate
    iterating over groups of atoms, a ``HierarchyScheme`` can be used to collect
    atoms into ``HierarchyElements``, groups of atoms that share the same
    values for certain metadata elements. Metadata elements are stored in the
    ``Atom.properties`` attribute.

    Hierarchy schemes are not updated dynamically; if a ``Molecule`` with
    hierarchy schemes changes, :meth:`Molecule.update_hierarchy_schemes()` must
    be called before the scheme is iterated over again or else the grouping
    may be incorrect.

    A ``HierarchyScheme`` contains the information needed to perceive
    ``HierarchyElement`` objects from a ``Molecule`` containing atoms with
    metadata.
    """

    def __init__(
        self,
        parent: FrozenMolecule,
        uniqueness_criteria: Union[Tuple[str], List[str]],
        iterator_name: str,
    ):
        """
        Create a new hierarchy scheme for iterating over groups of atoms.

        Parameters
        ----------

        parent
            The ``Molecule`` to which this scheme belongs.
        uniqueness_criteria
            The names of ``Atom`` metadata entries that define this scheme. An
            atom belongs to a ``HierarchyElement`` only if its metadata has the
            same values for these criteria as the other atoms in the
            ``HierarchyElement``.
        iterator_name
            The name of the iterator that will be exposed to access the hierarchy
            elements generated by this scheme
        """
        if (type(uniqueness_criteria) is not list) and (
            type(uniqueness_criteria) is not tuple
        ):
            raise TypeError(
                f"'uniqueness_criteria' kwarg must be a list or a tuple of strings,"
                f" received {repr(uniqueness_criteria)} "
                f"(type {type(uniqueness_criteria)}) instead."
            )

        for criterion in uniqueness_criteria:
            if type(criterion) is not str:
                raise TypeError(
                    f"Each item in the 'uniqueness_criteria' kwarg must be a string,"
                    f" received {repr(criterion)} "
                    f"(type {type(criterion)}) instead."
                )

        if type(iterator_name) is not str:
            raise TypeError(
                f"'iterator_name' kwarg must be a string, received {repr(iterator_name)} "
                f"(type {type(iterator_name)}) instead."
            )
        self.parent = parent
        self.uniqueness_criteria = uniqueness_criteria
        self.iterator_name = iterator_name

        self.hierarchy_elements: List[HierarchyElement] = list()

    def to_dict(self):
        """
        Serialize this object to a basic dict of strings, ints, and floats
        """
        return_dict = dict()
        return_dict["uniqueness_criteria"] = self.uniqueness_criteria
        return_dict["iterator_name"] = self.iterator_name
        return_dict["hierarchy_elements"] = [
            e.to_dict() for e in self.hierarchy_elements
        ]
        return return_dict

    def perceive_hierarchy(self):
        """
        Prepare the parent ``Molecule`` for iteration according to this scheme.

        Groups the atoms of the parent of this ``HierarchyScheme`` according to
        their metadata, and creates ``HierarchyElement`` objects suitable for
        iteration over the parent. Atoms missing the metadata fields in
        this object's ``uniqueness_criteria`` tuple will have those spots
        populated with the string ``'None'``.

        This method overwrites the scheme's ``hierarchy_elements`` attribute in
        place. Each ``HierarchyElement`` in the scheme's `hierarchy_elements`
        attribute is `static` --- that is, it is updated only when
        `perceive_hierarchy()` is called, and `not` on-the-fly when atom
        metadata is modified.
        """
        from collections import defaultdict

        self.hierarchy_elements = list()
        # Determine which atoms should get added to which HierarchyElements
        hier_eles_to_add = defaultdict(list)
        for atom in self.parent.atoms:
            atom_key = list()
            for field_key in self.uniqueness_criteria:
                if field_key in atom.metadata:
                    atom_key.append(atom.metadata[field_key])
                else:
                    atom_key.append("None")

            hier_eles_to_add[tuple(atom_key)].append(atom)

        # Create the actual HierarchyElements
        for atom_key, atoms_to_add in hier_eles_to_add.items():
            atom_indices = [p.molecule_atom_index for p in atoms_to_add]
            self.add_hierarchy_element(atom_key, atom_indices)

        self.sort_hierarchy_elements()

    def add_hierarchy_element(self, identifier, atom_indices):
        """
        Instantiate a new HierarchyElement belonging to this HierarchyScheme.

        This is the main way to instantiate new HierarchyElements.

        Parameters
        ----------
        identifier : tuple of str and int
            Tuple of metadata values (not keys) that define the uniqueness
            criteria for this element
        atom_indices : iterable int
            The indices of atoms in ``scheme.parent`` that are in this
            element
        """
        new_hier_ele = HierarchyElement(self, identifier, atom_indices)
        self.hierarchy_elements.append(new_hier_ele)
        return new_hier_ele

    def sort_hierarchy_elements(self):
        """
        Semantically sort the HierarchyElements belonging to this object, according to
        their identifiers.
        """
        # hard-code the sort_func value here, since it's hard to serialize safely
        def sort_func(x):
            return version.parse(".".join([str(i) for i in x.identifier]))

        self.hierarchy_elements.sort(key=sort_func)

    def __str__(self):
        return (
            f"HierarchyScheme with uniqueness_criteria '{self.uniqueness_criteria}', iterator_name "
            f"'{self.iterator_name}', and {len(self.hierarchy_elements)} elements"
        )

    def __repr__(self):
        return self.__str__()


class HierarchyElement:
    """An element in a metadata hierarchy scheme, such as a residue or chain."""

    def __init__(self, scheme, identifier, atom_indices):
        """
        Create a new hierarchy element.

        Parameters
        ----------

        scheme : HierarchyScheme
            The scheme to which this ``HierarchyElement`` belongs
        id : tuple of str and int
            Tuple of metadata values (not keys) that define the uniqueness
            criteria for this element
        atom_indices : iterable int
            The indices of particles in ``scheme.parent`` that are in this
            element
        """
        self.scheme = scheme
        self.identifier = identifier
        self.atom_indices = deepcopy(atom_indices)
        for id_component, uniqueness_component in zip(
            identifier, scheme.uniqueness_criteria
        ):
            setattr(self, uniqueness_component, id_component)

    def to_dict(self):
        """
        Serialize this object to a basic dict of strings, ints, and floats
        """
        return_dict = dict()
        return_dict["identifier"] = self.identifier
        return_dict["atom_indices"] = self.atom_indices
        return return_dict

    @property
    def atoms(self):
        for atom_index in self.atom_indices:
            yield self.parent.atoms[atom_index]

    def atom(self, index: int):
        """
        Get atom with a specified index.

        Parameters
        ----------
        index : int

        Returns
        -------
        atom : openff.toolkit.topology.molecule.Atom
        """
        return self.parent.atoms[self.atom_indices[index]]

    @property
    def parent(self) -> FrozenMolecule:
        return self.scheme.parent

    def __str__(self):
        return (
            f"HierarchyElement {self.identifier} of iterator '{self.scheme.iterator_name}' containing "
            f"{len(self.atom_indices)} atom(s)"
        )

    def __repr__(self):
        return self.__str__()
