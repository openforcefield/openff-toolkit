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
import pathlib
import warnings
from collections import UserDict
from copy import deepcopy
from functools import cmp_to_key
from typing import (
    TYPE_CHECKING,
    Any,
    DefaultDict,
    Generator,
    Iterable,
    Literal,
    Optional,
    Sequence,
    TextIO,
    Type,
    TypeAlias,
    TypeVar,
    Union,
    overload,
)

import networkx as nx
import numpy as np
from openff.units import Quantity, unit
from openff.units.elements import MASSES, SYMBOLS
from openff.utilities.exceptions import MissingOptionalDependencyError

from openff.toolkit.utils.constants import DEFAULT_AROMATICITY_MODEL
from openff.toolkit.utils.exceptions import (
    AtomMappingWarning,
    BondExistsError,
    HierarchyIteratorNameConflictError,
    HierarchySchemeNotFoundException,
    HierarchySchemeWithIteratorNameAlreadyRegisteredException,
    IncompatibleShapeError,
    IncompatibleTypeError,
    IncompatibleUnitError,
    InvalidAtomMetadataError,
    InvalidBondOrderError,
    InvalidConformerError,
    MissingConformersError,
    MissingPartialChargesError,
    MoleculeParseError,
    MultipleMoleculesInPDBError,
    RemapIndexError,
    SmilesParsingError,
    UnsupportedFileTypeError,
)
from openff.toolkit.utils.serialization import Serializable
from openff.toolkit.utils.toolkits import (
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
    import IPython.display
    import nglview

    from openff.toolkit.topology._mm_molecule import _SimpleAtom, _SimpleMolecule

# TODO: Can we have the `ALLOWED_*_MODELS` list automatically appear in the docstrings below?
# TODO: Should `ALLOWED_*_MODELS` be objects instead of strings?

# TODO: Allow all OpenEye aromaticity models to be used with OpenEye names?
#       Only support OEAroModel_MDL in RDKit version?

TKR: TypeAlias = Union[ToolkitRegistry, ToolkitWrapper]
FM = TypeVar("FM", bound="FrozenMolecule")
P = TypeVar("P", bound="Particle")
A = TypeVar("A", bound="Atom")
B = TypeVar("B", bound="Bond")


class MoleculeDeprecationWarning(UserWarning):
    """Warning for deprecated portions of the Molecule API."""


class Particle(Serializable):
    """
    Base class for all particles in a molecule.

    A particle object could be an ``Atom`` or similar.

    .. warning :: This API is experimental and subject to change.
    """

    @property
    def molecule(self) -> "FrozenMolecule":
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
    def molecule_particle_index(self) -> int:
        """
        Returns the index of this particle in its molecule
        """
        return self._molecule.atoms.index(self)

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
    def from_dict(cls: type[P], d: dict) -> P:
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

       * Do we want to support the addition of arbitrary additional properties,
         such as floating point quantities (e.g. ``charge``), integral
         quantities (such as ``id`` or ``serial`` index in a PDB file),
         or string labels (such as Lennard-Jones types)?

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

        Parameters
        ----------
        atomic_number : int
            Atomic number of the atom. Must be non-negative and non-zero.
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
        if not isinstance(atomic_number, int):
            raise ValueError(f"atomic number must be int, found {type(atomic_number)}")
        if atomic_number <= 0:
            raise ValueError(f"atomic number must be positive, given {atomic_number}.")

        self._atomic_number = atomic_number

        # Use the setter here, since it will handle either ints or Quantities
        # and it is designed to quickly process ints
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
    def add_bond(self, bond: "Bond"):
        """Adds a bond that this atom is involved in

        .. todo :: Is this how we want to keep records?

        Parameters
        ----------
        bond: an openff.toolkit.topology.molecule.Bond
            A bond involving this atom
        """

        self._bonds.append(bond)

    def to_dict(self) -> dict[str, Union[str, int, bool, dict[Any, Any]]]:
        """Return a dict representation of the atom."""
        # TODO: Should this be implicit in the atom ordering when saved?
        # atom_dict['molecule_atom_index'] = self._molecule_atom_index
        return {
            "atomic_number": self._atomic_number,
            "formal_charge": self._formal_charge.m,  # Trust that the unit is e
            "is_aromatic": self._is_aromatic,
            "stereochemistry": self._stereochemistry,
            "name": self._name,
            "metadata": dict(self._metadata),
        }

    @classmethod
    def from_dict(cls: Type[A], atom_dict: dict) -> A:
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
            self._formal_charge = Quantity(other, unit.elementary_charge)
        elif isinstance(other, Quantity):
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
                    f"Found object of type {type(other)}."
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

    @partial_charge.setter
    def partial_charge(self, charge):
        if self.molecule.partial_charges is None:
            raise MissingPartialChargesError(
                "Cannot set individual atom's partial charge if it is in a molecule with no partial charges. "
                "Instead, use the `Molecule.partial_charges` setter. If this behavior is important to you, "
                "please raise an issue describing your use case."
            )

        if not isinstance(charge, (Quantity, float)):
            raise ValueError(
                "Cannot set partial charge with an object that is not a openff.unit.Quantity or float. "
                f"Found object of type {type(charge)}."
            )

        if isinstance(charge, float):
            charge = Quantity(charge, unit.elementary_charge)

        if not isinstance(charge.m, float):
            raise ValueError(
                "Cannot set partial charge with an object that is not a wrapped int or float. "
                f"Found unit-wrapped {type(charge.m)}."
            )

        molecule_partial_charges = self.molecule.partial_charges
        molecule_partial_charges[self.molecule_atom_index] = charge

        self.molecule.partial_charges = molecule_partial_charges

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
    def stereochemistry(self, value: Literal["CW", "CCW", None]):
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
    def mass(self) -> Quantity:
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
        if type(other) is not str:
            raise ValueError(
                f"In setting atom name. Expected str, received {other} (type {type(other)})."
            )
        self._name = other

    @property
    def bonds(self):
        """
        The list of ``Bond`` objects this atom is involved in.

        """
        return self._bonds

    @property
    def bonded_atoms(self) -> Generator["Atom", None, None]:
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

    def is_in_ring(
        self,
        toolkit_registry: ToolkitRegistry = GLOBAL_TOOLKIT_REGISTRY,
    ) -> bool:
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
    def molecule_atom_index(self) -> int:
        """
        The index of this Atom within the the list of atoms in the parent ``Molecule``.
        """
        if self._molecule is None:
            raise ValueError("This Atom does not belong to a Molecule object")
        if "_molecule_atom_index" in self.__dict__:
            return self._molecule_atom_index  # type: ignore[has-type]
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
        assert type(atom1) is Atom
        assert type(atom2) is Atom
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

    def to_dict(self) -> dict[str, Union[int, bool, str, float]]:
        """
        Return a dict representation of the bond.

        """
        return {
            "atom1": self.atom1.molecule_atom_index,
            "atom2": self.atom2.molecule_atom_index,
            "bond_order": self._bond_order,
            "is_aromatic": self._is_aromatic,
            "stereochemistry": self._stereochemistry,
            "fractional_bond_order": self._fractional_bond_order,
        }

    @classmethod
    def from_dict(cls: Type[B], molecule: FM, d: dict) -> B:  # type: ignore[override]
        """Create a Bond from a dict representation."""
        # TODO: This is not used anywhere (`Molecule._initialize_bonds_from_dict()` just calls grabs
        #       the two atoms and calls `Molecule._add_bond`). Remove or change that?
        # TODO: There is no point in feeding in a `molecule` argument since `Bond.__init__` already
        #       requires (and checks) that the two atoms are part of the same molecule
        d["atom1"] = molecule.atoms[d["atom1"]]
        d["atom2"] = molecule.atoms[d["atom2"]]

        return cls(
            atom1=d["atom1"],
            atom2=d["atom2"],
            bond_order=d["bond_order"],
            is_aromatic=d["is_aromatic"],
            stereochemistry=d["stereochemistry"],
            fractional_bond_order=d["fractional_bond_order"],
        )

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    @property
    def atom1_index(self) -> int:
        return self.molecule.atoms.index(self._atom1)

    @property
    def atom2_index(self) -> int:
        return self.molecule.atoms.index(self._atom2)

    @property
    def atoms(self):
        return (self._atom1, self._atom2)

    @property
    def bond_order(self):
        return self._bond_order

    @bond_order.setter
    def bond_order(self, value):
        if isinstance(value, int):
            self._bond_order = value
        else:
            raise InvalidBondOrderError(
                "Only integer bond orders may be passed to `Bond.bond_order` setter. "
                "For aromatic bonds, instead kekulize the input structure and use "
                "the resulting integer bond orders. If performing partial bond "
                "order-based parameter interpolation, consider using "
                "`Bond.fractional_bond_order`."
            )

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
        # TODO: This is an impossible state (the constructor requires that atom1 and atom2
        #       are in a molecule, the same molecule, and sets that as self._molecule).
        #       Should we remove this?
        assert (
            self._molecule is None
        ), "Bond.molecule is already set and can only be set once"
        self._molecule = value

    @property
    def molecule_bond_index(self) -> int:
        """
        The index of this Bond within the the list of bonds in ``Molecules``.

        """
        if self._molecule is None:
            # TODO: This is unreachable; see `Bond.molecule` setter
            raise ValueError("This Atom does not belong to a Molecule object")
        return self._molecule.bonds.index(self)

    def is_in_ring(
        self,
        toolkit_registry: ToolkitRegistry = GLOBAL_TOOLKIT_REGISTRY,
    ) -> bool:
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
        file_format: Optional[str] = None,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo: bool = False,
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
        >>> molecule_copy = FrozenMolecule(serialized_molecule)

        """

        self._cached_smiles: dict[str, str] = dict()

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

            if isinstance(other, FrozenMolecule) and not loaded:
                self._copy_initializer(other)
                loaded = True
            if isinstance(other, Molecule) and not loaded:
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
            if (
                isinstance(other, (str, pathlib.Path))
                or hasattr(other, "read")
                and not loaded
            ):
                try:
                    mol = Molecule.from_file(
                        other,
                        file_format=file_format,
                        toolkit_registry=toolkit_registry,
                        allow_undefined_stereo=allow_undefined_stereo,
                    )  # returns a list only if multiple molecules are found
                    if type(mol) is list:
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
        """``True`` if the molecule has unique atom names, ``False`` otherwise."""
        return _has_unique_atom_names(self)

    def generate_unique_atom_names(self):
        """
        Generate unique atom names from the element symbol and count.

        Names are generated from the elemental symbol and the number of times
        that element is found in the molecule. The character 'x' is appended to
        these generated names to reduce the odds that they clash with an atom
        name or type imported from another source. For example, generated atom
        names might begin 'C1x', 'H1x', 'O1x', 'C2x', etc.
        """
        return _generate_unique_atom_names(self)

    def _validate(self):
        """
        Validate the molecule, ensuring it has unique atom names

        """
        if not self.has_unique_atom_names:
            self.generate_unique_atom_names()

    def strip_atom_stereochemistry(
        self,
        smarts: str,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
    ):
        """Delete stereochemistry information for certain atoms, if it is present.
        This method can be used to "normalize" molecules imported from different cheminformatics
        toolkits, which differ in which atom centers are considered stereogenic.

        Parameters
        ----------
        smarts: str
            Tagged SMARTS with a single atom with index 1. Any matches for this atom will have any assigned
            stereocheistry information removed.
        toolkit_registry : a :class:`ToolkitRegistry` or :class:`ToolkitWrapper` object, optional,
            default=GLOBAL_TOOLKIT_REGISTRY
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for I/O operations

        """
        matches = self.chemical_environment_matches(
            smarts,
            toolkit_registry=toolkit_registry,
        )

        for match in set(matches):
            atom_idx = match[0]
            self.atoms[atom_idx].stereochemistry = None

    ####################################################################################################
    # Safe serialization
    ####################################################################################################

    def to_dict(self) -> dict:
        """
        Return a dictionary representation of the molecule.

        .. todo ::

           * Document the representation standard.
           * How do we do version control with this standard?

        Returns
        -------
        molecule_dict : dict
            A dictionary representation of the molecule.

        """
        from openff.toolkit.utils.utils import serialize_numpy

        # typing.TypedDict might make this cleaner
        # https://mypy.readthedocs.io/en/latest/typed_dict.html#typeddict
        molecule_dict: dict[
            str,
            Union[
                None,
                str,
                bytes,
                dict[str, Any],
                list[str],
                list[bytes],
                list[HierarchyElement],
            ],
        ] = dict()
        molecule_dict["name"] = self._name

        molecule_dict["atoms"] = [atom.to_dict() for atom in self._atoms]
        molecule_dict["bonds"] = [bond.to_dict() for bond in self._bonds]

        # TODO: This assumes everything in _properties can safely be deepcopied
        molecule_dict["properties"] = deepcopy(self._properties)
        if hasattr(self, "_cached_properties"):
            molecule_dict["cached_properties"] = deepcopy(self._cached_properties)

        if self._conformers is None:
            molecule_dict["conformers"] = None
        else:
            molecule_dict["conformers_unit"] = "angstrom"
            molecule_dict["conformers"] = [
                serialize_numpy(conf.m_as(unit.angstrom))[0]
                for conf in self._conformers
            ]

        if self._partial_charges is None:
            molecule_dict["partial_charges"] = None
            molecule_dict["partial_charge_unit"] = None

        else:
            molecule_dict["partial_charges"], _ = serialize_numpy(
                self._partial_charges.m_as(unit.elementary_charge)
            )
            molecule_dict["partial_charge_unit"] = "elementary_charge"

        molecule_dict["hierarchy_schemes"] = dict()
        for iter_name, hier_scheme in self._hierarchy_schemes.items():
            molecule_dict["hierarchy_schemes"][iter_name] = hier_scheme.to_dict()  # type: ignore[index]

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
    def from_dict(cls: type[FM], molecule_dict: dict) -> FM:
        """
        Create a new Molecule from a dictionary representation

        Parameters
        ----------
        molecule_dict : dict
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
        molecule_dict : dict
            A dictionary representation of the molecule.
        """
        # TODO: Provide useful exception messages if there are any failures

        self._initialize()
        self.name = molecule_dict["name"]
        for atom_dict in molecule_dict["atoms"]:
            self._add_atom(**atom_dict, invalidate_cache=False)

        for bond_dict in molecule_dict["bonds"]:
            bond_dict["atom1"] = int(bond_dict["atom1"])
            bond_dict["atom2"] = int(bond_dict["atom2"])
            self._add_bond(**bond_dict, invalidate_cache=False)

        if molecule_dict["partial_charges"] is None:
            self._partial_charges = None
        else:
            from openff.toolkit.utils.utils import deserialize_numpy

            self._partial_charges = Quantity(
                deserialize_numpy(molecule_dict["partial_charges"], (self.n_atoms,)),
                unit.Unit(molecule_dict["partial_charge_unit"]),
            )

        if molecule_dict["conformers"] is None:
            self._conformers = None
        else:
            from openff.toolkit.utils.utils import deserialize_numpy

            self._conformers = [
                Quantity(
                    deserialize_numpy(ser_conf, (self.n_atoms, 3)),
                    unit.Unit(molecule_dict["conformers_unit"]),
                )
                for ser_conf in molecule_dict["conformers"]
            ]

        self._properties = deepcopy(molecule_dict["properties"])

        for iter_name, hierarchy_scheme_dict in molecule_dict[
            "hierarchy_schemes"
        ].items():
            # It's important that we do NOT call `add_hierarchy_scheme` here, since we
            # need to deserialize these HierarchyElements exactly as they were serialized,
            # even if that conflicts with the current values in atom metadata.
            new_hier_scheme = HierarchyScheme(
                self,
                tuple(hierarchy_scheme_dict["uniqueness_criteria"]),
                iter_name,
            )
            self._hierarchy_schemes[iter_name] = new_hier_scheme

            for element_dict in hierarchy_scheme_dict["hierarchy_elements"]:
                new_hier_scheme.add_hierarchy_element(
                    tuple(element_dict["identifier"]), element_dict["atom_indices"]
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
        self._bonds = list()  # list of bonds between Atom objects
        self._properties = {}  # Attached properties to be preserved
        # self._cached_properties = None # Cached properties (such as partial charges) can be recomputed as needed
        self._partial_charges = None
        self._conformers = None  # Optional conformers
        self._hill_formula = None  # Cached Hill formula
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
        self._initialize_from_dict(other.to_dict())

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

    def add_default_hierarchy_schemes(self, overwrite_existing: bool = True):
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

        See also
        --------
        HierarchyScheme, Molecule.add_hierarchy_scheme,
        Molecule.update_hierarchy_schemes, Molecule.perceive_residues,
        """
        self._add_chain_hierarchy_scheme(overwrite_existing=overwrite_existing)
        self._add_residue_hierarchy_scheme(overwrite_existing=overwrite_existing)

    def _add_chain_hierarchy_scheme(self, overwrite_existing: bool = True):
        """Add ``chain`` hierarchy scheme."""
        if overwrite_existing:
            if "chains" in self._hierarchy_schemes.keys():
                self.delete_hierarchy_scheme("chains")

        self.add_hierarchy_scheme(("chain_id",), "chains")

    def _add_residue_hierarchy_scheme(self, overwrite_existing: bool = True):
        """Add ``residue`` hierarchy scheme."""
        if overwrite_existing:
            if "residues" in self._hierarchy_schemes.keys():
                self.delete_hierarchy_scheme("residues")

        self.add_hierarchy_scheme(
            ("chain_id", "residue_number", "insertion_code", "residue_name"), "residues"
        )

    def add_hierarchy_scheme(
        self,
        uniqueness_criteria: Iterable[str],
        iterator_name: str,
    ) -> "HierarchyScheme":
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
            elements generated by this scheme. Must not match an existing attribute
            of the ``Molecule``, i.e. ``atoms``, ``angles``, etc.

        Returns
        -------
        new_hier_scheme : openff.toolkit.topology.HierarchyScheme
            The newly created HierarchyScheme

        See also
        --------
        Molecule.add_default_hierarchy_schemes, Molecule.hierarchy_schemes,
        Molecule.delete_hierarchy_scheme,  Molecule.update_hierarchy_schemes,
        HierarchyScheme,
        """
        if iterator_name in self._hierarchy_schemes:
            msg = (
                f'Can not add iterator with name "{iterator_name}" to this molecule, as iterator '
                f"name is already used by {self._hierarchy_schemes[iterator_name]}"
            )
            raise HierarchySchemeWithIteratorNameAlreadyRegisteredException(msg)
        elif iterator_name in dir(self):
            raise HierarchyIteratorNameConflictError(
                f"Can not add iterator with name {iterator_name} to this molecule as an "
                "attribute with that name already exists."
            )
        new_hier_scheme = HierarchyScheme(
            self,
            uniqueness_criteria,
            iterator_name,
        )
        self._hierarchy_schemes[iterator_name] = new_hier_scheme
        self.update_hierarchy_schemes([iterator_name])
        return new_hier_scheme

    @property
    def hierarchy_schemes(self) -> dict[str, "HierarchyScheme"]:
        """
        The hierarchy schemes available on the molecule.

        Hierarchy schemes allow iteration over groups of atoms according to
        their metadata. For more information, see
        :class:`~openff.toolkit.topology.molecule.HierarchyScheme`.

        Returns
        -------
        A dict of the form {str: HierarchyScheme}
            The HierarchySchemes associated with the molecule.

        See also
        --------
        Molecule.add_hierarchy_scheme, Molecule.delete_hierarchy_scheme,
        Molecule.update_hierarchy_schemes, Topology.hierarchy_iterator,
        HierarchyScheme
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

        See also
        --------
        Molecule.add_hierarchy_scheme, Molecule.update_hierarchy_schemes,
        Molecule.hierarchy_schemes, HierarchyScheme
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

        See also
        --------
        Molecule.add_hierarchy_scheme, Molecule.delete_hierarchy_schemes,
        Molecule.hierarchy_schemes, HierarchyScheme
        """
        if iter_names is None:
            iter_names = self._hierarchy_schemes.keys()

        for iter_name in iter_names:
            hierarchy_scheme = self._hierarchy_schemes[iter_name]
            hierarchy_scheme.perceive_hierarchy()

    def __getattr__(self, name: str) -> list["HierarchyElement"]:
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
        isomeric: bool = True,
        explicit_hydrogens: bool = True,
        mapped: bool = False,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
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
        # Figure out which toolkit should be used to create the SMILES
        if isinstance(toolkit_registry, ToolkitRegistry):
            to_smiles_method = toolkit_registry.resolve("to_smiles")
        elif isinstance(toolkit_registry, ToolkitWrapper):
            to_smiles_method = toolkit_registry.to_smiles  # type: ignore[attr-defined]
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
        cls: Type[FM],
        inchi: str,
        allow_undefined_stereo: bool = False,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
        name: str = "",
    ) -> FM:
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

        name : str, default=""
            An optional name for the output molecule

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
                name=name,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            molecule = toolkit.from_inchi(  # type: ignore[attr-defined]
                inchi,
                _cls=cls,
                allow_undefined_stereo=allow_undefined_stereo,
                name=name,
            )
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to from_inchi. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        return molecule

    def to_inchi(
        self,
        fixed_hydrogens: bool = False,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
    ) -> str:
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
            inchi = toolkit.to_inchi(self, fixed_hydrogens=fixed_hydrogens)  # type: ignore[attr-defined]
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to to_inchi. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        return inchi

    def to_inchikey(
        self,
        fixed_hydrogens: bool = False,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
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
            inchi_key = toolkit.to_inchikey(self, fixed_hydrogens=fixed_hydrogens)  # type: ignore[attr-defined]
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to to_inchikey. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        return inchi_key

    @classmethod
    def from_smiles(
        cls: Type[FM],
        smiles: str,
        hydrogens_are_explicit: bool = False,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo: bool = False,
        name: str = "",
    ) -> FM:
        """
        Construct a ``Molecule`` from a SMILES representation

        The order of atoms in the ``Molecule`` is unspecified and may change
        from version to version or with different toolkits. SMILES atom
        indices (also known as atom maps) are not used to order atoms; instead,
        they are stored in the produced molecule's properties attribute,
        accessible via ``molecule.properties["atom_map"]``. The atom map is
        stored as a dictionary mapping molecule atom indices to SMILES atom
        maps. To order atoms according to SMILES atom indices, see
        :py:meth:`Molecule.from_mapped_smiles`, which helpfully raises an
        exception if any atom map is missing, duplicated, or out-of-range,
        or else :py:meth:`Molecule.remap` for arbitrary remaps.

        Parameters
        ----------
        smiles
            The SMILES representation of the molecule.
        hydrogens_are_explicit
            If ``True``, forbid the cheminformatics toolkit from inferring
            hydrogen atoms not explicitly specified in the SMILES.
        toolkit_registry
            The cheminformatics toolkit to use to interpret the SMILES.
        allow_undefined_stereo
            Whether to accept SMILES with undefined stereochemistry. If
            ``False``, an exception will be raised if a SMILES with undefined
            stereochemistry is passed into this function.
        name : str, default=""
            An optional name for the output molecule


        Raises
        ------
        RadicalsNotSupportedError
            If any atoms in the input molecule contain radical electrons.

        Examples
        --------

        Create a ``Molecule`` representing toluene from SMILES:

        >>> molecule = Molecule.from_smiles('Cc1ccccc1')

        Create a ``Molecule`` representing phenol from SMILES with the oxygen
        at atom index 0 (SMILES indices begin at 1):

        >>> molecule = Molecule.from_smiles('c1ccccc1[OH:1]')
        >>> molecule = molecule.remap(
        ...     {k: v - 1 for k, v in molecule.properties["atom_map"].items()},
        ...     partial=True,
        ... )
        >>> assert molecule.atom(0).symbol == "O"

        See Also
        --------
        from_mapped_smiles, remap

        """
        if isinstance(toolkit_registry, ToolkitRegistry):
            molecule = toolkit_registry.call(
                "from_smiles",
                smiles,
                hydrogens_are_explicit=hydrogens_are_explicit,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
                name=name,
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            toolkit = toolkit_registry
            molecule = toolkit.from_smiles(  # type: ignore[attr-defined]
                smiles,
                hydrogens_are_explicit=hydrogens_are_explicit,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
                name=name,
            )
        else:
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to from_smiles. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        if "atom_map" in molecule._properties:
            if len(molecule._properties["atom_map"]) == molecule.n_atoms:
                warnings.warn(
                    "Warning! Fully mapped SMILES pattern passed to `from_smiles`. The atom map is "
                    "stored as a property in `Molecule._properties`, but these indices are NOT "
                    "used to determine atom ordering. To use these indices for atom ordering, use "
                    "`Molecule.from_mapped_smiles`.",
                    AtomMappingWarning,
                    stacklevel=2,
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
        mol1: Union["FrozenMolecule", "_SimpleMolecule", nx.Graph],
        mol2: Union["FrozenMolecule", "_SimpleMolecule", nx.Graph],
        return_atom_map: bool = False,
        aromatic_matching: bool = True,
        formal_charge_matching: bool = True,
        bond_order_matching: bool = True,
        atom_stereochemistry_matching: bool = True,
        bond_stereochemistry_matching: bool = True,
        strip_pyrimidal_n_atom_stereo: bool = True,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
    ) -> tuple[bool, Optional[dict[int, int]]]:
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
            Return a ``dict`` containing the atomic mapping, otherwise ``None``.
            Only processed if inputs are isomorphic, will always return ``None`` if
            inputs are not isomorphic.

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
            [dict[int,int]] ordered by mol1 indexing {mol1_index: mol2_index}
            If molecules are not isomorphic given input arguments, will return None instead of dict.
        """
        import networkx as nx

        _cls = FrozenMolecule

        if isinstance(mol1, nx.Graph) and isinstance(mol2, nx.Graph):
            pass

        elif isinstance(mol1, nx.Graph):
            assert isinstance(mol2, _cls)

        elif isinstance(mol2, nx.Graph):
            assert isinstance(mol1, _cls)

        else:
            # static methods (by definition) know nothing about their class,
            # so the class to compare to must be hard-coded here
            if not (isinstance(mol1, _cls) and isinstance(mol2, _cls)):
                return False, None

        def _object_to_n_atoms(obj):
            if isinstance(obj, FrozenMolecule):
                return obj.n_atoms
            elif isinstance(obj, nx.Graph):
                return obj.number_of_nodes()
            else:
                raise TypeError(
                    "are_isomorphic accepts a NetworkX Graph or OpenFF "
                    + f"(Frozen)Molecule, not {type(obj)}"
                )

        # Quick number of atoms check. Important for large molecules
        if _object_to_n_atoms(mol1) != _object_to_n_atoms(mol2):
            return False, None

        # If the number of atoms match, check the Hill formula
        if Molecule._object_to_hill_formula(mol1) != Molecule._object_to_hill_formula(
            mol2
        ):
            return False, None

        # Do a quick check to see whether the inputs are totally identical (including being in the same atom order)
        if isinstance(mol1, FrozenMolecule) and isinstance(mol2, FrozenMolecule):
            if mol1._is_exactly_the_same_as(mol2):
                if return_atom_map:
                    return True, {i: i for i in range(mol1.n_atoms)}
                else:
                    return True, None

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
            edge_match_func = None  # type: ignore

        # Here we should work out what data type we have, also deal with lists?
        def to_networkx(data: Union[FrozenMolecule, nx.Graph]) -> nx.Graph:
            """For the given data type, return the networkx graph"""
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

    def is_isomorphic_with(
        self, other: Union["FrozenMolecule", "_SimpleMolecule", nx.Graph], **kwargs
    ) -> bool:
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
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
        n_conformers: int = 10,
        rms_cutoff: Optional[Quantity] = None,
        clear_existing: bool = True,
        make_carboxylic_acids_cis: bool = True,
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
        rms_cutoff : openff.unit.Quantity-wrapped float, in units of distance, optional, default=None
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
            return toolkit.generate_conformers(  # type: ignore[attr-defined]
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

    def _make_carboxylic_acids_cis(
        self, toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY
    ):
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
        self._conformers = [Quantity(conf, unit.angstrom) for conf in conformers]

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

    def assign_partial_charges(
        self,
        partial_charge_method: str,
        strict_n_conformers: bool = False,
        use_conformers: Optional[Iterable[Quantity]] = None,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
        normalize_partial_charges: bool = True,
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

        By default, the conformers on the input molecule are not used
        in the charge calculation. Instead, any conformers needed for
        the charge calculation are generated by this method. If this
        behavior is undesired, specific conformers can be provided via the
        ``use_conformers`` argument.

        ELF10 methods will neither fail nor warn when fewer than the
        expected number of conformers could be generated, as many small
        molecules are too rigid to provide a large number of conformers. Note
        that only the ``"am1bccelf10"`` partial charge method uses ELF
        conformer selection; the ``"am1bcc"`` method only uses a single
        conformer. This may confuse users as the `ToolkitAM1BCC`_ SMIRNOFF
        tag in a force field file defines that AM1BCC-ELF10 should be used
        if the OpenEye Toolkits are available.

        For more supported charge methods and their details, see the
        corresponding methods in each toolkit wrapper:

        - :meth:`OpenEyeToolkitWrapper.assign_partial_charges \
          <openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.assign_partial_charges>`
        - :meth:`RDKitToolkitWrapper.assign_partial_charges \
          <openff.toolkit.utils.toolkits.RDKitToolkitWrapper.assign_partial_charges>`
        - :meth:`AmberToolsToolkitWrapper.assign_partial_charges \
          <openff.toolkit.utils.toolkits.AmberToolsToolkitWrapper.assign_partial_charges>`
        - :meth:`BuiltInToolkitWrapper.assign_partial_charges \
          <openff.toolkit.utils.toolkits.BuiltInToolkitWrapper.assign_partial_charges>`

        .. _ToolkitAM1BCC: https://openforcefield.github.io/standards/standards/smirnoff/\
            #toolkitam1bcc-temporary-support-for-toolkit-based-am1-bcc-partial-charges

        Parameters
        ----------
        partial_charge_method : string
            The partial charge calculation method to use for partial charge
            calculation.
        strict_n_conformers : bool, default=False
            Whether to raise an exception if an invalid number of conformers is
            provided for the given charge method. If this is False and an
            invalid number of conformers is found, a warning will be raised.
        use_conformers : Arrays with shape (n_atoms, 3) and dimensions of distance
            Coordinates to use for partial charge calculation. If ``None``, an
            appropriate number of conformers will be generated.
        toolkit_registry
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for the
            calculation.
        normalize_partial_charges : bool, default=True
            Whether to offset partial charges so that they sum to the total
            formal charge of the molecule. This is used to prevent accumulation
            of rounding errors when the partial charge assignment method returns
            values at limited precision.

        Examples
        --------

        Generate AM1 Mulliken partial charges. Conformers for the AM1
        calculation are generated automatically:

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.assign_partial_charges('am1-mulliken')

        To use pre-generated conformations, use the ``use_conformers`` argument:

        >>> molecule = Molecule.from_smiles('CCCCCC')
        >>> molecule.generate_conformers(n_conformers=1)
        >>> molecule.assign_partial_charges(
        ...     'am1-mulliken',
        ...     use_conformers=molecule.conformers
        ... )

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
        WARN_LARGE_MOLECULES: set[str] = {
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
                    stacklevel=2,
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
        bond_order_model: Optional[str] = None,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
        use_conformers: Optional[Iterable[Quantity]] = None,
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
        use_conformers : iterable of openff.unit.Quantity(np.array) with shape (n_atoms, 3) and dimension of distance,
            optional, default=None
            The conformers to use for fractional bond order calculation. If ``None``, an appropriate number
            of conformers will be generated by an available ``ToolkitWrapper``.

        Examples
        --------

        >>> from openff.toolkit import Molecule
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
            return toolkit.assign_fractional_bond_orders(  # type: ignore[attr-defined]
                self, bond_order_model=bond_order_model, use_conformers=use_conformers
            )
        else:
            raise InvalidToolkitRegistryError(
                f"Invalid toolkit_registry passed to assign_fractional_bond_orders. "
                f"Expected ToolkitRegistry or ToolkitWrapper. Got {type(toolkit_registry)}."
            )

    def _invalidate_cached_properties(self) -> None:
        """
        Indicate that the chemical entity has been altered.
        """
        self._conformers = None
        self._partial_charges = None
        self._propers: set[tuple[Atom, Atom, Atom, Atom]] = set()
        self._impropers: set[tuple[Atom, Atom, Atom, Atom]] = set()

        self._hill_formula = None
        self._cached_smiles = dict()
        # TODO: Clear fractional bond orders
        self._ordered_connection_table_hash = None
        for atom in self.atoms:
            if "_molecule_atom_index" in atom.__dict__:
                del atom.__dict__["_molecule_atom_index"]

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
        self,
        ignore_functional_groups: Optional[list[str]] = None,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
    ) -> list[Bond]:
        """
        Find all bonds classed as rotatable ignoring any matched to the ``ignore_functional_groups`` list.

        Parameters
        ----------
        ignore_functional_groups: optional, list[str], default=None,
            A list of bond SMARTS patterns to be ignored when finding rotatable bonds.

        toolkit_registry: openff.toolkit.utils.toolkits.ToolkitRegistry
            or openff.toolkit.utils.toolkits.ToolkitWrapperl, optional, default=None
            :class:`ToolkitRegistry` or :class:`ToolkitWrapper` to use for SMARTS matching

        Returns
        -------
        bonds: list[openff.toolkit.topology.molecule.Bond]
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
                    raise ValueError(
                        "Argument ignore_functional_groups must be iterable or str. "
                        f"Found type {type(ignore_functional_groups)=}"
                    )

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
        atomic_number: int,
        formal_charge: int,
        is_aromatic: bool,
        stereochemistry: Optional[str] = None,
        name: Optional[str] = None,
        metadata=None,
        invalidate_cache: bool = True,
    ) -> int:
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
            raise ValueError(
                "Invalid inputs to molecule._add_bond. Expected ints or Atoms. "
                f"Received {atom1} (type {type(atom1)}) and {atom2} (type {type(atom2)}) "
            )
        # TODO: Check to make sure bond does not already exist
        if atom1_atom.is_bonded_to(atom2_atom):
            raise BondExistsError(
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
        coordinates: openff.unit.Quantity(np.array) with shape (n_atoms, 3) and dimension of distance
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

        if isinstance(coordinates, Quantity):
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

        tmp_conf = Quantity(
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
        partial_charges : a openff.unit.Quantity - wrapped numpy array [1 x n_atoms] or None
            The partial charges on the molecule's atoms. Returns None if no charges have been specified.
        """
        return self._partial_charges

    @partial_charges.setter
    def partial_charges(self, charges):
        """
        Set the atomic partial charges for this molecule.

        Parameters
        ----------
        charges : None or a openff.unit.Quantity - wrapped numpy array [1 x n_atoms]
            The partial charges to assign to the molecule. If not None, must be in units compatible with
            openff.unit.elementary_charge

        """
        if charges is None:
            self._partial_charges = None
            return

        if not hasattr(charges, "shape"):
            raise IncompatibleTypeError(
                "Unsupported type passed to partial_charges setter. "
                f"Found object of type {type(charges)}. "
                "Expected openff.units.unit.Quantity"
            )

        if not charges.shape == (self.n_atoms,):
            raise IncompatibleShapeError(
                "Unsupported shape passed to partial_charges setter. "
                f"Found shape {charges.shape}, expected {(self.n_atoms,)}"
            )

        if isinstance(charges, Quantity):
            if charges.units in unit.elementary_charge.compatible_units():
                self._partial_charges = charges.astype(float)
            else:
                raise IncompatibleUnitError(
                    "Unsupported unit passed to partial_charges setter. "
                    f"Found unit {charges.units}, expected {unit.elementary_charge}"
                )

        elif hasattr(charges, "unit"):
            from openmm import unit as openmm_unit

            if not isinstance(charges, openmm_unit.Quantity):
                raise IncompatibleUnitError(
                    "Unsupported type passed to partial_charges setter. "
                    f"Found object of type {type(charges)}."
                )

            else:
                from openff.units.openmm import from_openmm

                converted = from_openmm(charges)
                if converted.units in unit.elementary_charge.compatible_units():
                    self._partial_charges = converted.astype(float)
                else:
                    raise IncompatibleUnitError(
                        "Unsupported unit passed to partial_charges setter. "
                        f"Found unit {converted.units}, expected {unit.elementary_charge}"
                    )

        else:
            raise IncompatibleTypeError(
                "Unsupported type passed to partial_charges setter. "
                f"Found object of type {type(charges)}, "
                "expected openff.units.unit.Quantity"
            )

    @property
    def n_atoms(self) -> int:
        """
        The number of Atom objects.
        """
        return len(self._atoms)

    @property
    def n_bonds(self) -> int:
        """
        The number of Bond objects in the molecule.
        """
        return len(self._bonds)

    @property
    def n_angles(self) -> int:
        """Number of angles in the molecule."""
        self._construct_angles()
        return len(self._angles)

    @property
    def n_propers(self) -> int:
        """Number of proper torsions in the molecule."""
        self._construct_torsions()
        return len(self._propers)

    @property
    def n_impropers(self) -> int:
        """Number of possible improper torsions in the molecule."""
        self._construct_torsions()
        return len(self._impropers)

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
    def bonds(self) -> list[Bond]:
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
    def angles(self) -> set[tuple[Atom, Atom, Atom]]:
        """
        Get an iterator over all i-j-k angles.
        """
        self._construct_angles()
        return self._angles

    @property
    def torsions(self) -> set[tuple[Atom, Atom, Atom, Atom]]:
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
    def propers(self) -> set[tuple[Atom, Atom, Atom, Atom]]:
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
    def impropers(self) -> set[tuple[Atom, Atom, Atom, Atom]]:
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

        return self._impropers

    @property
    def smirnoff_impropers(self) -> set[tuple[Atom, Atom, Atom, Atom]]:
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
        return {
            improper
            for improper in self.impropers
            if len(self._bonded_atoms[improper[1]]) == 3
        }

    @property
    def amber_impropers(self) -> set[tuple[Atom, Atom, Atom, Atom]]:
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
        self._construct_torsions()

        return {
            (improper[1], improper[0], improper[2], improper[3])
            for improper in self.smirnoff_impropers
        }

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
            tuples (len 2) of atom that are separated by ``n`` bonds.

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
            raise ValueError("Molecule name must be a string")

    @property
    def properties(self) -> dict[str, Any]:
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
        if self._hill_formula is None:
            atom_nums = [atom.atomic_number for atom in self.atoms]
            self._hill_formula = _atom_nums_to_hill_formula(atom_nums)

        return self._hill_formula

    @staticmethod
    def _object_to_hill_formula(obj: Union["FrozenMolecule", "nx.Graph"]) -> str:
        """Take a Molecule or NetworkX graph and generate its Hill formula.
        This provides a backdoor to the old functionality of Molecule.to_hill_formula, which
        was a static method that duck-typed inputs of Molecule or graph objects."""
        import networkx as nx

        if isinstance(obj, FrozenMolecule):
            return obj.to_hill_formula()
        elif isinstance(obj, nx.Graph):
            return _networkx_graph_to_hill_formula(obj)
        else:
            raise TypeError(
                "_object_to_hill_formula accepts a NetworkX Graph or OpenFF "
                + f"(Frozen)Molecule, not {type(obj)}"
            )

    def chemical_environment_matches(
        self,
        query: str,
        unique: bool = False,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
    ):
        """Find matches in the molecule for a SMARTS string

        Parameters
        ----------
        query : str
            SMARTS string (with one or more tagged atoms).
        unique : bool, default=False
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
        # Use specified cheminformatics toolkit to determine matches with specified aromaticity model
        # TODO: Simplify this by requiring a toolkit registry for the molecule?
        # TODO: Do we have to pass along an aromaticity model?
        if isinstance(toolkit_registry, ToolkitRegistry):
            matches = toolkit_registry.call(
                "find_smarts_matches",
                self,
                query,
                unique=unique,
                raise_exception_types=[],
            )
        elif isinstance(toolkit_registry, ToolkitWrapper):
            matches = toolkit_registry.find_smarts_matches(  # type: ignore[attr-defined]
                self,
                query,
                unique=unique,
            )
        else:
            raise InvalidToolkitRegistryError(
                "'toolkit_registry' must be either a ToolkitRegistry or a ToolkitWrapper"
            )

        return matches

    @classmethod
    def from_iupac(
        cls: Type[FM],
        iupac_name: str,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo: bool = False,
        **kwargs,
    ) -> FM:
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
            molecule = toolkit.from_iupac(  # type: ignore[attr-defined]
                iupac_name,
                allow_undefined_stereo=allow_undefined_stereo,
                _cls=cls,
                **kwargs,
            )
        else:
            raise InvalidToolkitRegistryError(
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
            raise InvalidToolkitRegistryError(
                "Invalid toolkit_registry passed to to_iupac. Expected ToolkitRegistry or ToolkitWrapper. "
                f"Got {type(toolkit_registry)}"
            )

        # TODO: Can `to_iupac` fail if given a well-behaved OFFMol/OEMol?
        result = to_iupac_method(self)
        return result

    @classmethod
    def from_topology(cls: Type[FM], topology) -> FM:
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

        >>> from openff.toolkit import Molecule, Topology
        >>> topology = Topology.from_molecules(Molecule.from_smiles('[CH4]'))
        >>> molecule = Molecule.from_topology(topology)

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

        >>> from openff.toolkit import Molecule
        >>> molecule = Molecule.from_iupac('imatinib')
        >>> topology = molecule.to_topology()

        """
        from openff.toolkit.topology import Topology

        return Topology.from_molecules(self)

    @classmethod
    def from_file(
        cls: Type[FM],
        file_path: Union[str, pathlib.Path, TextIO],
        file_format=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo: bool = False,
    ) -> Union[FM, list[FM]]:
        """
        Create one or more molecules from a file

        .. todo::

           * Extend this to also include some form of .offmol Open Force Field Molecule format?
           * Generalize this to also include file-like objects?

        Parameters
        ----------
        file_path : str, pathlib.Path, or file-like object
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

        >>> from openff.toolkit import Molecule
        >>> from openff.toolkit.utils.utils import get_data_file_path
        >>> sdf_file_path = get_data_file_path("molecules/toluene.sdf")
        >>> molecule = Molecule.from_file(sdf_file_path)

        """

        if file_format is None:
            if isinstance(file_path, pathlib.Path):
                file_path: str = file_path.as_posix()  # type: ignore[no-redef]
            if not isinstance(file_path, str):
                raise ValueError(
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

        if isinstance(file_path, (str, pathlib.Path)):
            if isinstance(file_path, pathlib.Path):
                file_path = file_path.as_posix()
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
            raise MoleculeParseError(f"Unable to read molecule from file: {file_path}")
        elif len(mols) == 1:
            return mols[0]

        return mols

    @classmethod
    @requires_package("openmm")
    def from_polymer_pdb(
        cls: Type[FM],
        file_path: Union[str, pathlib.Path, TextIO],
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        name: str = "",
    ) -> FM:
        """
        Loads a polymer from a PDB file.

        Also see :py:meth:`Topology.from_multicomponent_pdb`, which can do
        everything this method can and more.

        Currently only supports proteins with canonical amino acids that are
        either uncapped or capped by ACE/NME groups, but may later be extended
        to handle other common polymers, or accept user-defined polymer
        templates. Only one polymer chain may be present in the PDB file, and it
        must be the only molecule present.

        Connectivity and bond orders are assigned by matching SMARTS codes for
        the supported residues against atom names. The PDB file must include
        all atoms with the correct standard atom names described in the
        `PDB Chemical Component Dictionary <https://www.wwpdb.org/data/ccd>`_.
        Residue names are used to assist trouble-shooting failed assignments,
        but are not used in the actual assignment process.

        Metadata such as residues, chains, and atom names are recorded in the
        ``Atom.metadata`` attribute, which is a dictionary mapping from
        strings like "residue_name" to the appropriate value. ``from_polymer_pdb``
        returns a molecule that can be iterated over with the ``.residues`` and
        ``.chains`` attributes, as well as the usual ``.atoms``.

        Parameters
        ----------
        file_path : str or file object
            PDB information to be passed to OpenMM PDBFile object for loading
        toolkit_registry = ToolkitWrapper or ToolkitRegistry. Default = None
            Either a ToolkitRegistry, ToolkitWrapper
        name : str, default=""
            An optional name for the output molecule

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule

        Raises
        ------

        UnassignedChemistryInPDBError
            If an atom or bond could not be assigned; the exception will
            provide a detailed diagnostic of what went wrong.

        MultipleMoleculesInPDBError
            If all atoms and bonds could be assigned, but the PDB includes
            multiple chains or molecules.

        """
        import io

        import openmm.unit as openmm_unit
        from openmm.app import PDBFile

        warnings.warn(
            "`Molecule.from_polymer_pdb` is deprecated in favor of `Topology.from_pdb`, the recommended "
            "method for loading PDB files. This method will be removed in a future release of the OpenFF Toolkit.",
            MoleculeDeprecationWarning,
            stacklevel=2,
        )

        if isinstance(toolkit_registry, ToolkitWrapper):
            toolkit_registry = ToolkitRegistry([toolkit_registry])

        if isinstance(file_path, (str, io.TextIOWrapper)):
            pass
        elif isinstance(file_path, pathlib.Path):
            file_path = file_path.as_posix()
        else:
            raise ValueError(f"Unexpected type {type(file_path)}")

        pdb = PDBFile(file_path)

        # Kludgy fix for the fact that RDKitToolkitWrapper uses new substructure spec.
        # Hopefully this will be short-lived as we can deprecate this method entirely in favor of
        # Topology.from_pdb, which only uses the RDKit backend.
        resolved_method = toolkit_registry.resolve("_polymer_openmm_topology_to_offmol")
        if "RDKit" in str(resolved_method):
            substructure_file_path = get_data_file_path(
                "proteins/aa_residues_substructures_explicit_bond_orders_with_caps_explicit_connectivity.json"
            )
        else:
            substructure_file_path = get_data_file_path(
                "proteins/aa_residues_substructures_explicit_bond_orders_with_caps.json"
            )

        with open(substructure_file_path, "r") as subfile:
            substructure_dictionary = json.load(subfile)

        offmol = toolkit_registry.call(
            "_polymer_openmm_topology_to_offmol",
            cls,
            pdb.topology,
            substructure_dictionary,
        )

        coords = Quantity(
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
            offmol.atoms[i].metadata["insertion_code"] = atom.residue.insertionCode
            offmol.atoms[i].metadata["chain_id"] = atom.residue.chain.id
        offmol.add_default_hierarchy_schemes()

        if offmol._has_multiple_molecules():
            raise MultipleMoleculesInPDBError(
                "This PDB has multiple molecules. The OpenFF Toolkit requires "
                + "that only one molecule is present in a PDB. Try splitting "
                + "each molecule into its own PDB with another tool, and "
                + "load any small molecules with Molecule.from_pdb_and_smiles."
            )

        offmol.name = name

        return offmol

    def _has_multiple_molecules(self) -> bool:
        import networkx as nx

        graph = self.to_networkx()
        num_disconnected_subgraphs = sum(1 for _ in nx.connected_components(graph))
        return num_disconnected_subgraphs > 1

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
                Quantity(np.zeros((self.n_atoms, 3), dtype=float), unit.angstrom)
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
        if type(file_path) is str:
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
        molecules: list[openff.toolkit.topology.Molecule]
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
        undefined_only: bool = False,
        max_isomers: int = 20,
        rationalise: bool = True,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
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
        molecules: list[openff.toolkit.topology.Molecule]
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
            molecules = toolkit_registry.enumerate_stereoisomers(  # type: ignore[attr-defined]
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
        molecules: list[openff.toolkit.topology.Molecule],
            A list of the protomers of the input molecules not including the input.
        """

        toolkit = OpenEyeToolkitWrapper()
        molecules = toolkit.enumerate_protomers(molecule=self, max_states=max_states)

        return molecules

    @classmethod
    @RDKitToolkitWrapper.requires_toolkit()
    def from_rdkit(
        cls: Type[FM],
        rdmol,
        allow_undefined_stereo: bool = False,
        hydrogens_are_explicit: bool = False,
    ) -> FM:
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

        >>> from openff.toolkit import Molecule
        >>> from rdkit import Chem
        >>> rdmol = Chem.MolFromSmiles("CCO")
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
        aromaticity_model : str, optional, default="OEAroModel_MDL"
            The aromaticity model to use. Only OEAroModel_MDL is supported.

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
    def from_openeye(
        cls: Type[FM],
        oemol,
        allow_undefined_stereo: bool = False,
    ) -> "Molecule":
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

        >>> from openff.toolkit import Molecule
        >>> from openeye import oechem
        >>> oemol = oechem.OEMol()
        >>> oechem.OESmilesToMol(oemol, '[H]C([H])([H])C([H])([H])O[H]')
        True
        >>> molecule = Molecule.from_openeye(oemol)

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
        cls: Type[FM],
        mapped_smiles: str,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo: bool = False,
    ) -> FM:
        """
        Create a ``Molecule`` from a SMILES string, ordering atoms from mappings

        SMILES strings support mapping integer indices to each atom by ending a
        bracketed atom declaration with a colon followed by a 1-indexed
        integer:

        .. code:
            "[H:3][C:1](=[O:2])[H:4]"

        This method creates a ``Molecule`` from such a SMILES string whose atoms
        are ordered according to the mapping. Each atom must be mapped exactly
        once; any duplicate, missing, or out-of-range mappings will cause the
        method to fail.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        mapped_smiles: str
            A mapped SMILES string with explicit hydrogens.
        toolkit_registry
            Cheminformatics toolkit to use for SMILES-to-molecule conversion
        allow_undefined_stereo
            If false, raise an exception if the SMILES contains undefined
            stereochemistry.

        Returns
        ----------
        offmol
            An OpenFF molecule instance.

        Raises
        --------
        SmilesParsingError
            If the given SMILES had no indexing picked up by the toolkits, or if
            the indexing is missing indices.
        RemapIndexError
            If the mapping has duplicate or out-of-range indices.

        Examples
        --------

        Create a mapped chlorofluoroiodomethane molecule and check the atoms
        are placed accordingly:

        >>> molecule = Molecule.from_mapped_smiles(
        ...     "[Cl:2][C@:1]([F:3])([I:4])[H:5]"
        ... )
        >>> assert molecule.atom(0).symbol == "C"
        >>> assert molecule.atom(1).symbol == "Cl"
        >>> assert molecule.atom(2).symbol == "F"
        >>> assert molecule.atom(3).symbol == "I"
        >>> assert molecule.atom(4).symbol == "H"

        See Also
        --------
        from_smiles, remap
        """

        # create the molecule from the smiles and check we have the right number of indexes
        # in the mapped SMILES
        warnings.filterwarnings("ignore", category=AtomMappingWarning)

        offmol = cls.from_smiles(
            mapped_smiles,
            hydrogens_are_explicit=True,
            toolkit_registry=toolkit_registry,
            allow_undefined_stereo=allow_undefined_stereo,
        )

        # https://stackoverflow.com/a/53763710
        # this might be better: https://docs.python.org/3/library/warnings.html#warnings.catch_warnings
        warnings.filterwarnings("default", category=AtomMappingWarning)

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
        # the order is mapping = dict[current_index: new_index]
        # first renumber the mapping dict indexed from 0, currently from 1 as 0 indicates no mapping in toolkits
        adjusted_mapping = dict((current, new - 1) for current, new in mapping.items())

        return offmol.remap(adjusted_mapping, current_to_new=True)

    @classmethod
    @requires_package("qcelemental")
    def from_qcschema(
        cls: type[FM],
        qca_record,
        client=None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        allow_undefined_stereo: bool = False,
    ) -> FM:
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
        >>> offmol = Molecule.from_qcschema(
        ...     client.query_molecules(molecular_formula="C7H12N2O4")[0]
        ... )

        Get Molecule from a QCArchive optimization entry:

        >>> from qcportal import FractalClient
        >>> client = FractalClient()
        >>> optds = client.get_collection(
        ...     "OptimizationDataset",
        ...     "SMIRNOFF Coverage Set 1"
        ... )
        >>> offmol = Molecule.from_qcschema(optds.get_entry('coc(o)oc-0'))

        Same as above, but with conformer(s) from initial molecule(s) by
        providing client to database:

        >>> offmol = Molecule.from_qcschema(
        ...     optds.get_entry('coc(o)oc-0'),
        ...     client=client
        ... )

        Raises
        -------
        AttributeError
            If the record dict can not be made from ``qca_record``, or if the
            provided ``client`` could not retrieve the initial molecule.
        KeyError
            If the record does not contain the
            ``canonical_isomeric_explicit_hydrogen_mapped_smiles``.
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

            geometry = Quantity(
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
    def from_pdb_and_smiles(
        cls: type[FM],
        file_path,
        smiles,
        allow_undefined_stereo: bool = False,
        name: str = "",
    ) -> FM:
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
        name : str, default=""
            An optional name for the output molecule

        Returns
        --------
        molecule : openff.toolkit.Molecule
            An OFFMol instance with ordering the same as used in the PDB file.

        Raises
        ------
        InvalidConformerError
            If the SMILES and PDB molecules are not isomorphic.
        """
        warnings.warn(
            "`Molecule.from_pdb_and_smiles` is deprecated in favor of `Topology.from_pdb`, the recommended "
            "method for loading PDB files. This method will be removed in a future release of the OpenFF Toolkit.",
            MoleculeDeprecationWarning,
            stacklevel=2,
        )

        toolkit = RDKitToolkitWrapper()
        return toolkit.from_pdb_and_smiles(
            file_path, smiles, allow_undefined_stereo, _cls=cls, name=name
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

    def remap(
        self,
        mapping_dict: dict[int, int],
        current_to_new: bool = True,
        partial: bool = False,
    ):
        """
        Reorder the atoms in the molecule according to the given mapping dict.

        The mapping dict must be a dictionary mapping atom indices to atom
        indices. Each atom index must be an integer in the half-open interval
        ``[0, n_atoms)``; ie, it must be a valid index into the ``self.atoms``
        list. All atom indices in the molecule must be mapped from and to
        exactly once unless ``partial=True`` is given, in which case they must
        be mapped no more than once. Missing (unless ``partial=True``),
        out-of-range (including non-integer), or duplicate indices are not
        allowed in the ``mapping_dict`` and will lead to an exception.

        By default, the mapping dict's keys are the source indices and its
        values are destination indices, but this can be changed with the
        ``current_to_new`` argument.

        The keys of the ``self.properties["atom_map"]`` property are updated for
        the new ordering. Other values of the properties dictionary are
        transferred unchanged.

        .. warning :: This API is experimental and subject to change.

        Parameters
        ----------
        mapping_dict
            A dictionary of the mapping between indices. The mapping should be
            indexed starting from 0 for both the source and destination; note
            that SMILES atom mapping is typically 1-based.
        current_to_new
            If this is ``True``, then ``mapping_dict`` is of the form
            ``{current_index: new_index}``; otherwise, it is of the form
            ``{new_index: current_index}``.
        partial
            If ``False`` (the default), an exception will be raised if any atom
            is lacking a destination in the atom map. Note that if this is
            ``True``, atoms without entries in the mapping dict may be moved in
            addition to those in the dictionary. Note that partial maps must
            still be in-range and not include duplicates.

        Returns
        -------
        new_molecule :  openff.toolkit.topology.molecule.Molecule
            A copy of the molecule in the new order.

        Raises
        ------
        RemapIndexError
            When an out-of-range, duplicate, or missing index is found in the
            ``mapping_dict``.

        See Also
        --------
        from_mapped_smiles
        """

        # make sure the size of the mapping matches the current molecule
        if len(mapping_dict) > self.n_atoms or (
            len(mapping_dict) < self.n_atoms and not partial
        ):
            raise RemapIndexError(
                f"The number of mapping indices ({len(mapping_dict)}) does not "
                + f"match the number of atoms in this molecule ({self.n_atoms})"
            )

        # make two mapping dicts we need new to old for atoms
        # and old to new for bonds
        if current_to_new:
            cur_to_new = mapping_dict
            new_to_cur = dict(zip(mapping_dict.values(), mapping_dict.keys()))
        else:
            new_to_cur = mapping_dict
            cur_to_new = dict(zip(mapping_dict.values(), mapping_dict.keys()))

        # Make sure that there were no duplicate indices
        if len(new_to_cur) != len(cur_to_new):
            raise RemapIndexError(
                "There must be no duplicate source or destination indices in"
                + " mapping_dict"
            )

        if any(
            not (isinstance(i, int) and 0 <= i < self.n_atoms)
            for i in [*new_to_cur] + [*cur_to_new]
        ):
            raise RemapIndexError(
                f"All indices in a mapping_dict for a molecule with {self.n_atoms}"
                + f" atoms must be integers between 0 and {self.n_atoms-1}"
            )

        # If a partial map is allowed, complete it
        if partial and len(mapping_dict) < self.n_atoms:
            # Get a set of all the unspecified destination indices
            available_indices = {i for i in range(self.n_atoms) if i not in new_to_cur}
            # Find the atoms that can be left unmoved and don't move them
            for i in range(self.n_atoms):
                if i not in cur_to_new and i not in new_to_cur:
                    available_indices.remove(i)
                    cur_to_new[i] = i
                    new_to_cur[i] = i
            # Fill in the remaining indices
            for i in range(self.n_atoms):
                if i not in cur_to_new:
                    j = available_indices.pop()
                    cur_to_new[i] = j
                    new_to_cur[j] = i

        new_molecule = self.__class__()
        new_molecule.name = self.name

        try:
            # add the atoms list
            for i in range(self.n_atoms):
                # get the old atom info
                old_atom = self._atoms[new_to_cur[i]]
                new_molecule._add_atom(**old_atom.to_dict())
        # this is the first time we access the mapping; catch an index error
        # here corresponding to mapping that starts from 0 or higher
        except (KeyError, IndexError):
            raise RemapIndexError(
                f"The mapping supplied is missing a destination index for atom {i}"
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

        # remap the conformers, there can be more than one
        if self.conformers is not None:
            for conformer in self.conformers:
                new_conformer = np.zeros((self.n_atoms, 3))
                for i in range(self.n_atoms):
                    new_conformer[i] = conformer[new_to_cur[i]].m_as(unit.angstrom)
                new_molecule._add_conformer(new_conformer * unit.angstrom)

        # move any properties across
        new_molecule._properties = deepcopy(self._properties)

        # remap the atom map
        if "atom_map" in new_molecule.properties and isinstance(
            new_molecule.properties["atom_map"], dict
        ):
            new_molecule.properties["atom_map"] = {
                cur_to_new.get(k, k): v
                for k, v in new_molecule.properties["atom_map"].items()
            }

        return new_molecule

    def to_openeye(
        self,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
        aromaticity_model: str = DEFAULT_AROMATICITY_MODEL,
    ):
        """
        Create an OpenEye molecule

        Requires the OpenEye toolkit to be installed.

        .. todo ::

           * Use stored conformer positions instead of an argument.
           * Should the aromaticity model be specified in some other way?

        Parameters
        ----------
        aromaticity_model : str, optional, default="OEAroModel_MDL"
            The aromaticity model to use. Only OEAroModel_MDL is supported.

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

        if isinstance(toolkit_registry, ToolkitWrapper):
            return toolkit_registry.to_openeye(  # type: ignore[attr-defined]
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
        if not hasattr(self, "_angles"):
            self._construct_bonded_atoms_list()
            self._angles = set()
            for atom1 in self._atoms:
                for atom2 in self._bonded_atoms[atom1]:
                    for atom3 in self._bonded_atoms[atom2]:
                        if atom1 == atom3:
                            continue
                        if atom1.molecule_atom_index < atom3.molecule_atom_index:
                            self._angles.add((atom1, atom2, atom3))
                        else:
                            self._angles.add((atom3, atom2, atom1))

    def _construct_torsions(self) -> None:
        """
        Construct sets containing the atoms improper and proper torsions

        Impropers are constructed with the central atom listed second
        """
        if not hasattr(self, "_torsions"):
            self._construct_bonded_atoms_list()

            self._propers = set()
            self._impropers = set()

            for atom1 in self._atoms:
                for atom2 in self._bonded_atoms[atom1]:
                    for atom3 in self._bonded_atoms[atom2]:
                        if atom1 == atom3:
                            continue
                        for atom4 in self._bonded_atoms[atom3]:
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

                        for atom3i in self._bonded_atoms[atom2]:
                            if atom3i == atom3:
                                continue
                            if atom3i == atom1:
                                continue

                            improper = (atom1, atom2, atom3, atom3i)
                            self._impropers.add(improper)

            self._torsions = self._propers | self._impropers

    def _construct_bonded_atoms_list(self) -> None:
        """
        Construct list of all atoms each atom is bonded to.

        """
        # TODO: Add this to cached_properties
        if not hasattr(self, "_bonded_atoms"):
            self._bonded_atoms: dict[Atom, set[Atom]] = dict()
            for atom in self._atoms:
                self._bonded_atoms[atom] = set()
            for bond in self._bonds:
                atom1 = self.atoms[bond.atom1_index]
                atom2 = self.atoms[bond.atom2_index]
                self._bonded_atoms[atom1].add(atom2)
                self._bonded_atoms[atom2].add(atom1)

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
        return atom2 in self._bonded_atoms[atom1]

    def get_bond_between(self, i: Union[int, Atom], j: Union[int, Atom]) -> Bond:
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
        See FrozenMolecule.__init__

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
        atomic_number: int,
        formal_charge: int,
        is_aromatic: bool,
        stereochemistry: Optional[str] = None,
        name: Optional[str] = None,
        metadata: Optional[dict[str, Union[int, str]]] = None,
    ) -> int:
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
        >>> bond_idx = molecule.add_bond(C, H1, 1, False)
        >>> bond_idx = molecule.add_bond(C, H2, 1, False)
        >>> bond_idx = molecule.add_bond(C, H3, 1, False)
        >>> bond_idx = molecule.add_bond(C, H4, 1, False)
        >>> molecule.to_smiles(explicit_hydrogens=False)
        'C'

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
        atom1: Union[int, "Atom"],
        atom2: Union[int, "Atom"],
        bond_order: int,
        is_aromatic: bool,
        stereochemistry: Optional[str] = None,
        fractional_bond_order: Optional[float] = None,
    ) -> int:
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

        Examples
        --------
        For an example of use, see :py:meth:`add_atom`.
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

    def add_conformer(self, coordinates: Quantity) -> int:
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

    @overload
    def visualize(
        self,
        backend: Literal["rdkit"],
    ) -> "IPython.display.SVG":
        ...

    @overload
    def visualize(
        self,
        backend: Literal["openeye"],
    ) -> "IPython.display.Image":
        ...

    @overload
    def visualize(
        self,
        backend: Literal["nglview"],
    ) -> "nglview.NGLWidget":
        ...

    def visualize(
        self,
        backend: str = "rdkit",
        width: int = 500,
        height: int = 300,
        show_all_hydrogens: bool = True,
    ) -> Union["IPython.display.SVG", "IPython.display.Image", "nglview.NGLWidget"]:
        """
        Render a visualization of the molecule in Jupyter

        Parameters
        ----------
        backend : str, optional, default='rdkit'
            The visualization engine to use. Choose from:

            - ``"rdkit"``
            - ``"openeye"``
            - ``"nglview"`` (requires conformers)

        width : int, default=500
            Width of the generated representation (only applicable to
            ``backend="openeye"`` or ``backend="rdkit"``)
        height : int, default=300
            Width of the generated representation (only applicable to
            ``backend="openeye"`` or ``backend="rdkit"``)
        show_all_hydrogens : bool, default=True
            Whether to explicitly depict all hydrogen atoms. (only applicable to
            ``backend="openeye"`` or ``backend="rdkit"``)

        Returns
        -------
        object
            Depending on the backend chosen:

            - rdkit → IPython.display.SVG
            - openeye → IPython.display.Image
            - nglview → nglview.NGLWidget

        """
        import inspect

        from openff.toolkit.utils.toolkits import OPENEYE_AVAILABLE, RDKIT_AVAILABLE

        backend = backend.lower()

        if backend == "nglview":
            try:
                import nglview as nv
            except ImportError:
                raise MissingOptionalDependencyError("nglview")

            signature = inspect.signature(Molecule.visualize).parameters
            if (width != signature["width"].default) or (
                height != signature["height"].default
            ):
                warnings.warn(
                    f"Arguments `width` and `height` are ignored with {backend=}."
                    f"Found non-default values {width=} and {height=}",
                    stacklevel=2,
                )

            if self.conformers is None:
                raise MissingConformersError(
                    "Visualizing with NGLview requires that the molecule has "
                    f"conformers, found {self.conformers=}"
                )

            else:
                from openff.toolkit.utils._viz import MoleculeNGLViewTrajectory

                try:
                    widget = nv.NGLWidget(
                        MoleculeNGLViewTrajectory(
                            molecule=self,
                            ext="MOL2",
                        )
                    )
                except ValueError:
                    widget = nv.NGLWidget(
                        MoleculeNGLViewTrajectory(
                            molecule=self,
                            ext="PDB",
                        )
                    )

                widget.clear_representations()
                widget.add_representation(
                    "licorice",
                    sele="*" if show_all_hydrogens else "NOT hydrogen",
                    radius=0.25,
                    multipleBond=True,
                )

                return widget

        if backend == "rdkit":
            if RDKIT_AVAILABLE:
                from IPython.display import SVG
                from rdkit.Chem.Draw import (  # type: ignore[import-untyped]
                    rdDepictor,
                    rdMolDraw2D,
                )
                from rdkit.Chem.rdmolops import RemoveHs  # type: ignore[import-untyped]

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
                    "trying to use OpenEye for visualization.",
                    stacklevel=2,
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

    def perceive_residues(
        self,
        substructure_file_path: Optional[str] = None,
        strict_chirality: bool = True,
    ):
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
                for smarts in inner_dict.keys():
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
                self.atoms[atom_idx].metadata["residue_number"] = str(residue_num + 1)
                self.atoms[atom_idx].metadata["insertion_code"] = " "
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
        raise ValueError("The graph must be a NetworkX graph.")

    atom_nums = list(dict(graph.nodes(data="atomic_number", default=1)).values())
    return _atom_nums_to_hill_formula(atom_nums)


def _atom_nums_to_hill_formula(atom_nums: list[int]) -> str:
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
    Union[tuple[Atom, Atom], tuple["_SimpleAtom", "_SimpleAtom"]], None, None
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
        tuples (len 2) of atom that are separated by ``n`` bonds.
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

    See also
    --------
    Molecule.add_default_hierarchy_schemes, Molecule.add_hierarchy_scheme,
    Molecule.hierarchy_schemes, Molecule.delete_hierarchy_scheme,
    Molecule.update_hierarchy_schemes, Molecule.perceive_residues,
    Topology.hierarchy_iterator, HierarchyElement
    """

    def __init__(
        self,
        parent: FrozenMolecule,
        uniqueness_criteria: Iterable[str],
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

        self.hierarchy_elements: list[HierarchyElement] = list()

    def to_dict(self) -> dict:
        """
        Serialize this object to a basic dict of strings, ints, and floats
        """
        return_dict: dict[str, str | Sequence[str | int | dict]] = dict()
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

    def add_hierarchy_element(
        self,
        identifier: tuple[Union[str, int]],
        atom_indices: Sequence[int],
    ) -> "HierarchyElement":
        """
        Instantiate a new HierarchyElement belonging to this HierarchyScheme.

        This is the main way to instantiate new HierarchyElements.

        Parameters
        ----------
        identifier : tuple of str and int
            tuple of metadata values (not keys) that define the uniqueness
            criteria for this element
        atom_indices : sequence of int
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

        def compare_hier_identifiers(a, b):
            """A comparison function which can compare hierarchy elements.
            Expects identifiers to be tuples of string and int.
            Attempts to cast strings to int. Assumes that ints are "greater than" strings.

            Returns -1 if a < b, 0 if a==b, and 1 if a>b.

            See https://docs.python.org/3/howto/sorting.html#comparison-functions
            """

            # Iterate over identifier components for comparison
            for val1, val2 in zip(a.identifier, b.identifier):
                # Try converting any strings to ints
                try:
                    val1 = int(val1)
                except ValueError:
                    pass
                try:
                    val2 = int(val2)
                except ValueError:
                    pass

                # If val1 and val2 are the same type, use built-in comparison
                if type(val1) is type(val2):
                    if val1 < val2:
                        return -1
                    elif val1 > val2:
                        return 1
                    else:
                        continue

                # Otherwise, assume that ints are "greater than" strings.
                else:
                    if type(val1) is int:
                        return 1
                    elif type(val2) is int:
                        return -1
            # If we've finished comparing the values in the identifiers without
            # finding one to be greater than the other, then these two identifiers
            # must be equal.
            return 0

        self.hierarchy_elements.sort(key=cmp_to_key(compare_hier_identifiers))

    def __str__(self):
        return (
            f"HierarchyScheme with uniqueness_criteria '{self.uniqueness_criteria}', iterator_name "
            f"'{self.iterator_name}', and {len(self.hierarchy_elements)} elements"
        )

    def __repr__(self):
        return self.__str__()


class HierarchyElement:
    """An element in a metadata hierarchy scheme, such as a residue or chain."""

    def __init__(
        self,
        scheme: HierarchyScheme,
        identifier: tuple[Union[str, int]],
        atom_indices: Sequence[int],
    ):
        """
        Create a new hierarchy element.

        Parameters
        ----------

        scheme : HierarchyScheme
            The scheme to which this ``HierarchyElement`` belongs
        identifier : tuple of str and int
            tuple of metadata values (not keys) that define the uniqueness
            criteria for this element
        atom_indices : sequence of int
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

    def to_dict(self) -> dict[str, Union[tuple[Union[str, int]], Sequence[int]]]:
        """
        Serialize this object to a basic dict of strings and lists of ints.
        """
        return {
            "identifier": self.identifier,
            "atom_indices": self.atom_indices,
        }

    @property
    def n_atoms(self) -> int:
        """
        The number of atoms in this hierarchy element.
        """
        return len(self.atom_indices)

    @property
    def atoms(self) -> Generator["Atom", None, None]:
        """
        Iterator over the atoms in this hierarchy element.
        """
        for atom_index in self.atom_indices:
            yield self.parent.atoms[atom_index]

    def atom(self, index: int) -> Atom:
        """
        Get the atom with the specified index.
        """
        return self.parent.atoms[self.atom_indices[index]]

    @property
    def parent(self) -> FrozenMolecule:
        """
        The parent molecule for this hierarchy element
        """
        return self.scheme.parent

    def __str__(self):
        return (
            f"HierarchyElement {self.identifier} of iterator '{self.scheme.iterator_name}' containing "
            f"{len(self.atom_indices)} atom(s)"
        )

    def __repr__(self):
        return self.__str__()

    @property
    def has_unique_atom_names(self) -> bool:
        """``True`` if the element has unique atom names, ``False`` otherwise."""
        return _has_unique_atom_names(self)

    def generate_unique_atom_names(self):
        """
        Generate unique atom names from the element symbol and count.

        Names are generated from the elemental symbol and the number of times
        that element is found in the hierarchy element. The character 'x' is
        appended to these generated names to reduce the odds that they clash
        with an atom name or type imported from another source. For example,
        generated atom names might begin 'C1x', 'H1x', 'O1x', 'C2x', etc.
        """
        return _generate_unique_atom_names(self)


def _has_unique_atom_names(obj: Union[FrozenMolecule, HierarchyElement]) -> bool:
    """``True`` if the object has unique atom names, ``False`` otherwise."""
    unique_atom_names = set([atom.name for atom in obj.atoms])
    if len(unique_atom_names) < obj.n_atoms:
        return False
    return True


def _generate_unique_atom_names(obj: Union[FrozenMolecule, HierarchyElement]):
    """
    Generate unique atom names from the element symbol and count.

    Names are generated from the elemental symbol and the number of times that
    element is found in the hierarchy element or molecule. The character 'x' is
    appended to these generated names to reduce the odds that they clash with
    an atom name or type imported from another source. For example, generated
    atom names might begin 'C1x', 'H1x', 'O1x', 'C2x', etc.
    """
    from collections import defaultdict

    element_counts: DefaultDict[str, int] = defaultdict(int)
    for atom in obj.atoms:
        symbol = atom.symbol
        element_counts[symbol] += 1
        # TODO: It may be worth exposing this as a user option, i.e. to avoid multiple ligands
        # parameterized with OpenFF clashing because they have atom names like O1x, H3x, etc.
        # i.e. an optional argument could enable a user to `generate_unique_atom_names(blah="y")
        # to have one ligand be O1y, etc.
        # https://github.com/openforcefield/openff-toolkit/pull/1096#pullrequestreview-767227391
        atom.name = symbol + str(element_counts[symbol]) + "x"
