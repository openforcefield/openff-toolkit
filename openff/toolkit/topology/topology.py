"""
Class definitions to represent a molecular system and its chemical components

.. todo::

   * Create MoleculeImage, ParticleImage, AtomImage here. (Or ``MoleculeInstance``?)
   * Create ``MoleculeGraph`` to represent fozen set of atom elements and bonds that can used as a key for compression
   * Add hierarchical way of traversing Topology (chains, residues)
   * Make all classes hashable and serializable.
   * JSON/BSON representations of objects?
   * Use `attrs <http://www.attrs.org/>`_ for object setter boilerplate?

"""

import re
from collections import defaultdict
from collections.abc import MutableMapping
from contextlib import nullcontext
from copy import deepcopy
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Generator,
    Iterable,
    Iterator,
    Literal,
    Optional,
    TextIO,
    Union,
)

import numpy as np
from networkx import Graph
from numpy.typing import NDArray
from openff.units import ensure_quantity
from typing_extensions import TypeAlias

from openff.toolkit import Quantity, unit
from openff.toolkit.topology import Molecule
from openff.toolkit.topology._mm_molecule import (
    _SimpleAtom,
    _SimpleBond,
    _SimpleMolecule,
)
from openff.toolkit.topology.molecule import FrozenMolecule, HierarchyElement
from openff.toolkit.utils import quantity_to_string, string_to_quantity
from openff.toolkit.utils.constants import (
    ALLOWED_AROMATICITY_MODELS,
    DEFAULT_AROMATICITY_MODEL,
)
from openff.toolkit.utils.exceptions import (
    AtomNotInTopologyError,
    BondNotInTopologyError,
    ConstraintExsistsError,
    DuplicateUniqueMoleculeError,
    IncompatibleUnitError,
    InvalidAromaticityModelError,
    InvalidBoxVectorsError,
    InvalidPeriodicityError,
    MissingConformersError,
    MissingUniqueMoleculesError,
    MoleculeNotInTopologyError,
    NotBondedError,
    VirtualSitesUnsupportedError,
    WrongShapeError,
)
from openff.toolkit.utils.serialization import Serializable
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY
from openff.toolkit.utils.utils import get_data_file_path, requires_package

if TYPE_CHECKING:
    import mdtraj
    import openmm.app
    from nglview import NGLWidget
    from openmm.unit import Quantity as OMMQuantity

    from openff.toolkit.topology.molecule import Atom, Bond
    from openff.toolkit.utils import ToolkitRegistry, ToolkitWrapper

TKR: TypeAlias = Union["ToolkitRegistry", "ToolkitWrapper"]
MoleculeLike: TypeAlias = Union["Molecule", "FrozenMolecule", "_SimpleMolecule"]


class _TransformedDict(MutableMapping):
    """A dictionary that transform and sort keys.

    The function __keytransform__ can be inherited to apply an arbitrary
    key-altering function before accessing the keys.

    The function __sortfunc__ can be inherited to specify a particular
    order over which to iterate over the dictionary.

    """

    def __init__(self, *args, **kwargs):
        self.store = dict()
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

    @classmethod
    def _return_possible_index_of(cls, key, possible, permutations):
        """
        Returns canonical ordering of ``key``, given a dictionary of ordered
        ``permutations`` and a list of allowed orders ``possible``.

        Parameters
        ----------
        key: tuple of int
            A key of indices
        possible: list of tuples of int
            List of possible keys
        permutations: dict[tuple[int, ...], int]
            Dictionary of canonical orders
        """
        key = tuple(key)
        possible = [tuple(p) for p in possible]
        impossible = [p for p in possible if p not in permutations]
        if impossible:
            raise ValueError(f"Impossible permutations {impossible} for key {key}!")
        possible_permutations = [k for k in permutations if k in possible]
        for i, permutation in enumerate(possible_permutations):
            if key == permutation:
                return i
        raise ValueError(f"key {key} not in possible {possible}")


# TODO: Encapsulate this atom ordering logic directly into Atom/Bond/Angle/Torsion classes?
class ValenceDict(_TransformedDict):
    """Enforce uniqueness in atom indices."""

    @staticmethod
    def key_transform(key):
        """Reverse tuple if first element is larger than last element."""
        key = tuple(key)
        if key[0] > key[-1]:
            key = tuple(reversed(key))
        return key

    @classmethod
    def index_of(cls, key, possible=None):
        """
        Generates a canonical ordering of the equivalent permutations of ``key`` (equivalent rearrangements of indices)
        and identifies which of those possible orderings this particular ordering is. This method is useful when
        multiple SMARTS patterns might match the same atoms, but local molecular symmetry or the use of
        wildcards in the SMARTS could make the matches occur in arbitrary order.

        This method can be restricted to a subset of the canonical orderings, by providing
        the optional ``possible`` keyword argument. If provided, the index returned by this method will be
        the index of the element in ``possible`` after undergoing the same canonical sorting as above.

        Parameters
        ----------
        key : iterable of int
            A valid key for ValenceDict
        possible : iterable of iterable of int, optional. default=``None``
            A subset of the possible orderings that this match might take.

        Returns
        -------
        index : int
        """
        refkey = cls.key_transform(key)
        permutations = {refkey: 0, refkey[::-1]: 1}
        if possible is not None:
            return cls._return_possible_index_of(
                key, possible=possible, permutations=permutations
            )
        else:
            return permutations[tuple(key)]

    def __keytransform__(self, key):
        return self.key_transform(key)


class SortedDict(_TransformedDict):
    """Enforce uniqueness of atom index tuples, without any restrictions on atom reordering."""

    def __keytransform__(self, key):
        """Sort tuple from lowest to highest."""
        key = tuple(sorted(key))
        return key


class UnsortedDict(_TransformedDict):
    pass


class TagSortedDict(_TransformedDict):
    """
    A dictionary where keys, consisting of tuples of atom indices, are kept unsorted, but only allows one permutation
    of a key to exist. Certain situations require that atom indices are not transformed in any way, such as when the
    tagged order of a match is needed downstream. For example a parameter using charge increments needs the ordering of
    the tagged match, and so transforming the atom indices in any way will cause that information to be lost.

    Because deduplication is needed, we still must avoid the expected situation
    where we must not allow two permutations of the same atoms to coexist. For example,
    a parameter may have matched the indices (2, 0, 1), however a parameter with higher
    priority also matches the same indices, but the tagged order is (1, 0, 2). We need
    to make sure both keys don't exist, or else two parameters will apply to
    the same atoms. We allow the ability to query using either permutation and get
    identical behavior. The defining feature here, then, is that the stored indices are
    in tagged order, but one can supply any permutation and it will resolve to the
    stored value/parameter.

    As a subtle behavior, one must be careful if an external key is used that was not
    supplied from the TagSortedDict object itself. For example:

        >>> x = TagSortedDict({(2, 5, 0): 100})
        >>> y = x[(5, 0, 2)]

    The variable y will be 100, but this does not mean the tagged order is (5, 0, 2)
    as it was supplied from the external tuple. One should either use keys only from
    __iter__ (e.g.  `for k in x`) or one must transform the key:

        >>> key = (5, 0, 2)
        >>> y = x[key]
        >>> key = x.key_transform(key)

    Where `key` will now be `(2, 5, 0)`, as it is the key stored. One can overwrite
    this key with the new one as expected:

        >>> key = (5, 0, 2)
        >>> x[key] = 50

    Now the key `(2, 5, 0)` will no longer exist, as it was replaced with `(5, 0, 2)`.
    """

    def __init__(self, *args, **kwargs):
        # Because keytransform is O(n) due to needing to check the sorted keys,
        # we cache the sorted keys separately to make keytransform O(1) at
        # the expense of storage. This is also better in the long run if the
        # key is long and repeatedly sorting isn't a negligible cost.

        # Set this before calling super init, since super will call the get/set
        # methods implemented here as it populates self via args/kwargs,
        # which will automatically populate _sorted
        self._sorted = SortedDict()

        super().__init__(*args, **kwargs)

    def __setitem__(self, key, value):
        """
        Set the key to value, but only allow one permutation of key to exist. The
        key argument will replace the old permutation:value if it exists.
        """
        key = tuple(key)
        tr_key = self.__keytransform__(key)
        if key != tr_key:
            # this means our new key is a permutation of an existing, so we should
            # replace it
            del self.store[tr_key]
        self.store[key] = value
        # save the sorted version for faster keytransform
        self._sorted[key] = key

    def __keytransform__(self, key):
        """Give the key permutation that is currently stored"""
        # we check if there is a permutation clash by saving the sorted version of
        # each key. If the sorted version of the key exists, then the return value
        # corresponds to the explicit permutation we are storing in self (the public
        # facing key). This permutation may or may not be the same as the key argument
        # supplied. If the key is not present, then no transformation should be done
        # and we should return the key as is.

        # As stated in __init__, the alternative is to, on each call, sort the saved
        # permutations and check if it is equal to the sorted supplied key. In this
        # sense, self._sorted is a cache/lookup table.
        key = tuple(key)
        return self._sorted.get(key, key)

    def key_transform(self, key):
        key = tuple(key)
        return self.__keytransform__(key)

    def clear(self):
        """
        Clear the contents
        """
        self.store.clear()
        self._sorted.clear()


class ImproperDict(_TransformedDict):
    """Symmetrize improper torsions."""

    @staticmethod
    def key_transform(key):
        """
        Reorder tuple in numerical order except for element[1] which is the central atom; it retains its position.
        """
        # Ensure key is a tuple
        key = tuple(key)
        assert len(key) == 4, "Improper keys must be 4 atoms"
        # Retrieve connected atoms
        connectedatoms = [key[0], key[2], key[3]]
        # Sort connected atoms
        connectedatoms.sort()
        # Re-store connected atoms
        key = tuple([connectedatoms[0], key[1], connectedatoms[1], connectedatoms[2]])
        return key

    @classmethod
    def index_of(cls, key, possible=None):
        """
        Generates a canonical ordering of the equivalent permutations of ``key`` (equivalent rearrangements of indices)
        and identifies which of those possible orderings this particular ordering is. This method is useful when
        multiple SMARTS patterns might match the same atoms, but local molecular symmetry or the use of wildcards in
        the SMARTS could make the matches occur in arbitrary order.

        This method can be restricted to a subset of the canonical orderings, by providing the optional ``possible``
        keyword argument. If provided, the index returned by this method will be the index of the element in
        ``possible`` after undergoing the same canonical sorting as above.

        Parameters
        ----------
        key : iterable of int
            A valid key for ValenceDict
        possible : iterable of iterable of int, optional. default=``None``
            A subset of the possible orderings that this match might take.

        Returns
        -------
        index : int
        """
        key = tuple(key)
        assert len(key) == 4
        refkey = cls.key_transform(key)
        permutations = {
            (refkey[0], refkey[1], refkey[2], refkey[3]): 0,
            (refkey[0], refkey[1], refkey[3], refkey[2]): 1,
            (refkey[2], refkey[1], refkey[0], refkey[3]): 2,
            (refkey[2], refkey[1], refkey[3], refkey[0]): 3,
            (refkey[3], refkey[1], refkey[0], refkey[2]): 4,
            (refkey[3], refkey[1], refkey[2], refkey[0]): 5,
        }
        if possible is not None:
            return cls._return_possible_index_of(
                key, possible=possible, permutations=permutations
            )
        else:
            return permutations[key]

    def __keytransform__(self, key):
        return __class__.key_transform(key)


class Topology(Serializable):
    """
    A Topology is a chemical representation of a system containing one or more molecules appearing in a specified
    order.

    .. warning :: This API is experimental and subject to change.

    Examples
    --------

    Import some utilities

    >>> from openmm import app
    >>> from openff.toolkit._tests.utils import get_data_file_path, get_packmol_pdb_file_path
    >>> pdb_filepath = get_packmol_pdb_file_path('cyclohexane_ethanol_0.4_0.6')
    >>> monomer_names = ('cyclohexane', 'ethanol')

    Create a Topology object from a PDB file and sdf files defining the molecular contents

    >>> from openff.toolkit import Molecule, Topology
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
        # Assign cheminformatics models
        self._aromaticity_model = DEFAULT_AROMATICITY_MODEL

        # Initialize storage
        self._initialize()

        # TODO: Try to construct Topology copy from `other` if specified
        if isinstance(other, Topology):
            self.copy_initializer(other)
        elif isinstance(other, FrozenMolecule):
            self.from_molecules([other])
        elif isinstance(other, dict):
            self._initialize_from_dict(other)

    def _initialize(self):
        """
        Initializes a blank Topology.
        """
        self._constrained_atom_pairs = dict()
        self._box_vectors = None
        self._molecules = list()
        self._cached_chemically_identical_molecules = None

    def __iadd__(self, other):
        """Add two Topology objects in-place.

        This method has the following effects:
        * The cache (_cached_chemically_identical_molecules) is reset
        * The constrained atom pairs (Topology.constrained_atom_pairs) is reset
        * Box vectors are **not** updated.

        """
        if self.aromaticity_model != other.aromaticity_model:
            raise InvalidAromaticityModelError(
                "Mismatch in aromaticity models. Trying to add a Topology with aromaticity model "
                f"{other.aromaticity_model} to a Topology with aromaticity model "
                f"{self.aromaticity_model}"
            )

        self._cached_chemically_identical_molecules = None

        constrained_atom_pairs_to_add = other.constrained_atom_pairs
        atom_index_offset = self.n_atoms

        for molecule in other.molecules:
            self._add_molecule_keep_cache(molecule)
        self._invalidate_cached_properties()

        for key, value in constrained_atom_pairs_to_add.items():
            new_key = tuple(index + atom_index_offset for index in key)
            self._constrained_atom_pairs[new_key] = value

        return self

    def __add__(self, other):
        """Add two Topology objects. See Topology.__iadd__ for details."""
        combined = deepcopy(self)
        combined += other

        return combined

    @property
    def unique_molecules(self) -> Iterator[Molecule]:
        """
        Get a list of chemically unique molecules in this Topology.

        Molecules are considered unique if they fail an isomorphism check with
        default values (see :meth:`Molecule.is_isomorphic_with`).
        The order of molecules returned by this property is arbitrary.
        """
        for mol_idx in self.identical_molecule_groups.keys():
            yield deepcopy(self.molecule(mol_idx))

    @property
    def n_unique_molecules(self) -> int:
        """Returns the number of unique molecules in this Topology"""
        return len(self.identical_molecule_groups)

    @classmethod
    def from_molecules(
        cls,
        molecules: Union[MoleculeLike, list[MoleculeLike]],
    ) -> "Topology":
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
        if isinstance(molecules, (FrozenMolecule, _SimpleMolecule)):
            molecules = [molecules]

        # Create Topology and populate it with specified molecules
        topology = cls()
        for molecule in molecules:
            topology._add_molecule_keep_cache(molecule)
        topology._invalidate_cached_properties()

        return topology

    def assert_bonded(self, atom1: Union[int, "Atom"], atom2: Union[int, "Atom"]):
        """
        Raise an exception if the specified atoms are not bonded in the topology.

        Parameters
        ----------
        atom1, atom2 : openff.toolkit.topology.Atom or int
            The atoms or atom topology indices to check to ensure they are bonded

        """
        if isinstance(atom1, int) and isinstance(atom2, int):
            atom1 = self.atom(atom1)
            atom2 = self.atom(atom2)

        if not (self.is_bonded(atom1, atom2)):
            # TODO: Raise more specific exception.
            raise NotBondedError(
                f"Atoms {atom1} and {atom2} are not bonded in topology"
            )

    @property
    def aromaticity_model(self) -> str:
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
        if aromaticity_model not in ALLOWED_AROMATICITY_MODELS:
            raise InvalidAromaticityModelError(
                f"Read aromaticity model {aromaticity_model} which is not in the set of allowed aromaticity models:  "
                f"{ALLOWED_AROMATICITY_MODELS}."
            )

        self._aromaticity_model = aromaticity_model

    @property
    def box_vectors(self):
        """Return the box vectors of the topology, if specified

        Returns
        -------
        box_vectors : unit-wrapped numpy array of shape (3, 3)
            The unit-wrapped box vectors of this topology
        """
        return self._box_vectors

    @box_vectors.setter
    def box_vectors(self, box_vectors):
        """
        Sets the box vectors to be used for this topology.

        Parameters
        ----------
        box_vectors : unit-wrapped numpy array of shape (3, 3)
            The unit-wrapped box vectors

        """
        if box_vectors is None:
            self._box_vectors = None
            return
        if not hasattr(box_vectors, "units"):
            if hasattr(box_vectors, "unit"):
                # this is probably an openmm.unit.Quantity; we should gracefully import OpenMM but
                # the chances of this being an object with the two previous conditions met is low
                box_vectors = ensure_quantity(box_vectors, "openff")

            else:
                raise InvalidBoxVectorsError("Given unitless box vectors")

        # Unit.compatible_units() returns False with itself, for some reason
        if (box_vectors.units != unit.nm) and (
            box_vectors.units not in unit.nm.compatible_units()
        ):
            raise InvalidBoxVectorsError(
                f"Cannot set box vectors with quantities with unit {box_vectors.units}"
            )

        if hasattr(box_vectors, "shape"):
            if box_vectors.shape == (3,):
                # Cannot multiply in-place without ufunc support in Pint
                box_vectors = box_vectors * np.eye(3)
            if box_vectors.shape != (3, 3):
                raise InvalidBoxVectorsError(
                    f"Box vectors must be shape (3, 3). Found shape {box_vectors.shape}"
                )
        else:
            raise InvalidBoxVectorsError(
                f"Cannot set box vectors with object of type {type(box_vectors)}"
            )

        self._box_vectors = box_vectors

    @property
    def is_periodic(self):
        """Return whether or not this Topology is intended to be described with periodic
        boundary conditions."""
        return self.box_vectors is not None

    @is_periodic.setter
    def is_periodic(self, is_periodic):
        """
        Set the partial charge model used for all molecules in the topology.

        Parameters
        ----------
        is_periodic : bool
            Whether or not this Topology is periodici

        """
        if is_periodic is True and self.box_vectors is None:
            raise InvalidPeriodicityError(
                "Cannot set is_periodic to True without box vectors. Set box "
                "vectors directly instead."
            )
        if is_periodic is False and self.box_vectors is not None:
            raise InvalidPeriodicityError(
                "Cannot set is_periodic to False while box vectors are stored. "
                "First set box_vectors to None."
            )

    @property
    def constrained_atom_pairs(self) -> dict[tuple[int], Union[Quantity, bool]]:
        """Returns the constrained atom pairs of the Topology

        Returns
        -------
        constrained_atom_pairs : dict of tuple[int]: Union[unit.Quantity, bool]
             dictionary of the form {(atom1_topology_index, atom2_topology_index): distance}
        """
        return self._constrained_atom_pairs

    @property
    def n_molecules(self) -> int:
        """Returns the number of molecules in this Topology"""
        return len(self._molecules)

    @property
    def molecules(self) -> Generator[MoleculeLike, None, None]:
        """Returns an iterator over all the Molecules in this Topology

        Returns
        -------
        molecules : Iterable of Molecule
        """
        # Yields instead of returning the list itself. This prevents user from modifying the list
        # outside the Topology's knowledge. This is essential to make sure that atom index caches
        # invalidate themselves during appropriate events.
        for molecule in self._molecules:
            yield molecule

    def molecule(self, index):
        """
        Returns the molecule with a given index in this Topology.

        Returns
        -------
        molecule : openff.toolkit.topology.Molecule
        """
        return self._molecules[index]

    @property
    def n_atoms(self) -> int:
        """
        Returns the number of  atoms in in this Topology.

        Returns
        -------
        n_atoms : int
        """
        n_atoms = 0
        for molecule in self.molecules:
            n_atoms += molecule.n_atoms
        return n_atoms

    @property
    def atoms(self) -> Generator["Atom", None, None]:
        """Returns an iterator over the atoms in this Topology. These will be in ascending order of topology index.

        Returns
        -------
        atoms : Generator of Atom
        """
        for molecule in self._molecules:
            for atom in molecule.atoms:
                yield atom

    def atom_index(self, atom: "Atom") -> int:
        """
        Returns the index of a given atom in this topology

        Parameters
        ----------
        atom : openff.toolkit.topology.Atom

        Returns
        -------
        index : int
            The index of the given atom in this topology

        Raises
        ------
        AtomNotInTopologyError
            If the given atom is not in this topology
        """
        # If the atom's topology atom index isn't cached, calculate it for the whole topology
        if "_topology_atom_index" not in atom.__dict__:
            self._build_atom_index_cache()
        # After computing the topology atom indices for all atoms in this topology, if this
        # atom still doesn't have a topology atom index assigned, then it must not be in this topology
        if "_topology_atom_index" not in atom.__dict__:
            raise AtomNotInTopologyError("Atom not found in this Topology")

        return atom._topology_atom_index  # type: ignore[attr-defined]

    def molecule_index(self, molecule: MoleculeLike) -> int:
        """
        Returns the index of a given molecule in this topology

        Parameters
        ----------
        molecule

        Returns
        -------
        index : int
            The index of the given molecule in this topology

        Raises
        ------
        MoleculeNotInTopologyError
            If the given atom is not in this topology
        """
        for index, iter_molecule in enumerate(self.molecules):
            if molecule is iter_molecule:
                return index

        raise MoleculeNotInTopologyError("Molecule not found in this Topology")

    def molecule_atom_start_index(self, molecule: Molecule) -> int:
        """
        Returns the index of a molecule's first atom in this topology

        Parameters
        ----------
        molecule

        Returns
        -------
        index : int
        """
        return self.atom_index(molecule.atoms[0])

    @property
    def n_bonds(self) -> int:
        """
        Returns the number of Bonds in in this Topology.

        Returns
        -------
        n_bonds : int
        """
        n_bonds = 0
        for molecule in self.molecules:
            n_bonds += molecule.n_bonds
        return n_bonds

    @property
    def bonds(self) -> Generator["Bond", None, None]:
        """Returns an iterator over the bonds in this Topology

        Returns
        -------
        bonds : Iterable of Bond
        """
        for molecule in self.molecules:
            for bond in molecule.bonds:
                yield bond

    @property
    def n_angles(self) -> int:
        """int: number of angles in this Topology."""
        return sum(mol.n_angles for mol in self._molecules)

    @property
    def angles(self) -> Generator[tuple["Atom", ...], None, None]:
        """Iterator over the angles in this Topology. Returns a Generator of tuple[Atom]."""
        for molecule in self._molecules:
            for angle in molecule.angles:
                yield angle

    @property
    def n_propers(self) -> int:
        """The number of proper torsions in this Topology."""
        return sum(mol.n_propers for mol in self._molecules)

    @property
    def propers(self) -> Generator[tuple[Union["Atom", _SimpleAtom], ...], None, None]:
        """Iterable of tuple[Atom]: iterator over the proper torsions in this Topology."""
        for molecule in self.molecules:
            for proper in molecule.propers:
                yield proper

    @property
    def n_impropers(self) -> int:
        """The number of possible improper torsions in this Topology."""
        return sum(mol.n_impropers for mol in self._molecules)

    @property
    def impropers(self) -> Generator[tuple["Atom", ...], None, None]:
        """Generator of tuple[Atom]: iterator over the possible improper torsions in this Topology."""
        for molecule in self._molecules:
            for improper in molecule.impropers:
                yield improper

    @property
    def smirnoff_impropers(
        self,
    ) -> Generator[tuple[Union["Atom", _SimpleAtom], ...], None, None]:
        """
        Iterate over improper torsions in the molecule, but only those with
        trivalent centers, reporting the central atom second in each improper.

        Note that it's possible that a trivalent center will not have an improper assigned.
        This will depend on the force field that is used.

        Also note that this will return 6 possible atom orderings around each improper
        center. In current SMIRNOFF parameterization, three of these six
        orderings will be used for the actual assignment of the improper term
        and measurement of the angles. These three orderings capture the three unique
        angles that could be calculated around the improper center, therefore the sum
        of these three terms will always return a consistent energy.

        For more details on the use of three-fold ('trefoil') impropers, see
        https://openforcefield.github.io/standards/standards/smirnoff/#impropertorsions

        .. todo:: Offer a way to do the keytransform and get the final 3 orderings in this
                  method? How can we keep this logic synced up with the parameterization
                  machinery?

        Returns
        -------
        smirnoff_impropers : Generator of tuples of Atom
            An iterator of tuples, each containing the Atom objects comprising
            up a possible improper torsion. The central atom is listed second in
            each tuple.

        See Also
        --------
        impropers, amber_impropers

        """
        for molecule in self.molecules:
            for smirnoff_improper in molecule.smirnoff_impropers:
                yield smirnoff_improper

    @property
    def amber_impropers(
        self,
    ) -> Generator[tuple[Union["Atom", _SimpleAtom], ...], None, None]:
        """
        Iterate over improper torsions in the molecule, but only those with
        trivalent centers, reporting the central atom first in each improper.

        Note that it's possible that a trivalent center will not have an improper assigned.
        This will depend on the force field that is used.

        Also note that this will return 6 possible atom orderings around each improper
        center. In current AMBER parameterization, one of these six
        orderings will be used for the actual assignment of the improper term
        and measurement of the angle. This method does not encode the logic to
        determine which of the six orderings AMBER would use.

        Returns
        -------
        amber_impropers : Generator of tuples of Atom
            An iterator of tuples, each containing the Atom objects comprising
            up a possible improper torsion. The central atom is listed first in
            each tuple.

        See Also
        --------
        impropers, smirnoff_impropers
        """
        for molecule in self.molecules:
            for amber_improper in molecule.amber_impropers:
                yield amber_improper

    def nth_degree_neighbors(self, n_degrees: int):
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
        for molecule in self.molecules:
            for pair in molecule.nth_degree_neighbors(n_degrees=n_degrees):
                yield pair

    class _ChemicalEnvironmentMatch:
        """Represents the match of a given chemical environment query, storing
        both the matched topology atom indices and the indices of the corresponding
        reference molecule atoms, as well as a reference to the reference molecule.
        """

        @property
        def reference_atom_indices(self):
            """tuple of int: The indices of the corresponding reference molecule atoms."""
            return self._reference_atom_indices

        @property
        def reference_molecule(self):
            """topology.molecule.Molecule: The corresponding reference molecule."""
            return self._reference_molecule

        @property
        def topology_atom_indices(self):
            """tuple of int: The matched topology atom indices."""
            return self._topology_atom_indices

        def __init__(
            self, reference_atom_indices, reference_molecule, topology_atom_indices
        ):
            """Constructs a new _ChemicalEnvironmentMatch object

            Parameters
            ----------
            reference_atom_indices: tuple of int
                The indices of the corresponding reference molecule atoms.
            reference_molecule: topology.molecule.Molecule
                The corresponding reference molecule.
            topology_atom_indices: tuple of int
                The matched topology atom indices.
            """

            assert len(reference_atom_indices) == len(topology_atom_indices)

            self._reference_atom_indices = reference_atom_indices
            self._reference_molecule = reference_molecule

            self._topology_atom_indices = topology_atom_indices

    def chemical_environment_matches(
        self,
        query: str,
        aromaticity_model: str = "MDL",
        unique: bool = False,
        toolkit_registry: TKR = GLOBAL_TOOLKIT_REGISTRY,
    ):
        """
        Retrieve all matches for a given chemical environment query.

        TODO:

        * Do we want to generalize this to other kinds of queries too, like mdtraj DSL, pymol selections, atom index
          slices, etc? We could just call it topology.matches(query)

        Parameters
        ----------
        query : str
            SMARTS string (with one or more tagged atoms)
        aromaticity_model : str
            Override the default aromaticity model for this topology and use the specified aromaticity model instead.
            Allowed values: ['MDL']

        Returns
        -------
        matches : list of Topology._ChemicalEnvironmentMatch
            A list of tuples, containing the topology indices of the matching atoms.

        """

        # Render the query to a SMARTS string
        if isinstance(query, str):
            smarts = query
        else:
            raise ValueError(
                f"Don't know how to convert query '{query}' into SMARTS string"
            )

        # Perform matching on each unique molecule, unrolling the matches to all matching copies
        # of that molecule in the Topology object.
        matches = list()

        groupings = self.identical_molecule_groups

        for unique_mol_idx, group in groupings.items():
            unique_mol = self.molecule(unique_mol_idx)
            # Find all atomsets that match this definition in the reference molecule
            # This will automatically attempt to match chemically identical atoms in
            # a canonical order within the Topology
            mol_matches = unique_mol.chemical_environment_matches(
                smarts,
                unique=unique,
                toolkit_registry=toolkit_registry,
            )

            if len(mol_matches) == 0:
                continue

            for mol_instance_idx, atom_map in group:
                mol_instance = self.molecule(mol_instance_idx)
                # Loop over matches
                for match in mol_matches:
                    # Collect indices of matching `Atom`s
                    topology_atom_indices = []
                    for molecule_atom_index in match:
                        atom = mol_instance.atom(atom_map[molecule_atom_index])
                        topology_atom_indices.append(self.atom_index(atom))

                    environment_match = Topology._ChemicalEnvironmentMatch(
                        tuple(match), unique_mol, tuple(topology_atom_indices)
                    )

                    matches.append(environment_match)
        return matches

    @property
    def identical_molecule_groups(self) -> dict[int, list[tuple[int, dict[int, int]]]]:
        """
        Returns groups of chemically identical molecules, identified by index and atom map.

        Returns
        -------
        identical_molecule_groups
            A dict of the form ``{unique_mol_idx : [(topology_mol_idx, atom_map), ...].``
            Each key is the topology molecule index of the first instance of a
            unique chemical species. Iterating over the keys will yield all of
            the unique chemical species in the topology. Each value is a list
            of all instances of that chemical species in the topology. Each
            element of the list is a 2-tuple where the first element is the
            molecule index of the instance, and the second element maps the atom
            topology indices of the key molecule to the instance. The molecule
            instance corresponding to the key is included in the list, so the
            list is a complete list of all instances of that chemical species.

        Examples
        --------

        >>> from openff.toolkit import Molecule, Topology
        >>> # Create a water ordered as OHH
        >>> water1 = Molecule()
        >>> water1.add_atom(8, 0, False)
        0
        >>> water1.add_atom(1, 0, False)
        1
        >>> water1.add_atom(1, 0, False)
        2
        >>> _ = water1.add_bond(0, 1, 1, False)
        >>> _ = water1.add_bond(0, 2, 1, False)
        ...
        >>> # Create a different water ordered as HOH
        >>> water2 = Molecule()
        >>> water2.add_atom(1, 0, False)
        0
        >>> water2.add_atom(8, 0, False)
        1
        >>> water2.add_atom(1, 0, False)
        2
        >>> _ = water2.add_bond(0, 1, 1, False)
        >>> _ = water2.add_bond(1, 2, 1, False)
        ...
        >>> top = Topology.from_molecules([water1, water2])
        >>> top.identical_molecule_groups
        {0: [(0, {0: 0, 1: 1, 2: 2}), (1, {0: 1, 1: 0, 2: 2})]}
        """
        # Check whether this was run previously, and a cached result is available.
        if self._cached_chemically_identical_molecules is not None:
            return self._cached_chemically_identical_molecules

        # Convert molecule identity maps into groups of identical molecules
        identity_maps = self._identify_chemically_identical_molecules()
        groupings: dict[int, list[tuple[int, dict[int, int]]]] = defaultdict(list)
        for molecule_idx, (unique_mol, atom_map) in identity_maps.items():
            groupings[unique_mol] += [(molecule_idx, atom_map)]

        self._cached_chemically_identical_molecules = dict(groupings)

        return self._cached_chemically_identical_molecules

    def _identify_chemically_identical_molecules(
        self,
    ) -> dict[int, tuple[int, dict[int, int]]]:
        """
        Perform an all-by-all isomorphism check over molecules in the Topology.

        Efficiently performs an all-by-all isomorphism check for the molecules in
        this Topology. This method uses the strictest form of isomorphism
        checking, which will NOT match distinct kekule structures of multiple
        resonance forms of the same molecule, or different kekulizations of
        aromatic systems.

        Returns
        -------
        identical_molecules : {int: (int, {int: int})}
            A mapping from the index of each molecule in the topology to (the
            index of the first appearance of a chemically equivalent molecule in
            the topology, and a mapping from the atom indices of this molecule
            to the atom indices of that chemically equivalent molecule).
            ``identical_molecules[molecule_idx] = (
                unique_molecule_idx, {molecule_atom_idx, unique_molecule_atom_idx}
            )``
        """
        identity_maps: dict[int, tuple[int, dict[int, int]]] = dict()
        already_matched_mols = set()

        for mol1_idx in range(self.n_molecules):
            if mol1_idx in already_matched_mols:
                continue
            mol1 = self.molecule(mol1_idx)
            identity_maps[mol1_idx] = (
                mol1_idx,
                {i: i for i in range(mol1.n_atoms)},
            )
            for mol2_idx in range(mol1_idx + 1, self.n_molecules):
                if mol2_idx in already_matched_mols:
                    continue
                mol2 = self.molecule(mol2_idx)
                if isinstance(mol1, type(mol2)) or isinstance(mol2, type(mol1)):
                    are_isomorphic, atom_map = mol1.are_isomorphic(
                        mol1, mol2, return_atom_map=True
                    )
                else:
                    are_isomorphic = False

                if are_isomorphic:
                    identity_maps[mol2_idx] = (
                        mol1_idx,
                        atom_map,
                    )
                    already_matched_mols.add(mol2_idx)

        return identity_maps

    def _build_atom_index_cache(self):
        topology_molecule_atom_start_index = 0
        for molecule in self.molecules:
            for at in molecule.atoms:
                at._topology_atom_index = (
                    topology_molecule_atom_start_index + at.molecule_atom_index
                )
            topology_molecule_atom_start_index += molecule.n_atoms

    def _invalidate_cached_properties(self):
        self._cached_chemically_identical_molecules = None
        for atom in self.atoms:
            if "_topology_atom_index" in atom.__dict__:
                del atom.__dict__["_topology_atom_index"]

    def copy_initializer(self, other):
        import copy

        self.aromaticity_model = other.aromaticity_model
        self._constrained_atom_pairs = copy.deepcopy(other._constrained_atom_pairs)
        self._box_vectors = copy.deepcopy(other._box_vectors)
        self._molecules = copy.deepcopy(other._molecules)
        self._invalidate_cached_properties()

    def to_dict(self):
        from openff.toolkit.utils.utils import serialize_numpy

        return_dict = dict()
        return_dict["aromaticity_model"] = self._aromaticity_model
        return_dict["constrained_atom_pairs"] = dict()
        for constrained_atom_pair, distance in self._constrained_atom_pairs.items():
            return_dict["constrained_atom_pairs"][constrained_atom_pair] = (
                quantity_to_string(distance)
            )

        if self._box_vectors is None:
            return_dict["box_vectors"] = None
            return_dict["box_vectors_unit"] = None
        else:
            box_vectors_unitless = self.box_vectors.m_as(unit.nanometer)
            box_vectors_serialized, box_vectors_shape = serialize_numpy(
                box_vectors_unitless
            )
            if box_vectors_shape != (3, 3):
                raise RuntimeError(
                    f"Box vectors are assumed to be (3, 3); found shape {box_vectors_shape}"
                )
            return_dict["box_vectors"] = box_vectors_serialized
            return_dict["box_vectors_unit"] = "nanometer"
        return_dict["molecules"] = [mol.to_dict() for mol in self._molecules]
        return return_dict

    @classmethod
    def from_dict(cls, topology_dict: dict):
        """
        Create a new Topology from a dictionary representation

        Parameters
        ----------
        topology_dict : dict
            A dictionary representation of the topology.

        Returns
        -------
        topology : Topology
            A Topology created from the dictionary representation

        """
        topology = cls()
        topology._initialize_from_dict(topology_dict)
        return topology

    def _initialize_from_dict(self, topology_dict):
        from openff.toolkit.utils.utils import deserialize_numpy

        self._aromaticity_model = topology_dict["aromaticity_model"]
        for pair, distance in topology_dict["constrained_atom_pairs"]:
            deserialized_distance = string_to_quantity(distance)
            self.add_constraint(pair, deserialized_distance)

        if topology_dict["box_vectors"] is None:
            self._box_vectors = None
        else:
            # The box_vectors setters _should_ ensure (3, 3) shape, and
            # to_dict enforces this at serialization time
            box_vectors_unitless = deserialize_numpy(
                topology_dict["box_vectors"],
                (3, 3),
            )
            box_vectors_unit = getattr(unit, topology_dict["box_vectors_unit"])
            self.box_vectors = Quantity(box_vectors_unitless, box_vectors_unit)

        for molecule_dict in topology_dict["molecules"]:
            new_mol = Molecule.from_dict(molecule_dict)
            self._add_molecule_keep_cache(new_mol)
        self._invalidate_cached_properties()

    @staticmethod
    @requires_package("openmm")
    def _openmm_topology_to_networkx(openmm_topology):
        import networkx as nx

        # Convert all openMM mols to graphs
        omm_topology_G = nx.Graph()
        for atom in openmm_topology.atoms():
            if atom.element is None:
                raise VirtualSitesUnsupportedError(
                    f"Atom {atom} has element None and is probably a virtual site. Virtual sites "
                    "cannot be stored in `Topology` objects but can be stored in `Interchange` "
                    "objects produced by `ForceField`s with virtual site parameters."
                )
            omm_topology_G.add_node(
                atom.index,
                atomic_number=atom.element.atomic_number,
                atom_name=atom.name,
                residue_name=atom.residue.name,
                residue_id=atom.residue.id,
                insertion_code=atom.residue.insertionCode,
                chain_id=atom.residue.chain.id,
            )
        for bond in openmm_topology.bonds():
            omm_topology_G.add_edge(
                bond.atom1.index, bond.atom2.index, bond_order=bond.order
            )
        return omm_topology_G

    @classmethod
    @requires_package("openmm")
    def from_openmm(
        cls,
        openmm_topology: "openmm.app.Topology",
        unique_molecules: Optional[Iterable[FrozenMolecule]] = None,
        positions: Union[None, Quantity, "OMMQuantity"] = None,
    ) -> "Topology":
        """
        Construct an OpenFF Topology object from an OpenMM Topology object.

        This method guarantees that the order of atoms in the input OpenMM
        Topology will be the same as the ordering of atoms in the output OpenFF
        Topology. However it does not guarantee the order of the bonds will be
        the same.

        Hierarchy schemes are taken from the OpenMM topology, not from
        ``unique_molecules``.

        If any virtual sites are detected in the OpenMM topology,
        ``VirtualSitesUnsupportedError`` is raised because the ``Topology``
        object model does not store virtual sites.

        Parameters
        ----------
        openmm_topology
            The OpenMM Topology object to convert
        unique_molecules
            An iterable containing all the unique molecules in the topology.
            This is used to identify the molecules in the OpenMM topology and
            provide any missing chemical information. Each chemical species in
            the topology must be specified exactly once, though the topology
            may have any number of copies, including zero. The chemical
            elements of atoms and their bond connectivity will be used to match
            these reference molecules to the molecules appearing in the
            topology. If bond orders are specified in the topology, these will
            be used in matching as well.
        positions
            Positions for the atoms in the new topology.

        Returns
        -------
        topology
            An OpenFF Topology object, constructed from the molecules in
            ``unique_molecules``, with the same atom order as the input topology.

        Raises
        ------
        MissingUniqueMoleculesError
            If ``unique_molecules`` is ``None``
        DuplicateUniqueMoleculeError
            If the same connectivity graph is represented by two different
            molecules in ``unique_molecules``
        ValueError
            If a chemically impossible molecule is detected in the topology
        """
        import networkx as nx
        from openff.units.openmm import from_openmm

        from openff.toolkit.topology.molecule import Molecule

        # Check to see if the openMM system has defined bond orders, by looping over all Bonds in the Topology.
        omm_has_bond_orders = True
        for omm_bond in openmm_topology.bonds():
            if omm_bond.order is None:
                omm_has_bond_orders = False

        if unique_molecules is None:
            raise MissingUniqueMoleculesError(
                "Topology.from_openmm requires a list of Molecule objects "
                "passed as unique_molecules, but None was passed."
            )

        # Convert all unique mols to graphs
        topology = cls()
        graph_to_unq_mol: dict[Graph, FrozenMolecule] = {}
        for unq_mol in unique_molecules:
            unq_mol_graph = unq_mol.to_networkx()
            for existing_graph in graph_to_unq_mol.keys():
                if Molecule.are_isomorphic(
                    existing_graph,
                    unq_mol_graph,
                    return_atom_map=False,
                    aromatic_matching=False,
                    formal_charge_matching=False,
                    bond_order_matching=omm_has_bond_orders,
                    atom_stereochemistry_matching=False,
                    bond_stereochemistry_matching=False,
                )[0]:
                    msg = (
                        "Error: Two unique molecules have indistinguishable "
                        "graphs: {} and {}".format(
                            unq_mol, graph_to_unq_mol[existing_graph]
                        )
                    )
                    raise DuplicateUniqueMoleculeError(msg)
            graph_to_unq_mol[unq_mol_graph] = unq_mol

        omm_topology_G = cls._openmm_topology_to_networkx(openmm_topology)
        # For each connected subgraph (molecule) in the topology, find its match in unique_molecules
        topology_molecules_to_add = list()
        for omm_mol_G in (
            omm_topology_G.subgraph(c).copy()
            for c in nx.connected_components(omm_topology_G)
        ):
            match_found = False
            for unq_mol_G in graph_to_unq_mol.keys():
                isomorphic, mapping = Molecule.are_isomorphic(
                    omm_mol_G,
                    unq_mol_G,
                    return_atom_map=True,
                    aromatic_matching=False,
                    formal_charge_matching=False,
                    bond_order_matching=omm_has_bond_orders,
                    atom_stereochemistry_matching=False,
                    bond_stereochemistry_matching=False,
                )
                if isomorphic:
                    # Take the first valid atom indexing map
                    first_topology_atom_index = min(mapping.keys())  # type: ignore[union-attr]
                    topology_molecules_to_add.append(
                        (
                            first_topology_atom_index,
                            unq_mol_G,
                            mapping.items(),  # type: ignore[union-attr]
                            omm_mol_G,
                        )
                    )
                    match_found = True
                    break
            if match_found is False:
                from openff.toolkit.topology.molecule import (
                    _networkx_graph_to_hill_formula,
                )

                hill_formula = _networkx_graph_to_hill_formula(omm_mol_G)
                msg = f"No match found for molecule {hill_formula}. "
                probably_missing_conect = [
                    "C",
                    "H",
                    "O",
                    "N",
                    "P",
                    "S",
                    "F",
                    "Cl",
                    "Br",
                ]
                if hill_formula in probably_missing_conect:
                    msg += (
                        "This would be a very unusual molecule to try and parameterize, "
                        "and it is likely that the data source it was read from does not "
                        "contain connectivity information. If this molecule is coming from "
                        "PDB, please ensure that the file contains CONECT records. The PDB "
                        "format documentation (https://www.wwpdb.org/documentation/"
                        'file-format-content/format33/sect10.html) states "CONECT records '
                        "are mandatory for HET groups (excluding water) and for other bonds "
                        'not specified in the standard residue connectivity table."'
                    )
                raise ValueError(msg)

        # The connected_component_subgraph function above may have scrambled the molecule order, so sort molecules
        # by their first atom's topology index
        topology_molecules_to_add.sort(key=lambda x: x[0])
        for (
            first_index,
            unq_mol_G,
            top_to_ref_index,
            omm_mol_G,
        ) in topology_molecules_to_add:
            local_top_to_ref_index = dict(
                [
                    (top_index - first_index, ref_index)
                    for top_index, ref_index in top_to_ref_index
                ]
            )
            unq_mol = graph_to_unq_mol[unq_mol_G]
            remapped_mol = unq_mol.remap(local_top_to_ref_index, current_to_new=False)
            # Transfer hierarchy metadata from openmm mol graph to offmol metadata
            for off_atom_idx in range(remapped_mol.n_atoms):
                off_atom = remapped_mol.atom(off_atom_idx)

                omm_atom_idx = off_atom_idx + first_index

                off_atom.name = omm_mol_G.nodes[omm_atom_idx]["atom_name"]
                off_atom.metadata["residue_name"] = omm_mol_G.nodes[omm_atom_idx][
                    "residue_name"
                ]
                off_atom.metadata["residue_number"] = int(
                    omm_mol_G.nodes[omm_atom_idx]["residue_id"]
                )
                off_atom.metadata["insertion_code"] = omm_mol_G.nodes[omm_atom_idx][
                    "insertion_code"
                ]

                off_atom.metadata["chain_id"] = omm_mol_G.nodes[omm_atom_idx][
                    "chain_id"
                ]

            remapped_mol.add_default_hierarchy_schemes()
            topology._add_molecule_keep_cache(remapped_mol)
        topology._invalidate_cached_properties()

        if openmm_topology.getPeriodicBoxVectors() is not None:
            topology.box_vectors = from_openmm(openmm_topology.getPeriodicBoxVectors())

        if positions is not None:
            topology.set_positions(ensure_quantity(positions, "openff"))

        # TODO: How can we preserve metadata from the openMM topology when creating the OFF topology?
        return topology

    def _ensure_unique_atom_names(self, ensure_unique_atom_names: Union[str, bool]):
        """See `Topology.to_openmm`"""
        if not ensure_unique_atom_names:
            return
        for molecule in self._molecules:
            if isinstance(ensure_unique_atom_names, str) and hasattr(
                molecule, ensure_unique_atom_names
            ):
                for hier_elem in getattr(molecule, ensure_unique_atom_names):
                    if not hier_elem.has_unique_atom_names:
                        hier_elem.generate_unique_atom_names()
            elif not molecule.has_unique_atom_names:
                molecule.generate_unique_atom_names()

    @classmethod
    @requires_package("openmm")
    def from_pdb(
        cls,
        file_path: Union[str, Path, TextIO],
        unique_molecules: Optional[Iterable[Molecule]] = None,
        toolkit_registry=GLOBAL_TOOLKIT_REGISTRY,
        _custom_substructures: Optional[dict[str, list[str]]] = None,
        _additional_substructures: Optional[Iterable[Molecule]] = None,
    ):
        """
        Loads supported or user-specified molecules from a PDB file.

        Currently, canonical proteins, waters, and monoatomic ions are supported
        without CONECT records via residue and atom names, and molecules
        specified in the ``unique_molecules`` argument are supported when
        CONECT records are provided.

        .. warning::
            Molecules in the resulting Topology will adopt
            the geometric stereochemistry in the PDB, even if this conflicts
            with the stereochemistry specified in ``unique_molecules``.

        All molecules in the PDB file have the following requirements:

        * Polymer molecules must use the standard atom names described in the
          `PDB Chemical Component dictionary <https://www.wwpdb.org/data/ccd>`_
          (PDB CCD).
        * There must be no missing atoms (all hydrogens must be explicit).
        * All particles must correspond to an atomic nucleus (particles in the
          PDB representing "virtual sites" or "extra points" are not allowed).
        * All bonds must be specified by either CONECT records, or for
          polymers and water by the PDB CCD via the residue and atom name.
        * CONECT records must correspond only to chemical bonds (CONECT records
          representing an angle constraints are not allowed).
        * CONECT records may be redundant with connectivity defined by residue
          templates.

        Currently, the only supported polymers are proteins made of the 20
        canonical amino acids. For details on the polymer loading used here,
        see :py:meth:`Molecule.from_polymer_pdb`.

        Waters can be recognized in either of two ways:

        * By the residue name "HOH" and atom names "H1", "H2", and "O".
        * By ATOM records which include element information and CONECT records.

        Monoatomic ions supported by Sage are recognized (Na+, Li+, K+, Rb+,
        Cs+, F-, Cl-, Br-, and I-). To load other monoatomic ions, use the
        ``unique_molecules`` keyword argument.

        The ``unique_molecules`` keyword argument can be used to load arbitrary
        molecules from the PDB file. These molecules match a group of atoms in
        the PDB file when their atomic elements and connectivity are identical;
        elements and CONECT records must therefore be explicitly specified in
        the PDB file. Information missing from the PDB format, such as bond
        orders and formal charges, is then taken from the matching unique
        molecule. Unique molecule matches will overwrite bond order and formal
        charge assignments from other sources. Stereochemistry is assigned
        based on the geometry in the PDB, even if this differs from the
        stereochemistry in the unique molecule.

        A user-defined molecule in the PDB file must exactly match a unique
        molecule to successfully load it - substructures and superstructures
        will raise :py:exc:`UnassignedChemistryInPDBError`. Unique molecules
        need not be present in the PDB.

        Metadata such as residues, chains, and atom names are recorded in the
        ``Atom.metadata`` attribute, which is a dictionary mapping from the
        strings ``"residue_name"``, ``"residue_number"``, ``"insertion_code"``,
        and ``"chain_id"`` to the appropriate value. The topology returned by
        this method can expose residue and chain iterators which can be
        accessed using :py:meth:`Topology.hierarchy_iterator`, such as
        ``top.hierarchy_iterator("residues")`` and ``top.hierarchy_iterator
        ("chains")``.

        A CRYST line in the PDB, if present, will be interpreted as periodic
        box vectors in Angstroms.

        Parameters
        ----------
        file_path : str, Path, or file object
            PDB information to be passed to OpenMM PDBFile object for loading
        unique_molecules : Iterable of Molecule. Default = None
            OpenFF Molecule objects corresponding to the molecules in the input
            PDB. See above for details.
        toolkit_registry : ToolkitRegistry. Default = None
            The ToolkitRegistry to use as the cheminformatics backend.
        _custom_substructures: dict[str, list[str]], Default = {}
            Experimental and unstable. Dictionary where keys are the names of new substructures
            (cannot overlap with existing amino acid names) and the values are the new substructure
            entries that follow the same format as those used in the amino acid substructure library
        _additional_substructures : Iterable of Molecule, Default = None
            Experimental and unstable. Molecule with atom.metadata["substructure_atom"] =
            True or False for all atoms. Currently only stable for independent, standalone
            molecules not bonded to a larger protein/molecule. (For that use _custom_substructures)

        Returns
        -------
        topology : openff.toolkit.topology.Topology

        Raises
        ------

        UnassignedChemistryInPDBError
            If an atom or bond could not be assigned; the exception will
            provide a detailed diagnostic of what went wrong.

        Examples
        --------

        >>> from openff.toolkit import Topology
        >>> from openff.toolkit.utils import get_data_file_path
        >>> top = Topology.from_pdb(get_data_file_path("proteins/TwoMol_SER_CYS.pdb"))
        >>> # The molecules in the loaded topology are full-fledged OpenFF Molecule objects
        >>> for match in top.chemical_environment_matches('[O:1]=[C:2][N:3][H:4]'): print(match.topology_atom_indices)
        (1, 0, 6, 13)
        (9, 8, 17, 19)
        (24, 23, 29, 36)
        (32, 31, 40, 42)
        >>> [*top.hierarchy_iterator("residues")]
        [HierarchyElement ('A', '1', ' ', 'ACE') of iterator 'residues' containing 6 atom(s),
         HierarchyElement ('A', '2', ' ', 'SER') of iterator 'residues' containing 11 atom(s),
         HierarchyElement ('A', '3', ' ', 'NME') of iterator 'residues' containing 6 atom(s),
         HierarchyElement ('B', '1', ' ', 'ACE') of iterator 'residues' containing 6 atom(s),
         HierarchyElement ('B', '2', ' ', 'CYS') of iterator 'residues' containing 11 atom(s),
         HierarchyElement ('B', '3', ' ', 'NME') of iterator 'residues' containing 6 atom(s)]

         Polymer systems can also be supported if _custom_substructures are given as a dict[str, list[str]],
         where the keys are unique atom names and the values are lists of substructure smarts. The
         substructure smarts must follow the same format as given in
         "proteins/aa_residues_substructures_explicit_bond_orders_with_caps_explicit_connectivity.json":
         ”<bond>[#<atomic number>D<degree>+<formal charge>:<id>]<bond>” for monomer atoms and
         ”<bond>[*:<id>]” for adjacent neighboring atoms
         (NOTE: This functionality is experimental!)

        >>> PE_substructs = {
        ...     "PE": [
        ...         "[#6D4+0:2](-[#1D1+0:3])(-[#1D1+0:4])(-[#6D4+0:5](-[#1D1+0:6])(-[#1D1+0:7])-[*:8])-[*:1]",
        ...         "[#6D4+0:2](-[#1D1+0:3])(-[#1D1+0:4])(-[#6D4+0:5](-[#1D1+0:6])(-[#1D1+0:7])-[#1D1+0:8])-[*:1]",
        ...         "[#6D4+0:2](-[#1D1+0:3])(-[#1D1+0:4])(-[#6D4+0:5](-[#1D1+0:6])(-[#1D1+0:7])-[*:8])-[#1D1+0:1]",
        ...     ]
        ... }
        >>> top = Topology.from_pdb(
        ...          get_data_file_path("systems/test_systems/PE.pdb"),
        ...          _custom_substructures=PE_substructs,
        ...      )
        """
        import io
        import json

        import openmm.unit as openmm_unit
        from openmm.app import PDBFile

        if isinstance(file_path, (str, io.TextIOWrapper)):
            pass
        elif isinstance(file_path, Path):
            file_path = file_path.as_posix()
        else:
            raise ValueError(f"Unexpected type {type(file_path)}")

        if not _custom_substructures:
            _custom_substructures = dict()

        pdb = PDBFile(file_path)

        substructure_file_path = get_data_file_path(
            "proteins/aa_residues_substructures_explicit_bond_orders_with_caps_explicit_connectivity.json"
        )

        with open(substructure_file_path, "r") as subfile:
            substructure_dictionary = json.load(
                subfile
            )  # preserving order is useful later when saving metadata
        substructure_dictionary["HOH"] = {"[H:1][O:2][H:3]": ["H1", "O", "H2"]}
        substructure_dictionary["Li"] = {"[Li+1:1]": ["Li"]}
        substructure_dictionary["Na"] = {"[Na+1:1]": ["Na"]}
        substructure_dictionary["K"] = {"[K+1:1]": ["K"]}
        substructure_dictionary["Rb"] = {"[Rb+1:1]": ["Rb"]}
        substructure_dictionary["Cs"] = {"[Cs+1:1]": ["Cs"]}
        substructure_dictionary["F"] = {"[F-1:1]": ["F"]}
        substructure_dictionary["Cl"] = {"[Cl-1:1]": ["Cl"]}
        substructure_dictionary["Br"] = {"[Br-1:1]": ["Br"]}
        substructure_dictionary["I"] = {"[I-1:1]": ["I"]}

        if not (unique_molecules):
            unique_molecules = []

        substructure_dictionary["UNIQUE_MOLECULE"] = {}

        for unique_molecule in unique_molecules:
            mapped_smiles = unique_molecule.to_smiles(mapped=True)
            substructure_dictionary["UNIQUE_MOLECULE"][mapped_smiles] = [
                a.name for a in unique_molecule.atoms
            ]

        substructure_dictionary["ADDITIONAL_SUBSTRUCTURE"] = {}

        if _additional_substructures:
            for mol in _additional_substructures:
                label_mol = Molecule(mol)
                c = 0
                label_mol.properties["atom_map"] = {}
                for atom in label_mol.atoms:
                    if atom.metadata["substructure_atom"]:
                        label_mol.properties["atom_map"][atom.molecule_atom_index] = c
                        c += 1
                smi = label_mol.to_smiles(mapped=True)
                # remove unmapped atoms from mapped smiles. This will catch things like
                # `[H]` and `[Cl]` but not anything with 3 characters like `[H:1]`
                smi = re.sub("\[[A-Za-z]{1,2}\]", "[*]", smi)
                # Remove any orphaned () that remain
                smi = smi.replace("()", "")

                substructure_dictionary["ADDITIONAL_SUBSTRUCTURE"][smi] = []

        substructure_dictionary["ADDITIONAL_SUBSTRUCTURE_OVERLAP"] = {}

        coords_angstrom = np.array(
            [[*vec3.value_in_unit(openmm_unit.angstrom)] for vec3 in pdb.getPositions()]
        )

        topology = toolkit_registry.call(
            "_polymer_openmm_pdbfile_to_offtop",
            cls,
            pdb,
            substructure_dictionary,
            coords_angstrom,
            _custom_substructures,
        )

        for off_atom, atom in zip([*topology.atoms], pdb.topology.atoms()):
            off_atom.metadata["residue_name"] = atom.residue.name
            off_atom.metadata["residue_number"] = atom.residue.id
            off_atom.metadata["insertion_code"] = atom.residue.insertionCode
            off_atom.metadata["chain_id"] = atom.residue.chain.id
            off_atom.name = atom.name

        for offmol in topology.molecules:
            offmol.add_default_hierarchy_schemes()

        return topology

    @requires_package("openmm")
    def to_openmm(self, ensure_unique_atom_names: Union[str, bool] = "residues"):
        """
        Create an OpenMM Topology object.

        The atom metadata fields `residue_name`, `residue_number`, `insertion_code`, and `chain_id`
        are used to group atoms into OpenMM residues and chains.

        Contiguously-indexed atoms with the same `residue_name`, `residue_number`, `insertion_code`,
        and `chain_id` will be put into the same OpenMM residue.

        Contiguously-indexed residues with with the same `chain_id` will be put
        into the same OpenMM chain.

        This method will never make an OpenMM chain or residue larger than the
        OpenFF Molecule that it came from. In other words, no chain or residue
        will span two OpenFF Molecules.

        This method will **not** populate the OpenMM Topology with virtual sites.

        Parameters
        ----------
        ensure_unique_atom_names
            Whether to generate new atom names to ensure uniqueness within a
            molecule or hierarchy element.

            - If the name of a :class:`HierarchyScheme` is given as a string,
              new atom names will be generated so that each element of that
              scheme has unique atom names. Molecules without the given
              hierarchy scheme will be given unique atom names within that
              molecule.
            - If ``True``, new atom names will be generated so that atom names
              are unique within a molecule.
            - If ``False``, the existing atom names will be used.

        Returns
        -------
        openmm_topology : openmm.app.Topology
            An OpenMM Topology object
        """
        # TODO: MT needs to write a virtual sites section of the Interchange user guide.
        #       Once that exists, the last note in this docstring should link to that.
        from openmm import app

        from openff.toolkit.topology.molecule import Bond

        off_topology = Topology(self)
        omm_topology = app.Topology()

        off_topology._ensure_unique_atom_names(ensure_unique_atom_names)

        # Go through atoms in OpenFF to preserve the order.
        omm_atoms = []

        # For each atom in each molecule, determine which chain/residue it should be a part of
        for molecule in off_topology.molecules:
            # No chain or residue can span more than one OFF molecule, so reset these to None for the first
            # atom in each molecule.
            last_chain = None
            last_residue = None
            for atom in molecule.atoms:
                # If the residue name is undefined, assume a default of "UNK"
                if "residue_name" in atom.metadata:
                    atom_residue_name = atom.metadata["residue_name"]
                else:
                    atom_residue_name = "UNK"

                # If the residue number is undefined, assume a default of "0"
                if "residue_number" in atom.metadata:
                    atom_residue_number = atom.metadata["residue_number"]
                else:
                    atom_residue_number = "0"

                # If the insertion code  is undefined, assume a default of " "
                if "insertion_code" in atom.metadata:
                    atom_insertion_code = atom.metadata["insertion_code"]
                else:
                    atom_insertion_code = " "

                # If the chain ID is undefined, assume a default of "X"
                if "chain_id" in atom.metadata:
                    atom_chain_id = atom.metadata["chain_id"]
                else:
                    atom_chain_id = "X"

                # Determine whether this atom should be part of the last atom's chain, or if it
                # should start a new chain
                if last_chain is None:
                    chain = omm_topology.addChain(atom_chain_id)
                elif last_chain.id == atom_chain_id:
                    chain = last_chain
                else:
                    chain = omm_topology.addChain(atom_chain_id)
                # Determine whether this atom should be a part of the last atom's residue, or if it
                # should start a new residue
                if last_residue is None:
                    residue = omm_topology.addResidue(
                        atom_residue_name,
                        chain,
                        id=atom_residue_number,
                        insertionCode=atom_insertion_code,
                    )
                elif (
                    (last_residue.name == atom_residue_name)
                    and (int(last_residue.id) == int(atom_residue_number))
                    and (last_residue.insertionCode == atom_insertion_code)
                    and (chain.id == last_chain.id)
                ):
                    residue = last_residue
                else:
                    residue = omm_topology.addResidue(
                        atom_residue_name,
                        chain,
                        id=atom_residue_number,
                        insertionCode=atom_insertion_code,
                    )

                # Add atom.
                element = app.Element.getByAtomicNumber(atom.atomic_number)
                omm_atom = omm_topology.addAtom(atom.name, element, residue)

                # Make sure that OpenFF and OpenMM Topology atoms have the same indices.
                assert off_topology.atom_index(atom) == int(omm_atom.id) - 1
                omm_atoms.append(omm_atom)

                last_chain = chain
                last_residue = residue

            # Add all bonds.
            bond_types = {1: app.Single, 2: app.Double, 3: app.Triple}
            for bond in molecule.bonds:
                atom1, atom2 = bond.atoms
                atom1_idx, atom2_idx = off_topology.atom_index(
                    atom1
                ), off_topology.atom_index(atom2)
                if isinstance(bond, Bond):
                    if bond.is_aromatic:
                        bond_type = app.Aromatic
                    else:
                        bond_type = bond_types[bond.bond_order]
                    bond_order = bond.bond_order
                elif isinstance(bond, _SimpleBond):
                    bond_type = None
                    bond_order = None
                else:
                    raise RuntimeError(
                        "Unexpected bond type found while iterating over Topology.bonds."
                        f"Found {type(bond)}, allowed are Bond and _SimpleBond."
                    )

                omm_topology.addBond(
                    omm_atoms[atom1_idx],
                    omm_atoms[atom2_idx],
                    type=bond_type,
                    order=bond_order,
                )

        if off_topology.box_vectors is not None:
            from openff.units.openmm import to_openmm

            omm_topology.setPeriodicBoxVectors(to_openmm(off_topology.box_vectors))
        return omm_topology

    @requires_package("openmm")
    def to_file(
        self,
        file: Union[Path, str, TextIO],
        positions: Optional[Union["OMMQuantity", Quantity, NDArray]] = None,
        file_format: Literal["PDB"] = "PDB",
        keep_ids: bool = False,
        ensure_unique_atom_names: Union[str, bool] = "residues",
    ):
        """
        Save coordinates and topology to a PDB file.

        Reference: https://github.com/openforcefield/openff-toolkit/issues/502

        Notes:

        1. Atom numbering may not remain same, for example if the atoms
           in water are numbered as 1001, 1002, 1003, they would change
           to 1, 2, 3. This doesn't affect the topology or coordinates or
           atom-ordering in any way.
        2. Same issue with the amino acid names in the pdb file, they are
           not returned.

        Parameters
        ----------
        file
            A file-like object to write to, or a path to save the file to.
        positions : Array with shape ``(n_atoms, 3)`` and dimensions of length
            May be a...

            - ``openmm.unit.Quantity`` object which has atomic positions as a
              List of unit-tagged ``Vec3`` objects
            - ``openff.units.unit.Quantity`` object which wraps a
              ``numpy.ndarray`` with dimensions of length
            - (unitless) 2D ``numpy.ndarray``, in which it is assumed that the
              positions are in units of Angstroms.
            - ``None`` (the default), in which case the first conformer of
              each molecule in the topology will be used.

        file_format
            Output file format. Case insensitive. Currently only supported values
            are ``"PDB"``.
        keep_ids
            If ``True``, keep the residue and chain IDs specified in the Topology
            rather than generating new ones.
        ensure_unique_atom_names
            Whether to generate new atom names to ensure uniqueness within a
            molecule or hierarchy element.

            - If the name of a :class:`HierarchyScheme` is given as a string,
              new atom names will be generated so that each element of that
              scheme has unique atom names. Molecules without the given
              hierarchy scheme will be given unique atom names within that
              molecule.
            - If ``True``, new atom names will be generated so that atom names
              are unique within a molecule.
            - If ``False``, the existing atom names will be used.

            Note that this option cannot guarantee name uniqueness for formats
            like PDB that truncate long atom names.

        """
        import openmm.app
        import openmm.unit

        # Convert the topology to OpenMM
        openmm_top = self.to_openmm(ensure_unique_atom_names=ensure_unique_atom_names)

        # Get positions in OpenMM format
        if isinstance(positions, openmm.unit.Quantity):
            openmm_positions: openmm.unit.Quantity = positions
        elif isinstance(positions, Quantity):
            openmm_positions = positions.to_openmm()
        elif isinstance(positions, np.ndarray):
            openmm_positions = openmm.unit.Quantity(positions, openmm.unit.angstroms)
        elif positions is None:
            openmm_positions = self.get_positions().to_openmm()  # type: ignore[union-attr]
        else:
            raise ValueError(f"Could not process positions of type {type(positions)}.")

        if file_format.upper() == "PDB":
            import openmm

            # Convert the topology to OpenMM
            openmm_top = self.to_openmm(
                ensure_unique_atom_names=ensure_unique_atom_names
            )

            # Write PDB file
            ctx_manager: Union[nullcontext[TextIO], TextIO]  # MyPy needs some help here
            if isinstance(file, (str, Path)):
                ctx_manager = open(file, "w")
            else:
                ctx_manager = nullcontext(file)
            with ctx_manager as outfile:
                openmm.app.PDBFile.writeFile(
                    topology=openmm_top,
                    positions=openmm_positions,
                    file=outfile,
                    keepIds=keep_ids,
                )

        else:
            raise NotImplementedError("Topology.to_file supports only PDB")

    def get_positions(self) -> Optional[Quantity]:
        """
        Copy the positions of the topology into a new array.

        Topology positions are stored as the first conformer of each molecule.
        If any molecule has no conformers, this method returns ``None``. Note
        that modifying the returned array will not update the positions in the
        topology. To change the positions, use :meth:`Topology.set_positions`.

        See Also
        ========
        set_positions
        """
        conformers = []
        for molecule in self.molecules:
            try:
                conformer = molecule.conformers[0]
            except (IndexError, TypeError):
                return None

            conformer = conformer.m_as(unit.nanometer)

            conformers.append(conformer)
        positions = np.concatenate(conformers, axis=0)

        return Quantity(positions, unit.nanometer)

    def set_positions(self, array: Quantity):
        """
        Set the positions in a topology by copying from a single n×3 array.

        Note that modifying the original array will not update the positions
        in the topology; it must be passed again to ``set_positions()``.

        Parameters
        ==========

        array
            Positions for the topology. Should be a unit-wrapped array-like
            object with shape (n_atoms, 3) and dimensions of length.

        See Also
        ========
        get_positions
        """
        if not isinstance(array, Quantity):
            raise IncompatibleUnitError(
                "array should be an OpenFF Quantity with dimensions of length"
            )

        # Copy the array in nanometers and make it an OpenFF Quantity
        array = Quantity(np.asarray(array.to(unit.nanometer).magnitude), unit.nanometer)
        if array.shape != (self.n_atoms, 3):
            raise WrongShapeError(
                f"Array has shape {array.shape} but should have shape {self.n_atoms, 3}"
            )

        start = 0
        for molecule in self.molecules:
            stop = start + molecule.n_atoms
            if molecule.conformers is None:
                if isinstance(molecule, Molecule):
                    molecule._conformers = [array[start:stop]]
                else:
                    molecule.conformers = [array[start:stop]]  # type: ignore[misc]
            else:
                molecule.conformers[0:1] = [array[start:stop]]
            start = stop

    @classmethod
    @requires_package("mdtraj")
    def from_mdtraj(
        cls,
        mdtraj_topology: "mdtraj.Topology",
        unique_molecules: Optional[Iterable[MoleculeLike]] = None,
        positions: Union[None, "OMMQuantity", Quantity] = None,
    ):
        """
        Construct an OpenFF ``Topology`` from an MDTraj ``Topology``

        This method guarantees that the order of atoms in the input MDTraj
        Topology will be the same as the ordering of atoms in the output OpenFF
        Topology. However it does not guarantee the order of the bonds will be
        the same.

        Hierarchy schemes are taken from the MDTraj topology, not from
        ``unique_molecules``.

        Parameters
        ----------
        mdtraj_topology
            The MDTraj Topology object to convert
        unique_molecules
            An iterable containing all the unique molecules in the topology.
            This is used to identify the molecules in the MDTraj topology and
            provide any missing chemical information. Each chemical species in
            the topology must be specified exactly once, though the topology
            may have any number of copies, including zero. The chemical
            elements of atoms and their bond connectivity will be used to match
            these reference molecules to the molecules appearing in the
            topology. If bond orders are specified in the topology, these will
            be used in matching as well.
        positions
            Positions for the atoms in the new topology.

        Returns
        -------
        topology
            An OpenFF Topology object, constructed from the molecules in
            ``unique_molecules``, with the same atom order as the input topology.

        Raises
        ------
        MissingUniqueMoleculesError
            If ``unique_molecules`` is ``None``
        DuplicateUniqueMoleculeError
            If the same connectivity graph is represented by two different
            molecules in ``unique_molecules``
        ValueError
            If a chemically impossible molecule is detected in the topology
        """
        return cls.from_openmm(
            mdtraj_topology.to_openmm(),
            unique_molecules=unique_molecules,
            positions=positions,
        )

    # Avoid removing this method, even though it is private and would not be difficult for most
    # users to replace. Also avoid making it public as round-trips with MDTraj are likely
    # to not preserve necessary information.
    @requires_package("mdtraj")
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

    def get_bond_between(self, i: Union[int, "Atom"], j: Union[int, "Atom"]) -> "Bond":
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
        from openff.toolkit.topology import Atom

        if (type(i) is int) and (type(j) is int):
            atomi = self.atom(i)
            atomj = self.atom(j)
        elif (type(i) is Atom) and (type(j) is Atom):
            atomi = i
            atomj = j
        else:
            raise ValueError(
                "Invalid input passed to is_bonded(). Expected ints or `Atom`s, "
                "got {} and {}".format(type(i), type(j))
            )

        for bond in atomi.bonds:
            for atom in bond.atoms:
                if atom == atomi:
                    continue
                if atom == atomj:
                    return bond

        raise NotBondedError("No bond between atom {} and {}".format(i, j))

    def is_bonded(self, i: Union[int, "Atom"], j: Union[int, "Atom"]) -> bool:
        """Returns True if the two atoms are bonded

        Parameters
        ----------
        i, j : int or Atom
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

    def atom(self, atom_topology_index: int) -> "Atom":
        """
        Get the Atom at a given Topology atom index.

        Parameters
        ----------
        atom_topology_index : int
             The index of the Atom in this Topology

        Returns
        -------
        An openff.toolkit.topology.Atom
        """
        if not isinstance(atom_topology_index, int):
            raise ValueError(
                "Argument passed to `Topology.atom` must be an int. "
                f"Given argument of type {type(atom_topology_index)}."
            )

        if not (0 <= atom_topology_index < self.n_atoms):
            raise AtomNotInTopologyError(
                f"No atom with index {atom_topology_index} exists in this topology, "
                f"which contains {self.n_atoms} atoms."
            )

        this_molecule_start_index = 0
        next_molecule_start_index = 0
        for molecule in self.molecules:
            next_molecule_start_index += molecule.n_atoms
            if next_molecule_start_index > atom_topology_index:
                atom_molecule_index = atom_topology_index - this_molecule_start_index
                # NOTE: the index here should still be in the topology index order, NOT the reference molecule's
                return molecule.atom(atom_molecule_index)
            this_molecule_start_index += molecule.n_atoms

        raise AtomNotInTopologyError(
            f"No atom with index {atom_topology_index} exists in this topology, "
            f"which contains {self.n_atoms} atoms."
        )

        # Potentially more computationally efficient lookup ( O(largest_molecule_natoms)? )
        # start_index_2_top_mol is an ordered dict of [starting_atom_index] --> [topology_molecule]
        # search_range = range(atom_topology_index - largest_molecule_natoms, atom_topology_index)
        # search_index = atom_topology_index
        # Only efficient if start_index_2_top_mol.keys() is a set (constant time lookups)
        # while not(search_index in start_index_2_top_mol.keys()):
        #     search_index -= 1
        # topology_molecule = start_index_2_top_mol(search_index)
        # atom_molecule_index = atom_topology_index - search_index
        # return topology_molecule.atom(atom_molecule_index)

    def bond(self, bond_topology_index: int) -> "Bond":  # type: ignore[return]
        """
        Get the Bond at a given Topology bond index.

        Parameters
        ----------
        bond_topology_index : int
             The index of the Bond in this Topology

        Returns
        -------
        An openff.toolkit.topology.Bond
        """
        if not isinstance(bond_topology_index, int):
            raise ValueError(
                "Argument passed to `Topology.bond` must be an int. "
                f"Given argument of type {type(bond_topology_index)}."
            )

        if not (0 <= bond_topology_index < self.n_bonds):
            raise BondNotInTopologyError(
                f"No bond with index {bond_topology_index} exists in this topology, "
                f"which contains {self.n_bonds} bonds."
            )

        this_molecule_start_index = 0
        next_molecule_start_index = 0
        for molecule in self._molecules:
            next_molecule_start_index += molecule.n_bonds
            if next_molecule_start_index > bond_topology_index:
                bond_molecule_index = bond_topology_index - this_molecule_start_index
                return molecule.bond(bond_molecule_index)
            this_molecule_start_index += molecule.n_bonds

    def add_molecule(self, molecule: MoleculeLike) -> int:
        """Add a copy of the molecule to the topology"""
        idx = self._add_molecule_keep_cache(molecule)
        self._invalidate_cached_properties()
        return idx

    def _add_molecule_keep_cache(self, molecule: MoleculeLike) -> int:
        self._molecules.append(deepcopy(molecule))
        return len(self._molecules)

    def add_constraint(self, iatom, jatom, distance=True):
        """
        Mark a pair of atoms as constrained.

        Constraints between atoms that are not bonded (e.g., rigid waters) are permissible.

        Parameters
        ----------
        iatom, jatom : Atom
            Atoms to mark as constrained
            These atoms may be bonded or not in the Topology
        distance : unit-wrapped float, optional, default=True
            Constraint distance
            ``True`` if distance has yet to be determined
            ``False`` if constraint is to be removed

        """
        # Check that constraint hasn't already been specified.
        if (iatom, jatom) in self._constrained_atom_pairs:
            existing_distance = self._constrained_atom_pairs[(iatom, jatom)]
            if isinstance(existing_distance, Quantity) and distance is True:
                raise ConstraintExsistsError(
                    f"Atoms ({iatom},{jatom}) already constrained with distance {existing_distance} "
                    "but attempting to override with unspecified distance"
                )
            if (existing_distance is True) and (distance is True):
                raise ConstraintExsistsError(
                    f"Atoms ({iatom},{jatom}) already constrained with unspecified distance "
                    "but attempting to override with unspecified distance"
                )
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
        distance : unit-wrapped float or bool
            True if constrained but constraints have not yet been applied
            Distance if constraint has already been added to Topology

        """
        if (iatom, jatom) in self._constrained_atom_pairs:
            return self._constrained_atom_pairs[(iatom, jatom)]
        else:
            return False

    @requires_package("nglview")
    def visualize(self, ensure_correct_connectivity: bool = False) -> "NGLWidget":
        """
        Visualize with NGLView.

        Requires all molecules in this topology have positions.

        NGLView is a 3D molecular visualization library for use in Jupyter
        notebooks. Note that for performance reasons, by default the
        visualized connectivity is inferred from positions and may not reflect
        the connectivity in the ``Topology``.

        Parameters
        ==========

        ensure_correct_connectivity: bool, default=False
            If ``True``, the visualization will be guaranteed to reflect the
            connectivity in the ``Topology``. Note that this will severely
            degrade performance, especially for topologies with many atoms.

        Examples
        ========

        Visualize a complex PDB file

        >>> from openff.toolkit import Topology
        >>> from openff.toolkit.utils.utils import get_data_file_path
        >>> pdb_filename = get_data_file_path("systems/test_systems/T4_lysozyme_water_ions.pdb")
        >>> topology = Topology.from_pdb(pdb_filename)
        >>> topology.visualize()  # doctest: +SKIP

        """
        import nglview

        from openff.toolkit.utils._viz import TopologyNGLViewStructure

        if ensure_correct_connectivity:
            raise ValueError(
                "`ensure_correct_connectivity` not (yet) implemented "
                "(requires passing multi-molecule SDF files to NGLview)"
            )

        if self.get_positions() is None:
            raise MissingConformersError(
                "All molecules in this topology must have positions for it to be visualized in a widget."
            )

        widget = nglview.NGLWidget(
            TopologyNGLViewStructure(
                topology=self,
                ext="pdb",
            ),
            representations=[
                dict(type="unitcell", params=dict()),
            ],
        )

        widget.add_representation("line", sele="water")
        widget.add_representation("spacefill", sele="ion")
        widget.add_representation("cartoon", sele="protein")
        widget.add_representation(
            "licorice",
            sele="not water and not ion and not protein",
            radius=0.25,
            multipleBond=bool(ensure_correct_connectivity),
        )

        return widget

    def hierarchy_iterator(
        self,
        iter_name: str,
    ) -> Iterator[HierarchyElement]:
        """
        Iterate over all molecules with the given hierarchy scheme.

        Get an iterator over hierarchy elements from all of the molecules in
        this topology that provide the appropriately named iterator. This
        iterator will yield hierarchy elements sorted first by the order that
        molecules are listed in the Topology, and second by the specific
        sorting of hierarchy elements defined in each molecule. Molecules
        without the named iterator are not included.

        Parameters
        ----------
        iter_name
            The iterator name associated with the HierarchyScheme to retrieve
            (for example 'residues' or 'chains')

        Returns
        -------
        iterator of :class:`HierarchyElement`

        See also
        --------
        HierarchyScheme, HierarchyElement, Molecule.hierarchy_schemes
        """
        for molecule in self._molecules:
            if hasattr(molecule, iter_name):
                for item in getattr(molecule, iter_name):
                    yield item

    def __getattr__(self, name: str) -> list["HierarchyElement"]:
        """If a requested attribute is not found, check the hierarchy schemes"""
        # Avoid attempting to process dunder methods as hierarchy scheme iterator names
        if name.startswith("__"):
            raise AttributeError

        try:
            return [
                element
                for molecule in self.molecules
                for element in molecule.hierarchy_schemes[name].hierarchy_elements
            ]
        except (KeyError, AttributeError) as error:
            raise AttributeError(
                f"'{self.__class__.__name__}' object has no attribute '{name}'. If looking for a "
                "`HierarchyScheme` iterator, not all molecules in this topology have an interator "
                f"name {name} defined."
            ) from error
