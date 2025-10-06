"""
TypedMolecule TODOs
* Testing
* Detemine whether atomic numbers are required, or whether a TypedMol could have a mix of
  TypedAtom and TypedParticles, where the latter doesn't need atomic nubmers
* Determine whether typedmols should have isomorphism operations
* Should be able to go to hill formula (?), networkx, openmm
* Topology serialization will have trouble here - Won't know whether it's trying to
  deserialize a Molecule or a TypedMolecule.

"""

import functools
from collections.abc import Generator, Iterable
from typing import TYPE_CHECKING, NoReturn, Union

from openff.units.elements import MASSES, SYMBOLS

from openff.toolkit import unit
from openff.toolkit.topology.molecule import (
    AtomMetadataDict,
    HierarchyElement,
    HierarchyScheme,
    Molecule,
    _atom_nums_to_hill_formula,
    _has_unique_atom_names,
)
from openff.toolkit.utils.exceptions import UnsupportedMoleculeConversionError
from openff.toolkit.utils.utils import deserialize_numpy, serialize_numpy

if TYPE_CHECKING:
    import networkx as nx

    from openff.toolkit import Quantity, Topology
    from openff.toolkit.topology.molecule import FrozenMolecule


class _SimpleMolecule:
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self._conformers = None
        self._hierarchy_schemes = dict()

    def add_atom(self, atomic_number: int, **kwargs):
        atom = _SimpleAtom(atomic_number, self, **kwargs)
        self.atoms.append(atom)

    def add_bond(self, atom1, atom2, **kwargs):
        if isinstance(atom1, int) and isinstance(atom2, int):
            atom1_atom = self.atoms[atom1]
            atom2_atom = self.atoms[atom2]
        elif isinstance(atom1, _SimpleAtom) and isinstance(atom2, _SimpleAtom):
            atom1_atom = atom1
            atom2_atom = atom2
        else:
            raise ValueError(
                "Invalid inputs to molecule._add_bond. Expected ints or Atoms. "
                f"Received {atom1} (type {type(atom1)}) and {atom2} (type {type(atom2)}) "
            )
        bond = _SimpleBond(atom1_atom, atom2_atom, **kwargs)
        self.bonds.append(bond)

    @property
    def conformers(self) -> list["Quantity"] | None:
        return self._conformers

    def add_conformer(self, conformer):
        if self._conformers is None:
            self._conformers = list()
        self._conformers.append(conformer)

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_bonds(self) -> int:
        return len(self.bonds)

    @property
    def n_conformers(self) -> int:
        return 0 if self._conformers is None else len(self._conformers)

    def atom(self, index: int) -> "_SimpleAtom":
        return self.atoms[index]

    def atom_index(self, atom) -> int:
        return self.atoms.index(atom)

    def bond(self, index: int) -> "_SimpleBond":
        return self.bonds[index]

    def get_bond_between(self, atom1_index, atom2_index):
        atom1 = self.atom(atom1_index)
        atom2 = self.atom(atom2_index)
        for bond in atom1.bonds:
            for atom in bond.atoms:
                if atom is atom2:
                    return bond

    @property
    def angles(
        self,
    ) -> Generator[
        tuple[
            "_SimpleAtom",
            "_SimpleAtom",
            "_SimpleAtom",
        ],
        None,
        None,
    ]:
        for atom1 in self.atoms:
            for atom2 in atom1.bonded_atoms:
                for atom3 in atom2.bonded_atoms:
                    if atom1 is atom3:
                        continue
                    if atom1.molecule_atom_index < atom3.molecule_atom_index:
                        yield (atom1, atom2, atom3)
                    else:
                        # Do not return i.e. (2, 1, 0), only (0, 1, 2)
                        pass

    @property
    def propers(
        self,
    ) -> Generator[
        tuple[
            "_SimpleAtom",
            "_SimpleAtom",
            "_SimpleAtom",
            "_SimpleAtom",
        ],
        None,
        None,
    ]:
        for atom1 in self.atoms:
            for atom2 in atom1.bonded_atoms:
                for atom3 in atom2.bonded_atoms:
                    if atom1 is atom3:
                        continue
                    for atom4 in atom3.bonded_atoms:
                        if atom4 in (atom1, atom2):
                            continue

                        if atom1.molecule_atom_index < atom4.molecule_atom_index:
                            yield (atom1, atom2, atom3, atom4)
                        else:
                            # Do no duplicate
                            pass

    @property
    def impropers(
        self,
    ) -> Generator[
        tuple[
            "_SimpleAtom",
            "_SimpleAtom",
            "_SimpleAtom",
            "_SimpleAtom",
        ],
        None,
        None,
    ]:
        for atom1 in self.atoms:
            for atom2 in atom1.bonded_atoms:
                for atom3 in atom2.bonded_atoms:
                    if atom1 is atom3:
                        continue
                    for atom3i in atom2.bonded_atoms:
                        if atom3i == atom3:
                            continue
                        if atom3i == atom1:
                            continue

                        yield (atom1, atom2, atom3, atom3i)

    @property
    def smirnoff_impropers(
        self,
    ) -> Generator[
        tuple[
            "_SimpleAtom",
            "_SimpleAtom",
            "_SimpleAtom",
            "_SimpleAtom",
        ],
        None,
        None,
    ]:
        for improper in self.impropers:
            if len(list(improper[1].bonded_atoms)) == 3:
                yield improper

    @property
    def amber_impropers(
        self,
    ) -> Generator[
        tuple[
            "_SimpleAtom",
            "_SimpleAtom",
            "_SimpleAtom",
            "_SimpleAtom",
        ],
        None,
        None,
    ]:
        for improper in self.smirnoff_impropers:
            yield (improper[1], improper[0], improper[2], improper[3])

    @property
    def n_angles(self) -> int:
        return len(list(self.angles))

    @property
    def n_propers(self) -> int:
        return len(list(self.propers))

    @property
    def n_impropers(self) -> int:
        return len(list(self.impropers))

    @property
    def hill_formula(self) -> str:
        """
        Return the Hill formula of this molecule.

        If any atoms are non-elemental, defined by having an atomic number of 0,
        return INVALID instead.
        """
        return self.to_hill_formula()

    @property
    def hierarchy_schemes(self) -> dict[str, "HierarchyScheme"]:
        return self._hierarchy_schemes

    def to_hill_formula(self) -> str:
        atom_nums: list[int] = [atom.atomic_number for atom in self.atoms]

        if min(atom_nums) < 0:
            return "INVALID"

        return _atom_nums_to_hill_formula(atom_nums)

    def to_networkx(self) -> "nx.Graph":
        # TODO: Custom attribtues should probably be attached to the nodes (and possibly also
        #       the edges?). See for more:
        #       https://github.com/openforcefield/openff-toolkit/pull/1179#discussion_r808549385
        import networkx as nx

        graph: nx.classes.graph.Graph = nx.Graph()

        for atom_index, atom in enumerate(self.atoms):
            graph.add_node(
                atom_index,
                atomic_number=atom.atomic_number,
            )
        for bond in self.bonds:
            graph.add_edge(
                bond.atom1_index,
                bond.atom2_index,
            )

        return graph

    def to_topology(self) -> "Topology":
        from openff.toolkit.topology import Topology

        return Topology.from_molecules([self])

    def nth_degree_neighbors(self, n_degrees):
        from openff.toolkit.topology.molecule import (
            _nth_degree_neighbors_from_graphlike,
        )

        return _nth_degree_neighbors_from_graphlike(graphlike=self, n_degrees=n_degrees)

    @classmethod
    def _from_subgraph(cls, subgraph: "nx.Graph"):
        molecule = cls()

        # The subgraph stores indices that might not start at zero (i.e. topology indices)
        # but we need 0-index indices (like molecule indices) to add bonds
        offset = min(subgraph.nodes())

        for _, node_data in subgraph.nodes(data=True):
            # Hierarchy metadata like residue name needs to be passed in as a separate argument,
            # so we extract those values from the node data
            metadata = dict()
            for key, val in node_data.items():
                if key in [
                    "residue_name",
                    "residue_number",
                    "insertion_code",
                    "chain_id",
                ]:
                    metadata[key] = val
            # Then we remove the metadata items that we took from the node data
            for key in metadata.keys():
                del node_data[key]

            molecule.add_atom(metadata=metadata, **node_data)

        for topology_index1, topology_index2, _edge_data in subgraph.edges(data=True):
            molecule.add_bond(topology_index1 - offset, topology_index2 - offset)

        return molecule

    # TODO: Implement me - shouldn't need to be too different than Molecule.add_hierarchy_scheme
    #       since the extra chemical information is unrelated to hierarchy/metadata
    def add_hierarchy_scheme(
        self,
        uniqueness_criteria: Iterable[str],
        iterator_name: str,
    ) -> NoReturn:
        raise NotImplementedError()

    def to_dict(self) -> dict:
        molecule_dict = dict()
        special_serialization_logic = [
            "atoms",
            "bonds",
            "conformers",
            "hierarchy_schemes",
            "molecule",
        ]
        for attr_name, attr_val in self.__dict__.items():
            if attr_name.startswith("_"):
                continue
            if attr_name in special_serialization_logic:
                continue
            molecule_dict[attr_name] = attr_val

        molecule_dict["atoms"] = [atom.to_dict() for atom in self.atoms]

        molecule_dict["bonds"] = [bond.to_dict() for bond in self.bonds]

        if self._conformers is None:
            molecule_dict["conformers"] = None
        else:
            molecule_dict["conformers"] = []
            molecule_dict["conformers_unit"] = "angstrom"  # Have this defined as a class variable?
            for conf in self._conformers:
                conf_unitless = conf.m_as(unit.angstrom)
                conf_serialized, _ = serialize_numpy(conf_unitless)
                molecule_dict["conformers"].append(conf_serialized)

        molecule_dict["hierarchy_schemes"] = dict()
        for iter_name, hier_scheme in self.hierarchy_schemes.items():
            molecule_dict["hierarchy_schemes"][iter_name] = hier_scheme.to_dict()

        return molecule_dict

    @classmethod
    def from_dict(cls, molecule_dict):
        molecule = cls()

        atom_dicts = molecule_dict.pop("atoms")
        for atom_dict in atom_dicts:
            molecule.atoms.append(_SimpleAtom.from_dict(atom_dict))
            molecule.atoms[-1]._molecule = molecule

        bond_dicts = molecule_dict.pop("bonds")

        for bond_dict in bond_dicts:
            atom1_index = bond_dict["atom1_index"]
            atom2_index = bond_dict["atom2_index"]
            molecule.add_bond(atom1=molecule.atom(atom1_index), atom2=molecule.atom(atom2_index))

        conformers = molecule_dict.pop("conformers")
        if conformers is None:
            molecule._conformers = None
        else:
            conformers_unit = molecule_dict.pop("conformers_unit")
            molecule._conformers = list()
            for ser_conf in conformers:
                conformers_shape = (molecule.n_atoms, 3)
                conformer_unitless = deserialize_numpy(ser_conf, conformers_shape)
                conformer = unit.Quantity(conformer_unitless, conformers_unit)
                molecule._conformers.append(conformer)

        hier_scheme_dicts = molecule_dict.pop("hierarchy_schemes")
        for iterator_name, hierarchy_scheme_dict in hier_scheme_dicts.items():
            molecule._hierarchy_schemes[iterator_name] = HierarchyScheme(
                parent=molecule,
                uniqueness_criteria=tuple(hierarchy_scheme_dict["uniqueness_criteria"]),
                iterator_name=iterator_name,
            )

            for element_dict in hierarchy_scheme_dict["hierarchy_elements"]:
                molecule._hierarchy_schemes[iterator_name].add_hierarchy_element(
                    identifier=element_dict["identifier"],
                    atom_indices=element_dict["atom_indices"],
                )

        for key, val in molecule_dict.items():
            setattr(molecule, key, val)

        return molecule

    @classmethod
    def from_molecule(cls, molecule: Molecule):
        """Generate an MM molecule from an OpenFF Molecule."""
        mm_molecule = cls()
        for atom in molecule.atoms:
            mm_molecule.add_atom(
                atomic_number=atom.atomic_number,
                metadata=atom.metadata,
            )

        for bond in molecule.bonds:
            mm_molecule.add_bond(
                atom1=mm_molecule.atom(bond.atom1_index),
                atom2=mm_molecule.atom(bond.atom2_index),
            )

        if molecule.conformers is not None:
            for conformer in molecule.conformers:
                mm_molecule.add_conformer(conformer)

        for name, hierarchy_scheme in molecule.hierarchy_schemes.items():
            assert name == hierarchy_scheme.iterator_name
            mm_molecule.add_hierarchy_scheme(
                uniqueness_criteria=hierarchy_scheme.uniqueness_criteria,
                iterator_name=hierarchy_scheme.iterator_name,
            )

        return mm_molecule

    def to_molecule(self) -> NoReturn:
        raise UnsupportedMoleculeConversionError(
            "The information content of a _SimpleMolecule is insufficient for creating "
            "an OpenFF Molecule with sufficiently specified chemistry."
        )

    def is_isomorphic_with(self, other: Union["FrozenMolecule", "_SimpleMolecule", "nx.Graph"], **kwargs) -> bool:
        """
        Check for pseudo-isomorphism.

        This currently checks that the two molecules have
        * The same number of atoms
        * The same number of bonds
        * Atoms with matching atomic numbers, in order

        This currently does NOT checks that the two molecules have
        * Topologically identical bond graphs
        """
        return _SimpleMolecule.are_isomorphic(
            self,
            other,
            return_atom_map=False,
        )[0]

    @staticmethod
    def are_isomorphic(
        mol1: Union["FrozenMolecule", "_SimpleMolecule", "nx.Graph"],
        mol2: Union["FrozenMolecule", "_SimpleMolecule", "nx.Graph"],
        return_atom_map: bool = False,
    ) -> tuple[bool, dict[int, int] | None]:
        import networkx

        _cls = _SimpleMolecule

        if isinstance(mol1, networkx.Graph) and isinstance(mol2, networkx.Graph):
            graph1 = _SimpleMolecule._from_subgraph(mol1).to_networkx()
            graph2 = _SimpleMolecule._from_subgraph(mol2).to_networkx()

        elif isinstance(mol1, networkx.Graph):
            assert isinstance(mol2, _cls)
            graph1 = _SimpleMolecule._from_subgraph(mol1).to_networkx()
            graph2 = mol2.to_networkx()

        elif isinstance(mol2, networkx.Graph):
            assert isinstance(mol1, _cls)
            graph1 = mol1.to_networkx()
            graph2 = _SimpleMolecule._from_subgraph(mol2).to_networkx()

        else:
            # static methods (by definition) know nothing about their class,
            # so the class to compare to must be hard-coded here
            if not (isinstance(mol1, _cls) and isinstance(mol2, _cls)):
                return False, None

            graph1 = mol1.to_networkx()
            graph2 = mol2.to_networkx()

        def node_match_func(node1, node2):
            return node1["atomic_number"] == node2["atomic_number"]

        edge_match_func = None

        matcher = networkx.algorithms.isomorphism.GraphMatcher(
            graph1, graph2, node_match=node_match_func, edge_match=edge_match_func
        )

        if matcher.is_isomorphic():
            if return_atom_map:
                topology_atom_map = matcher.mapping

                return True, {key: topology_atom_map[key] for key in sorted(topology_atom_map)}

            else:
                return True, None
        else:
            return False, None

    def generate_unique_atom_names(self):
        """Generate unique atom names. See `Molecule.generate_unique_atom_names`."""
        from collections import defaultdict

        element_counts: defaultdict[str, int] = defaultdict(int)
        for atom in self.atoms:
            symbol = SYMBOLS[atom.atomic_number]
            element_counts[symbol] += 1

            atom.name = symbol + str(element_counts[symbol]) + "x"

    @property
    def has_unique_atom_names(self) -> bool:
        """``True`` if the molecule has unique atom names, ``False`` otherwise."""
        return _has_unique_atom_names(self)

    def __getattr__(self, name: str) -> list["HierarchyElement"]:
        """If a requested attribute is not found, check the hierarchy schemes"""
        try:
            return self.__dict__["_hierarchy_schemes"][name].hierarchy_elements
        except KeyError:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute {name!r}")

    def __deepcopy__(self, memo):
        return self.__class__.from_dict(self.to_dict())


_SimpleMolecule.add_hierarchy_scheme = Molecule.add_hierarchy_scheme  # type: ignore[assignment]
_SimpleMolecule.update_hierarchy_schemes = Molecule.update_hierarchy_schemes  # type: ignore[attr-defined]


class _SimpleAtom:
    def __init__(
        self,
        atomic_number: int,
        molecule: _SimpleMolecule | None = None,
        metadata: AtomMetadataDict | None = None,
        name: str = "",
        **kwargs,
    ):
        if metadata is None:
            self.metadata = AtomMetadataDict()
        else:
            self.metadata = AtomMetadataDict(metadata)

        # This _could_ be hidden away into _metadata
        self._name = name
        self._atomic_number = atomic_number
        self._molecule = molecule
        self._bonds: list[_SimpleBond] = list()
        for key, val in kwargs.items():
            setattr(self, key, val)

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, value: str):
        self._name = value

    @property
    def atomic_number(self) -> int:
        return self._atomic_number

    @atomic_number.setter
    def atomic_number(self, value):
        if not isinstance(value, int):
            raise ValueError("atomic_number must be an integer")
        if value < 0:
            raise ValueError("atomic_number must be non-negative. An atomic number of 0 is acceptable.")
        self._atomic_number = value

    @property
    def symbol(self) -> str:
        return SYMBOLS[self.atomic_number]

    @property
    def mass(self) -> "Quantity":
        return MASSES[self.atomic_number]

    @property
    def molecule(self):
        """The ``Molecule`` this atom is part of."""
        return self._molecule

    @property
    def bonds(self):
        return self._bonds

    def add_bond(self, bond):
        self._bonds.append(bond)

    @property
    def bonded_atoms(self):
        for bond in self._bonds:
            for atom in [bond.atom1, bond.atom2]:
                if atom is not self:
                    yield atom

    @functools.cached_property
    def molecule_atom_index(self) -> int:
        return self.molecule.atoms.index(self)

    def to_dict(self) -> dict[str, dict | str | int]:
        atom_dict: dict[str, dict | str | int] = dict()
        atom_dict["metadata"] = dict(self.metadata)
        atom_dict["atomic_number"] = self._atomic_number
        atom_dict["name"] = self._name

        keys_to_skip = ["metadata", "molecule", "bonds"]

        for attr_name, attr_val in self.__dict__.items():
            if attr_name.startswith("_"):
                continue
            if attr_name in keys_to_skip:
                continue
            atom_dict[attr_name] = attr_val

        return atom_dict

    @classmethod
    def from_dict(cls, atom_dict: dict):
        atom = cls(**atom_dict)
        return atom


class _SimpleBond:
    def __init__(self, atom1, atom2, **kwargs):
        self.molecule = atom1.molecule
        self.atom1 = atom1
        self.atom2 = atom2

        atom1.add_bond(self)
        atom2.add_bond(self)

        for key, val in kwargs.items():
            setattr(self, key, val)

    @property
    def atoms(self) -> list[_SimpleAtom]:
        return [self.atom1, self.atom2]

    @property
    def atom1_index(self) -> int:
        return self.atom1.molecule_atom_index

    @property
    def atom2_index(self) -> int:
        return self.atom2.molecule_atom_index

    def to_dict(self) -> dict[str, int]:
        return {
            "atom1_index": self.atom1.molecule_atom_index,
            "atom2_index": self.atom2.molecule_atom_index,
        }
