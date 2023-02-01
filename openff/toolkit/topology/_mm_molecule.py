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
from typing import TYPE_CHECKING, Dict, List, NoReturn, Optional, Union

from openff.units import unit
from openff.units.elements import MASSES, SYMBOLS

from openff.toolkit.topology.molecule import (
    AtomMetadataDict,
    Molecule,
    _atom_nums_to_hill_formula,
)
from openff.toolkit.utils.exceptions import UnsupportedMoleculeConversionError
from openff.toolkit.utils.utils import deserialize_numpy, serialize_numpy

if TYPE_CHECKING:
    import networkx as nx
    from openff.units.unit import Quantity


class _SimpleMolecule:
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.hierarchy_schemes = dict()
        self.conformers = None

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
            raise Exception(
                "Invalid inputs to molecule._add_bond. Expected ints or Atoms. "
                "Received {} (type {}) and {} (type {}) ".format(
                    atom1, type(atom1), atom2, type(atom2)
                )
            )
        bond = _SimpleBond(atom1_atom, atom2_atom, **kwargs)
        self.bonds.append(bond)

    def add_conformer(self, conformer):
        if self.conformers is None:
            self.conformers = list()
        self.conformers.append(conformer)

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_bonds(self) -> int:
        return len(self.bonds)

    @property
    def n_conformers(self) -> int:
        if self.conformers is None:
            return 0
        return len(self.conformers)

    def atom(self, index):
        return self.atoms[index]

    def atom_index(self, atom) -> int:
        return self.atoms.index(atom)

    def bond(self, index):
        return self.bonds[index]

    def get_bond_between(self, atom1_index, atom2_index):
        atom1 = self.atom(atom1_index)
        atom2 = self.atom(atom2_index)
        for bond in atom1.bonds:
            for atom in bond.atoms:
                if atom is atom2:
                    return bond

    @property
    def angles(self):
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
    def propers(self):
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
                            pass  # yield (atom4, atom3, atom2, atom1)

    @property
    def impropers(self):
        return {}

    @property
    def smirnoff_impropers(self):
        return {}

    @property
    def amber_impropers(self):
        return {}

    @property
    def hill_formula(self) -> str:
        """
        Return the Hill formula of this molecule.

        If any atoms are non-elemental, defined by having an atomic number of 0,
        return INVALID instead.
        """
        return self.to_hill_formula()

    def to_hill_formula(self) -> str:
        atom_nums: List[int] = [atom.atomic_number for atom in self.atoms]

        if min(atom_nums) < 0:
            return "INVALID"

        return _atom_nums_to_hill_formula(atom_nums)

    def to_networkx(self) -> "nx.Graph":
        # TODO: Custom attribtues should probably be attached to the nodes (and possibly also
        #       the edges?). See for more:
        #       https://github.com/openforcefield/openff-toolkit/pull/1179#discussion_r808549385
        import networkx as nx

        graph = nx.Graph()

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

        for index, node_data in subgraph.nodes(data=True):
            molecule.add_atom(atomic_number=node_data["atomic_number"])

        for topology_index1, topology_index2, edge_data in subgraph.edges(data=True):
            molecule.add_bond(topology_index1 - offset, topology_index2 - offset)

        return molecule

    def to_dict(self) -> Dict:
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

        atom_list = list()
        for atom in self.atoms:
            atom_list.append(atom.to_dict())
        molecule_dict["atoms"] = atom_list

        bond_list = list()
        for bond in self.bonds:
            bond_list.append(bond.to_dict())
        molecule_dict["bonds"] = bond_list

        if self.conformers is None:
            molecule_dict["conformers"] = None
        else:
            molecule_dict["conformers"] = []
            molecule_dict[
                "conformers_unit"
            ] = "angstrom"  # Have this defined as a class variable?
            for conf in self.conformers:
                conf_unitless = conf.m_as(unit.angstrom)
                conf_serialized, conf_shape = serialize_numpy((conf_unitless))
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

        bond_dicts = molecule_dict.pop("bonds")

        for bond_dict in bond_dicts:
            atom1_index = bond_dict["atom1_index"]
            atom2_index = bond_dict["atom2_index"]
            molecule.add_bond(
                atom1=molecule.atom(atom1_index), atom2=molecule.atom(atom2_index)
            )

        conformers = molecule_dict.pop("conformers")
        if conformers is None:
            molecule.conformers = None
        else:
            conformers_unit = molecule_dict.pop("conformers_unit")
            molecule.conformers = list()
            for ser_conf in conformers:
                conformers_shape = (molecule.n_atoms, 3)
                conformer_unitless = deserialize_numpy(ser_conf, conformers_shape)
                conformer = unit.Quantity(conformer_unitless, conformers_unit)
                molecule.conformers.append(conformer)

        hier_scheme_dicts = molecule_dict.pop("hierarchy_schemes")
        for iter_name, hierarchy_scheme_dict in hier_scheme_dicts.items():
            new_hier_scheme = cls.add_hierarchy_scheme(
                hierarchy_scheme_dict["uniqueness_criteria"],
                iter_name,
            )
            for element_dict in hierarchy_scheme_dict["hierarchy_elements"]:
                new_hier_scheme.add_hierarchy_element(
                    element_dict["identifier"], element_dict["particle_indices"]
                )
            molecule._expose_hierarchy_scheme(iter_name)

        for key, val in molecule_dict:
            setattr(molecule, key, val)

        return molecule

    @classmethod
    def from_molecule(cls, molecule: Molecule):
        """Generate an MM molecule from an OpenFF Molecule."""
        mm_molecule = cls()
        for atom in molecule.atoms:
            mm_molecule.add_atom(
                atomic_number=atom.atomic_number,
                meatadata=atom.metadata,
            )

        for bond in molecule.bonds:
            mm_molecule.add_bond(
                atom1=mm_molecule.atom(bond.atom1_index),
                atom2=mm_molecule.atom(bond.atom2_index),
            )

        mm_molecule.conformers = molecule.conformers

        return mm_molecule

    def to_molecule(self) -> NoReturn:
        raise UnsupportedMoleculeConversionError(
            "The information content of a _SimpleMolecule is insufficient for creating "
            "an OpenFF Molecule with sufficiently specified chemistry."
        )

    def _is_isomorphic_with(self, other) -> bool:
        """
        Untrustworthy check for pseudo-isomorphism.

        This currently checks that the two molecules have
        * The same number of atoms
        * The same number of bonds
        * Atoms with matching atomic numbers, in order

        This currently does NOT checks that the two molecules have
        * Topologically identical bond graphs
        """

        if self.n_atoms != other.n_atoms:
            return False

        if self.n_bonds != other.n_bonds:
            return False

        for this_atom, other_atom in zip(self.atoms, other.atoms):
            if this_atom.atomic_number != other_atom.atomic_number:
                return False

        return True

    def generate_unique_atom_names(self):
        """Generate unique atom names. See `Molecule.generate_unique_atom_names`."""
        from collections import defaultdict

        element_counts = defaultdict(int)
        for atom in self.atoms:
            symbol = SYMBOLS[atom.atomic_number]
            element_counts[symbol] += 1

            atom.name = symbol + str(element_counts[symbol]) + "x"


class _SimpleAtom:
    def __init__(self, atomic_number: int, molecule=None, metadata=None, **kwargs):
        if metadata is None:
            self.metadata = AtomMetadataDict()
        else:
            self.metadata = AtomMetadataDict(metadata)
        self._atomic_number = atomic_number
        self._molecule = molecule
        self._bonds: List[Optional[_SimpleBond]] = list()
        for key, val in kwargs.items():
            setattr(self, key, val)

    @property
    def atomic_number(self) -> int:
        return self._atomic_number

    @atomic_number.setter
    def atomic_number(self, value):
        if not isinstance(value, int):
            raise ValueError("atomic_number must be an integer")
        if value < 0:
            raise ValueError(
                "atomic_number must be non-negative. An atomic number "
                "of 0 is acceptable."
            )
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

    @property
    def molecule_atom_index(self) -> int:
        return self.molecule.atoms.index(self)

    def to_dict(self) -> Dict[str, Union[Dict, str, int]]:
        atom_dict: Dict[str, Union[Dict, str, int]] = dict()
        atom_dict["metadata"] = dict(self.metadata)
        atom_dict["atomic_number"] = self._atomic_number

        keys_to_skip = ["metadata", "molecule", "bonds"]

        for attr_name, attr_val in self.__dict__.items():
            if attr_name.startswith("_"):
                continue
            if attr_name in keys_to_skip:
                continue
            atom_dict[attr_name] = attr_val
        return atom_dict

    @classmethod
    def from_dict(cls, atom_dict: Dict):
        atom = cls(atomic_number=atom_dict["atomic_number"])
        # TODO: Metadata
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
    def atoms(self) -> List[_SimpleAtom]:
        return [self.atom1, self.atom2]

    @property
    def atom1_index(self) -> int:
        return self.atom1.molecule_atom_index

    @property
    def atom2_index(self) -> int:
        return self.atom2.molecule_atom_index

    def to_dict(self) -> Dict:
        bond_dict = dict()
        bond_dict["atom1_index"] = self.atom1.molecule_atom_index
        bond_dict["atom2_index"] = self.atom2.molecule_atom_index

        return bond_dict
