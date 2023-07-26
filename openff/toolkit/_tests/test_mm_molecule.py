import numpy as np
import pytest

from openff.toolkit._tests.create_molecules import create_ethanol
from openff.toolkit.topology._mm_molecule import _SimpleMolecule
from openff.toolkit.topology.molecule import Molecule


@pytest.fixture()
def water():
    water = _SimpleMolecule()
    water.add_atom(atomic_number=8)
    water.add_atom(atomic_number=1)
    water.add_atom(atomic_number=1)
    water.add_bond(0, 1)
    water.add_bond(0, 2)

    return water


@pytest.fixture()
def methane():
    methane = _SimpleMolecule()
    methane.add_atom(atomic_number=6)
    methane.add_atom(atomic_number=1)
    methane.add_atom(atomic_number=1)
    methane.add_atom(atomic_number=1)
    methane.add_atom(atomic_number=1)
    methane.add_bond(0, 1)
    methane.add_bond(0, 2)
    methane.add_bond(0, 3)
    methane.add_bond(0, 4)

    return methane


@pytest.fixture()
def methanol():
    methanol = _SimpleMolecule()
    methanol.add_atom(atomic_number=8)
    methanol.add_atom(atomic_number=6)
    methanol.add_atom(atomic_number=1)
    methanol.add_atom(atomic_number=1)
    methanol.add_atom(atomic_number=1)
    methanol.add_bond(0, 1)
    methanol.add_bond(1, 2)
    methanol.add_bond(1, 3)
    methanol.add_bond(1, 4)

    return methanol


@pytest.fixture()
def molecule_with_zero_atom():
    molecule = _SimpleMolecule()
    molecule.add_atom(atomic_number=6)
    molecule.add_atom(atomic_number=0)
    molecule.add_bond(0, 1)

    return molecule


@pytest.fixture()
def molecule_with_bogus_atom():
    molecule = _SimpleMolecule()
    molecule.add_atom(atomic_number=6)
    molecule.add_atom(atomic_number=-1)
    molecule.add_bond(0, 1)

    return molecule


class TestMMMolecule:
    def test_create_water(self):
        water = _SimpleMolecule()
        water.add_atom(atomic_number=8)
        water.add_atom(atomic_number=1)
        water.add_atom(atomic_number=1)
        water.add_bond(0, 1)
        water.add_bond(0, 2)
        assert water.n_atoms == 3
        assert water.n_bonds == 2
        assert water.n_conformers == 0

    def test_bond_getters(self, water):
        assert water.get_bond_between(0, 1) is not water.get_bond_between(0, 2)
        assert water.bond(0) is not water.bond(1)

        assert water.get_bond_between(0, 1) is water.bond(0)
        assert water.get_bond_between(0, 2) is water.bond(1)

    def test_hill_formula(
        self,
        water,
        molecule_with_zero_atom,
        molecule_with_bogus_atom,
    ):
        assert water.hill_formula == "H2O"
        assert molecule_with_zero_atom.hill_formula == "CX"
        assert molecule_with_bogus_atom.hill_formula == "INVALID"

    def test_to_networkx(self, water):
        graph = water.to_networkx()

        assert graph.number_of_nodes() == water.n_atoms
        assert graph.number_of_edges() == water.n_bonds

    def test_add_conformer(self, water):
        water_molecule = Molecule.from_smiles("O")
        water_molecule.generate_conformers(n_conformers=1)

        assert water.n_conformers == 0
        water.add_conformer(water_molecule.conformers[0])
        assert water.n_conformers == 1

        assert np.allclose(water.conformers[0], water_molecule.conformers[0])

    def test_dict_roundtrip(self, water):
        roundtrip = _SimpleMolecule.from_dict(water.to_dict())

        assert roundtrip.n_atoms == water.n_atoms
        assert roundtrip.n_bonds == water.n_bonds

        for atom_index in range(roundtrip.n_atoms):
            assert (
                roundtrip.atom(atom_index).atomic_number
                == water.atom(atom_index).atomic_number  # noqa
            )

    def test_dict_roundtrip_conformers(self, water):
        water_molecule = Molecule.from_smiles("O")
        water_molecule.generate_conformers(n_conformers=1)

        water.add_conformer(water_molecule.conformers[0])

        roundtrip = _SimpleMolecule.from_dict(water.to_dict())

        assert water.n_conformers == roundtrip.n_conformers
        assert np.allclose(water.conformers[0], roundtrip.conformers[0])

    def test_from_molecule(self):
        converted = _SimpleMolecule.from_molecule(Molecule.from_smiles("O"))

        assert converted.n_atoms == 3
        assert converted.n_bonds == 2

        expected_atomic_numbers = [8, 1, 1]
        for atom_index in range(converted.n_atoms):
            found = converted.atom(atom_index).atomic_number
            assert found == expected_atomic_numbers[atom_index]


class TestImpropers:
    @pytest.mark.parametrize(
        ("smiles", "n_impropers", "n_pruned"),
        [
            ("C", 24, 0),
            ("CC", 48, 0),
            ("N", 6, 6),
        ],
    )
    def test_pruned_impropers(self, smiles, n_impropers, n_pruned):
        """See equivalent test in TestMolecule."""
        molecule = _SimpleMolecule.from_molecule(
            Molecule.from_smiles(smiles),
        )

        assert molecule.n_impropers == n_impropers
        assert len(list(molecule.smirnoff_impropers)) == n_pruned
        assert len(list(molecule.amber_impropers)) == n_pruned

        amber_impropers = {*molecule.amber_impropers}

        for smirnoff_imp in molecule.smirnoff_impropers:
            assert (
                smirnoff_imp[1],
                smirnoff_imp[0],
                smirnoff_imp[2],
                smirnoff_imp[3],
            ) in amber_impropers


class TestIsomorphism:
    @pytest.fixture()
    def n_propanol(self):
        return _SimpleMolecule.from_molecule(
            Molecule.from_mapped_smiles(
                "[H:5][C:1]([H:6])([H:7])[C:2]([H:8])([H:9])[C:3]([H:10])([H:11])[O:4][H:12]"
            )
        )

    @pytest.fixture()
    def iso_propanol(self):
        return _SimpleMolecule.from_molecule(
            Molecule.from_mapped_smiles(
                "[H:5][C:1]([H:6])([H:7])[C:2]([H:8])([C:3]([H:10])([H:11])[H:12])[O:4][H:9]"
            )
        )

    @pytest.fixture()
    def o_dichlorobezene(self):
        return _SimpleMolecule.from_molecule(
            Molecule.from_mapped_smiles(
                "[C:1]=1([Cl:7])[C:2]([H:9])=[C:3]([H:10])[C:4]([H:11])=[C:5]([H:12])[C:6]([Cl:8])=1"
            )
        )

    @pytest.fixture()
    def m_dichlorobezene(self):
        return _SimpleMolecule.from_molecule(
            Molecule.from_mapped_smiles(
                "[C:1]([Cl:7])(=[C:2]1([H:9]))[C:3]([H:10])=[C:4]([H:11])[C:5]([H:12])=[C:6]([Cl:8])1"
            )
        )

    @pytest.mark.parametrize("as_graphs", [True, False])
    def test_are_isomorphic(self, water, methane, methanol, as_graphs):
        if as_graphs:
            water = water.to_networkx()
            methane = methane.to_networkx()
            methanol = methanol.to_networkx()

        assert _SimpleMolecule.are_isomorphic(water, water)[0]
        assert _SimpleMolecule.are_isomorphic(methane, methane)[0]
        assert _SimpleMolecule.are_isomorphic(methanol, methanol)[0]
        assert not _SimpleMolecule.are_isomorphic(water, methane)[0]
        assert not _SimpleMolecule.are_isomorphic(water, methanol)[0]
        assert not _SimpleMolecule.are_isomorphic(methane, methanol)[0]

    def test_is_isomorphic_with(self, water, methane, methanol):
        assert water.is_isomorphic_with(water)
        assert methane.is_isomorphic_with(methane)
        assert methanol.is_isomorphic_with(methanol)
        assert not water.is_isomorphic_with(methane)
        assert not water.is_isomorphic_with(methanol)
        assert not methane.is_isomorphic_with(methanol)

    def test_short_circuit_heterogeneous_input(self):
        assert not _SimpleMolecule.are_isomorphic(
            create_ethanol(),
            _SimpleMolecule.from_molecule(create_ethanol()),
        )[0]

        assert not _SimpleMolecule.are_isomorphic(
            _SimpleMolecule.from_molecule(create_ethanol()),
            create_ethanol(),
        )[0]

    def test_graph_and_molecule_inputs(self, methanol):
        graph = methanol.to_networkx()

        assert _SimpleMolecule.are_isomorphic(methanol, graph)[0]
        assert _SimpleMolecule.are_isomorphic(graph, methanol)[0]

    def test_propanol_isopropanol_not_isomorphic(self, n_propanol, iso_propanol):
        assert not _SimpleMolecule.are_isomorphic(n_propanol, iso_propanol)[0]

    def test_o_m_dichlorobenzene_not_isomorphic(
        self, o_dichlorobezene, m_dichlorobezene
    ):
        assert not _SimpleMolecule.are_isomorphic(o_dichlorobezene, m_dichlorobezene)[0]
