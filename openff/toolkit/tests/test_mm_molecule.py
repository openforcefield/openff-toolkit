import numpy as np
import pytest

from openff.toolkit.tests.create_molecules import create_ethanol
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


class TestIsomorphism:
    def test_are_isomorphic(self, water, methane, methanol):
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
