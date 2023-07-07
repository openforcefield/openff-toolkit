import pytest
from openff.toolkit import ForceField, Molecule, Topology
from openff.toolkit.topology._mm_molecule import _SimpleMolecule

from openff.interchange._tests import _BaseTest
from openff.interchange.components.toolkit import (
    _check_electrostatics_handlers,
    _combine_topologies,
    _get_14_pairs,
    _get_num_h_bonds,
    _simple_topology_from_openmm,
)


@pytest.fixture()
def simple_methane():
    return _SimpleMolecule.from_molecule(Molecule.from_smiles("C"))


@pytest.fixture()
def simple_water():
    return _SimpleMolecule.from_molecule(Molecule.from_smiles("O"))


def test_simple_topology_uniqueness(simple_methane, simple_water):
    topology = Topology.from_molecules(
        [
            simple_methane,
            simple_water,
            simple_methane,
            simple_methane,
            simple_water,
        ],
    )
    assert len(topology.identical_molecule_groups) == 2


class TestToolkitUtils(_BaseTest):
    @pytest.mark.parametrize(
        ("smiles", "num_pairs"),
        [
            ("C#C", 1),
            ("CCO", 12),
            ("C1=CC=CC=C1", 24),
            ("C=1=C=C1", 0),
            ("C=1=C=C=C1", 0),
            ("C=1(Cl)-C(Cl)=C1", 1),
            ("C=1=C(Cl)C(=C=1)Cl", 5),
        ],
    )
    def test_get_14_pairs(self, smiles, num_pairs):
        mol = Molecule.from_smiles(smiles)
        assert len([*_get_14_pairs(mol)]) == num_pairs
        assert len([*_get_14_pairs(mol.to_topology())]) == num_pairs

    def test_check_electrostatics_handlers(self, tip3p_missing_electrostatics_xml):
        # https://github.com/openforcefield/openff-toolkit/blob/0.10.2/openff/toolkit/data/test_forcefields/tip3p.offxml
        tip3p_missing_electrostatics = ForceField(tip3p_missing_electrostatics_xml)

        assert _check_electrostatics_handlers(tip3p_missing_electrostatics)

        tip3p_missing_electrostatics.deregister_parameter_handler("LibraryCharges")

        assert not _check_electrostatics_handlers(tip3p_missing_electrostatics)

    @pytest.mark.parametrize(
        ("smiles", "num_h_bonds"),
        [("C", 4), ("C#C", 2), ("O", 2)],
    )
    def test_get_num_h_bonds(self, smiles, num_h_bonds):
        topology = Molecule.from_smiles(smiles).to_topology()
        assert _get_num_h_bonds(topology) == num_h_bonds, smiles

    def test_combine_topologies(self):
        ethanol = Molecule.from_smiles("CCO")
        ethanol.name = "ETH"
        ethanol_topology = ethanol.to_topology()

        water = Molecule.from_smiles("O")
        water.name = "WAT"
        water_topology = water.to_topology()

        combined = _combine_topologies(ethanol_topology, water_topology)

        for attr in (
            "atoms",
            "bonds",
        ):
            attr = "n_" + attr
            assert getattr(combined, attr) == getattr(ethanol_topology, attr) + getattr(
                water_topology,
                attr,
            )

    def test_simple_topology_from_openmm(self):
        simple_topology = _simple_topology_from_openmm(
            Topology.from_molecules(
                [
                    Molecule.from_smiles("O"),
                    Molecule.from_smiles("CCO"),
                ],
            ).to_openmm(),
        )

        assert all(
            isinstance(molecule, _SimpleMolecule)
            for molecule in simple_topology.molecules
        )

        assert sorted(molecule.n_atoms for molecule in simple_topology.molecules) == [
            3,
            9,
        ]
