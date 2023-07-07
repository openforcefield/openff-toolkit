"""Test SMIRNOFF-GROMACS conversion."""
import pytest
from openff.toolkit import Molecule
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest
from openff.interchange.interop.gromacs.models.models import GROMACSMolecule
from openff.interchange.smirnoff._gromacs import (
    _convert,
    _convert_angles,
    _convert_bonds,
    _convert_settles,
)


class TestConvert(_BaseTest):
    def test_residue_names(self, sage):
        """Reproduce issue #642."""
        ligand = Molecule.from_smiles("CCO")
        ligand.generate_conformers(n_conformers=1)

        for atom in ligand.atoms:
            atom.metadata["residue_name"] = "LIG"

        system = _convert(
            Interchange.from_smirnoff(
                sage,
                [ligand],
            ),
        )

        for molecule_type in system.molecule_types.values():
            for atom in molecule_type.atoms:
                assert atom.residue_name == "LIG"


class TestSettles(_BaseTest):
    @pytest.fixture()
    def tip3p_interchange(self, tip3p, water):
        return tip3p.create_interchange(water.to_topology())

    def test_catch_other_water_ordering(self, tip3p):
        molecule = Molecule.from_mapped_smiles("[H:1][O:2][H:3]")
        interchange = tip3p.create_interchange(molecule.to_topology())

        with pytest.raises(Exception, match="OHH"):
            _convert_settles(
                GROMACSMolecule(name="foo"),
                interchange.topology.molecule(0),
                interchange,
            )

    def test_convert_settles(self, tip3p_interchange):
        molecule = GROMACSMolecule(name="foo")

        _convert_settles(
            molecule,
            tip3p_interchange.topology.molecule(0),
            tip3p_interchange,
        )

        assert len(molecule.settles) == 1

        settle = molecule.settles[0]

        assert settle.first_atom == 1
        assert settle.hydrogen_hydrogen_distance.m_as(unit.angstrom) == pytest.approx(
            1.5139006545247014,
        )
        assert settle.oxygen_hydrogen_distance.m_as(unit.angstrom) == pytest.approx(
            0.9572,
        )

        assert molecule.exclusions[0].first_atom == 1
        assert molecule.exclusions[0].other_atoms == [2, 3]
        assert molecule.exclusions[1].first_atom == 2
        assert molecule.exclusions[1].other_atoms == [3]

    def test_convert_no_settles_unconstrained_water(self, tip3p_interchange):
        tip3p_interchange.collections["Constraints"].key_map = dict()

        molecule = GROMACSMolecule(name="foo")

        _convert_settles(
            molecule,
            tip3p_interchange.topology.molecule(0),
            tip3p_interchange,
        )

        assert len(molecule.settles) == 0

    def test_convert_no_settles_no_constraints(self, tip3p_interchange):
        tip3p_interchange.collections.pop("Constraints")

        molecule = GROMACSMolecule(name="foo")

        _convert_settles(
            molecule,
            tip3p_interchange.topology.molecule(0),
            tip3p_interchange,
        )

        assert len(molecule.settles) == 0

    def test_no_bonds_or_angles_if_settle(self, tip3p_interchange):
        molecule = GROMACSMolecule(name="foo")

        for function in [_convert_settles, _convert_bonds, _convert_angles]:
            function(
                molecule,
                tip3p_interchange.topology.molecule(0),
                tip3p_interchange,
            )

        assert len(molecule.settles) == 1
        assert len(molecule.angles) == 0
        assert len(molecule.bonds) == 0
