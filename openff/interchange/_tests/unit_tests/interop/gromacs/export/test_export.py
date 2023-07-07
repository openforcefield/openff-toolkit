import pytest
from openff.toolkit import ForceField, Molecule
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest
from openff.interchange.interop.gromacs._import._import import from_files


class TestToGro(_BaseTest):
    def test_residue_names(self, sage):
        """Reproduce issue #642."""
        # This could maybe just test the behavior of _convert?
        ligand = Molecule.from_smiles("CCO")
        ligand.generate_conformers(n_conformers=1)

        for atom in ligand.atoms:
            atom.metadata["residue_name"] = "LIG"

        Interchange.from_smirnoff(
            sage,
            [ligand],
        ).to_gro("should_have_residue_names.gro")

        for line in open("should_have_residue_names.gro").readlines()[2:-2]:
            assert line[5:10] == "LIG  "


class TestSettles(_BaseTest):
    def test_settles_units(self, monkeypatch):
        """Reproduce issue #720."""
        monkeypatch.setenv("INTERCHANGE_EXPERIMENTAL", "1")

        molecule = Molecule.from_mapped_smiles("[H:2][O:1][H:3]")
        molecule.generate_conformers(n_conformers=1)
        molecule.name = "WAT"

        ForceField("openff-2.1.0.offxml").create_interchange(
            molecule.to_topology(),
        ).to_gromacs(
            prefix="settles",
        )

        system = from_files("settles.top", "settles.gro")

        settles = system.molecule_types["WAT"].settles[0]

        assert settles.oxygen_hydrogen_distance.m_as(unit.angstrom) == pytest.approx(
            0.981124525388,
        )
        assert settles.hydrogen_hydrogen_distance.m_as(unit.angstrom) == pytest.approx(
            1.513900654525,
        )
