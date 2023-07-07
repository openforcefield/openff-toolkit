import numpy
import pytest
from openff.toolkit import Molecule, Topology
from openff.units import unit

from openff.interchange import Interchange
from openff.interchange._tests import _BaseTest
from openff.interchange.interop.gromacs._interchange import (
    _convert_topology,
    to_interchange,
)
from openff.interchange.interop.gromacs.models.models import GROMACSSystem
from openff.interchange.smirnoff._gromacs import _convert


class TestToInterchange(_BaseTest):
    @pytest.fixture()
    def simple_interchange(self, sage_unconstrained) -> Interchange:
        topology = Topology()
        for index, (smiles, name) in enumerate(
            zip(
                [
                    "[H:2][C:1]([H:3])([H:4])[H:5]",
                    "[H:3][C:1](=[C:2]([H:5])[H:6])[H:4]",
                    "[H:2][O:1][H:3]",
                ],
                ["MET", "ETH", "WAT"],
            ),
        ):
            molecule = Molecule.from_mapped_smiles(smiles)
            molecule.generate_conformers(n_conformers=1)
            molecule.conformers[0] += numpy.array([index, 0, 0]) * unit.nanometer
            molecule.name = name

            topology.add_molecule(molecule)

        topology.box_vectors = [4, 4, 4] * unit.nanometer

        return sage_unconstrained.create_interchange(topology)

    @pytest.mark.slow()
    def test_convert_basic_system(self, monkeypatch, simple_interchange):
        monkeypatch.setenv("INTERCHANGE_EXPERIMENTAL", "1")

        converted = to_interchange(_convert(simple_interchange))

        assert numpy.allclose(converted.positions, simple_interchange.positions)
        assert numpy.allclose(converted.box, simple_interchange.box)

        for name in simple_interchange.collections:
            # Currently two bond potentials (for rigid water) are lost in the conversion;
            # not clear if this information should be retained.
            # * bond parameters would need to be written out to [ bonds ] even with [ settles ] populated
            # * GromacsMolecule.bonds would need to be populated even if it has a settles
            # * the collection would need to store the bond parameters despite also being constrained

            # Similar situation with angles

            assert name in converted.collections


class TestConvertTopology(_BaseTest):
    @pytest.fixture()
    def simple_system(self, sage_unconstrained) -> GROMACSSystem:
        molecule = Molecule.from_smiles("CCO")
        molecule.generate_conformers(n_conformers=1)
        topology = Topology.from_molecules([molecule])
        topology.box_vectors = [4, 4, 4] * unit.nanometer

        return _convert(sage_unconstrained.create_interchange(topology))

    @pytest.fixture()
    def water_dimer(self, sage_unconstrained):
        water = Molecule.from_mapped_smiles("[H:2][O:1][H:3]")
        water.name = "WAT"
        water.generate_conformers(n_conformers=1)
        topology = Topology.from_molecules([water, water])
        topology.box_vectors = [4, 4, 4] * unit.nanometer

        return _convert(sage_unconstrained.create_interchange(topology))

    def test_convert_basic_system(self, monkeypatch, simple_system):
        monkeypatch.setenv("INTERCHANGE_EXPERIMENTAL", "1")
        converted = _convert_topology(simple_system)

        assert converted.n_molecules == 1
        assert converted.n_atoms == 9
        assert converted.n_bonds == 8

    def test_water_with_settles_has_bonds_in_topology(self, water_dimer):
        assert "WAT" in water_dimer.molecule_types
        assert len(water_dimer.molecule_types["WAT"].settles) == 1

        converted = _convert_topology(water_dimer)

        assert converted.molecule(0).name == "WAT"
        assert converted.molecule(1).name == "WAT"
        assert converted.n_molecules == 2
        assert converted.n_atoms == 6
        assert converted.n_bonds == 4

        assert [atom.atomic_number for atom in converted.atoms] == [8, 1, 1, 8, 1, 1]
