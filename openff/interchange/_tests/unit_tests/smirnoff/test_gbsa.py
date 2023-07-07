import pytest
from openff.toolkit import Molecule
from openff.toolkit.typing.engines.smirnoff.parameters import GBSAHandler
from openff.units import unit

from openff.interchange.constants import kcal_mol_a2
from openff.interchange.models import TopologyKey
from openff.interchange.smirnoff._gbsa import SMIRNOFFGBSACollection


class TestGBSACollection:
    @pytest.fixture()
    def gbsa_handler(self) -> GBSAHandler:
        handler = GBSAHandler(
            version="0.3",
            gb_model="OBC1",
            solvent_dielectric="78.5",
            solute_dielectric="1.0",
            sa_model="ACE",
            surface_area_penalty=5.4 * kcal_mol_a2,
            solvent_radius=1.4 * unit.angstrom,
        )

        for parameter in [
            {"smirks": "[#1:1]", "radius": "0.12 * nanometer", "scale": 0.85},
            {"smirks": "[#1:1]~[#6]", "radius": "0.13 * nanometer", "scale": 0.85},
            {"smirks": "[#1:1]~[#8]", "radius": "0.08 * nanometer", "scale": 0.85},
            {"smirks": "[#1:1]~[#16]", "radius": "0.08 * nanometer", "scale": 0.85},
            {"smirks": "[#6:1]", "radius": "0.22 * nanometer", "scale": 0.72},
            {"smirks": "[#7:1]", "radius": "0.155 * nanometer", "scale": 0.79},
            {"smirks": "[#8:1]", "radius": "0.15 * nanometer", "scale": 0.85},
            {"smirks": "[#9:1]", "radius": "0.15 * nanometer", "scale": 0.88},
            {"smirks": "[#14:1]", "radius": "0.21 * nanometer", "scale": 0.8},
            {"smirks": "[#15:1]", "radius": "0.185 * nanometer", "scale": 0.86},
            {"smirks": "[#16:1]", "radius": "0.18 * nanometer", "scale": 0.96},
            {"smirks": "[#17:1]", "radius": "0.17 * nanometer", "scale": 0.8},
        ]:
            handler.add_parameter(parameter)

        return handler

    def test_create_basic(self, gbsa_handler):
        collection = SMIRNOFFGBSACollection.create(
            gbsa_handler,
            Molecule.from_mapped_smiles(
                "[H:4][C:1]([H:5])([H:6])[C:2]([H:7])([H:8])[O:3][H:9]",
            ).to_topology(),
        )

        assert collection.gb_model == "OBC1"
        assert collection.solvent_dielectric == 78.5
        assert collection.solute_dielectric == 1.0
        assert collection.sa_model == "ACE"
        assert collection.surface_area_penalty == 5.4 * kcal_mol_a2
        assert collection.solvent_radius == 1.4 * unit.angstrom

        assert len(collection.key_map) == 9
        assert len(collection.potentials) == 4

        oxygen_parameters = collection.potentials[
            collection.key_map[TopologyKey(atom_indices=(2,))]
        ].parameters

        assert oxygen_parameters["radius"] == unit.Quantity(0.15, unit.nanometer)
        assert oxygen_parameters["scale"] == unit.Quantity(0.85, unit.dimensionless)
