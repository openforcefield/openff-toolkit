import pytest
from openff.toolkit.typing.engines.smirnoff.parameters import BondHandler
from openff.units import unit

from openff.interchange._tests import _BaseTest
from openff.interchange.components.potentials import (
    Collection,
    Potential,
    WrappedPotential,
)


class TestWrappedPotential(_BaseTest):
    def test_interpolated_potentials(self):
        """Test the construction of and .parameters getter of WrappedPotential"""

        bt = BondHandler.BondType(
            smirks="[#6X4:1]~[#8X2:2]",
            id="bbo1",
            k_bondorder1="100.0 * kilocalories_per_mole/angstrom**2",
            k_bondorder2="200.0 * kilocalories_per_mole/angstrom**2",
            length_bondorder1="1.4 * angstrom",
            length_bondorder2="1.3 * angstrom",
        )

        pot1 = Potential(
            parameters={"k": bt.k_bondorder[1], "length": bt.length_bondorder[1]},
        )
        pot2 = Potential(
            parameters={"k": bt.k_bondorder[2], "length": bt.length_bondorder[2]},
        )

        interp_pot = WrappedPotential(data={pot1: 0.2, pot2: 0.8})
        assert interp_pot.parameters == {
            "k": 180 * unit.Unit("kilocalorie / angstrom ** 2 / mole"),
            "length": 1.32 * unit.angstrom,
        }

        # Ensure a single Potential object can be wrapped with similar behavior
        simple = WrappedPotential(data=pot2)
        assert simple.parameters == pot2.parameters


class TestCollectionSubclassing(_BaseTest):
    def test_dummy_collection(self):
        handler = Collection(
            type="foo",
            expression="m*x+b",
        )
        assert handler.type == "foo"
        assert handler.expression == "m*x+b"


class TestPotentialSerialization(_BaseTest):
    @pytest.fixture()
    def dummy_potential(self):
        return Potential(
            parameters={
                "a": unit.Quantity(1.0, unit.kilocalorie / unit.mole),
                "b": unit.Quantity(2.0, unit.angstrom),
            },
        )

    def test_json_roundtrip(self, dummy_potential):
        potential = Potential.parse_raw(dummy_potential.json())

        assert potential.parameters == dummy_potential.parameters
