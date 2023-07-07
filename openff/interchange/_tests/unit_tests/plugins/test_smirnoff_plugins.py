import pytest
from nonbonded_plugins.nonbonded import (
    BuckinghamHandler,
    DoubleExponentialHandler,
    SMIRNOFFBuckinghamCollection,
    SMIRNOFFDoubleExponentialCollection,
)
from openff.toolkit import ForceField, Molecule
from openff.toolkit.typing.engines.smirnoff.plugins import load_handler_plugins

from openff.interchange import Interchange
from openff.interchange.plugins import load_smirnoff_plugins


def test_load_handler_plugins():
    # This does not test Interchange per se but this behavior is necessary for the other tests to function
    available_plugins = load_handler_plugins()

    assert BuckinghamHandler in available_plugins
    assert DoubleExponentialHandler in available_plugins


def test_load_smirnoff_plugins():
    available_plugins = load_smirnoff_plugins()

    assert SMIRNOFFBuckinghamCollection in available_plugins
    assert SMIRNOFFDoubleExponentialCollection in available_plugins


class TestDoubleExponential:
    pytest.importorskip("deforcefields")

    @pytest.fixture()
    def de_force_field(self) -> ForceField:
        force_field = ForceField("de-force-1.0.0.offxml", load_plugins=True)

        # An early version of this force field includes a blank vdW handler
        # https://github.com/jthorton/de-forcefields/blob/aaaa4c721967310ec8e9e3176b97de858273b2d5/deforcefields/offxml/de-force-1.0.0.offxml#L325-L327
        force_field.deregister_parameter_handler("vdW")

        # This force field also includes a 4-site water model, remove that for now
        force_field.deregister_parameter_handler("VirtualSites")
        force_field.deregister_parameter_handler("LibraryCharges")
        force_field["DoubleExponential"].parameters.pop(-1)
        force_field["DoubleExponential"].parameters.pop(-1)
        force_field["Constraints"].parameters.pop(-1)
        force_field["Constraints"].parameters.pop(-1)

        return force_field

    def test_loadable(self):
        ForceField("de-force-1.0.0.offxml", load_plugins=True)

    def test_create_interchange(self, de_force_field):
        Interchange.from_smirnoff(
            de_force_field,
            Molecule.from_smiles("CO").to_topology(),
        )
