import pytest
from openff.toolkit import ForceField, Molecule
from openff.toolkit.typing.engines.smirnoff.parameters import (
    ChargeIncrementModelHandler,
)

from openff.interchange._tests import get_test_file_path


@pytest.fixture()
def cb8_host() -> Molecule:
    return Molecule.from_file(get_test_file_path("CB8.sdf"))


@pytest.fixture()
def no_charges() -> ForceField:
    sage = ForceField("openff_unconstrained-2.0.0.offxml")
    sage.deregister_parameter_handler("ToolkitAM1BCC")
    sage.register_parameter_handler(
        ChargeIncrementModelHandler(version=0.3, partial_charge_method="formal_charge"),
    )

    return sage
