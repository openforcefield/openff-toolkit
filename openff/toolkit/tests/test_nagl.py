import numpy
from openff.utilities.testing import skip_if_missing

from openff.toolkit.tests.create_molecules import create_ethanol
from openff.toolkit.tests.utils import requires_openeye
from openff.toolkit.utils._nagl_wrapper import _NAGLToolkitWrapper
from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper


@skip_if_missing("openff.nagl")
class TestNAGLToolkitWrapper:
    def test_version(self):
        assert "0.1.0" in _NAGLToolkitWrapper()._toolkit_version

    @requires_openeye
    def test_charges_sanity(self):
        molecule = create_ethanol()

        molecule.assign_partial_charges(
            partial_charge_method="am1bccelf10",
            toolkit_registry=OpenEyeToolkitWrapper(),
        )

        openeye_charges = molecule.partial_charges

        molecule = create_ethanol()

        molecule.assign_partial_charges(
            partial_charge_method="_nagl_am1bccelf10",
            toolkit_registry=_NAGLToolkitWrapper(),
        )

        nagl_charges = molecule.partial_charges

        assert numpy.allclose(openeye_charges, nagl_charges)
