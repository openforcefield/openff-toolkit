import re

import numpy
import pytest
from openff.units import unit
from openff.utilities.testing import skip_if_missing

from openff.toolkit import Molecule
from openff.toolkit._tests.create_molecules import (
    create_acetaldehyde,
    create_cis_1_2_dichloroethene,
    create_cyclohexane,
    create_ethanol,
    create_reversed_ethanol,
)
from openff.toolkit._tests.utils import requires_openeye
from openff.toolkit.utils._nagl_wrapper import _NAGLToolkitWrapper
from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper


@skip_if_missing("openff.nagl")
class TestNAGLToolkitWrapper:
    def test_version(self):
        from openff.nagl import __version__ as parsed_version

        assert parsed_version == _NAGLToolkitWrapper()._toolkit_version

    @requires_openeye
    @pytest.mark.parametrize(
        "molecule_function",
        [
            create_ethanol,
            create_cis_1_2_dichloroethene,
            create_acetaldehyde,
            create_cyclohexane,
        ],
    )
    def test_assign_partial_charges_basic(self, molecule_function):
        molecule = molecule_function()

        molecule.assign_partial_charges(
            partial_charge_method="am1bccelf10",
            toolkit_registry=OpenEyeToolkitWrapper(),
        )

        openeye_charges = molecule.partial_charges

        molecule = molecule_function()

        molecule.assign_partial_charges(
            partial_charge_method="_nagl_am1bccelf10",
            toolkit_registry=_NAGLToolkitWrapper(),
        )

        assert molecule.partial_charges is not None

        nagl_charges = molecule.partial_charges.m_as(unit.elementary_charge)
        assert nagl_charges.dtype == float

        numpy.testing.assert_allclose(
            openeye_charges.m_as(unit.elementary_charge), nagl_charges, atol=0.07
        )

    def test_atom_order_dependence(self):
        """Regression test against atom order dependence."""
        forward = create_ethanol()
        reverse = create_reversed_ethanol()

        for molecule in [forward, reverse]:
            molecule.assign_partial_charges(
                partial_charge_method="_nagl_am1bccelf10",
                toolkit_registry=_NAGLToolkitWrapper(),
            )
        numpy.testing.assert_allclose(
            forward.partial_charges,
            reverse.partial_charges[::-1],
            atol=1e-7,
        )

    def test_unsupported_charge_method(self):
        with pytest.raises(ValueError, match="Charge model hartree_fock not supported"):
            create_ethanol().assign_partial_charges(
                partial_charge_method="hartree_fock",
                toolkit_registry=_NAGLToolkitWrapper(),
            )

    def test_unsupported_molecule_element(self):
        si = Molecule.from_smiles("[Si+4]")
        with pytest.raises(ValueError, match="Molecule contains forbidden element 14"):
            si.assign_partial_charges(
                partial_charge_method="_nagl_am1bccelf10",
                toolkit_registry=_NAGLToolkitWrapper(),
            )

    def test_unsupported_molecule_bond(self):
        mol = Molecule.from_smiles("C=[Cl+1]")
        err = re.escape("Molecule contains forbidden SMARTS pattern [#17:1]#,:,=[*:2]")
        with pytest.raises(ValueError, match=err):
            mol.assign_partial_charges(
                partial_charge_method="_nagl_am1bccelf10",
                toolkit_registry=_NAGLToolkitWrapper(),
            )
