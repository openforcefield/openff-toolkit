import pathlib
import re

import numpy
import pytest
from openff.utilities import has_package, skip_if_missing

from openff.toolkit import Molecule, unit
from openff.toolkit._tests.create_molecules import (
    create_acetaldehyde,
    create_cis_1_2_dichloroethene,
    create_cyclohexane,
    create_ethanol,
    create_reversed_ethanol,
)
from openff.toolkit._tests.utils import requires_openeye
from openff.toolkit.utils import GLOBAL_TOOLKIT_REGISTRY
from openff.toolkit.utils.exceptions import (
    ToolkitUnavailableException,
)
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper

_DEFAULT_MODEL = "openff-gnn-am1bcc-0.1.0-rc.3.pt"
from openff.nagl_models._dynamic_fetch import BadFileSuffixError

try:
    from openff.nagl_models import list_available_nagl_models
except ModuleNotFoundError:

    def list_available_nagl_models():
        return list()


@skip_if_missing("openff.nagl")
class TestNAGLToolkitWrapper:
    def test_version(self):
        from openff.nagl import __version__ as parsed_version

        assert parsed_version == NAGLToolkitWrapper()._toolkit_version

    def test_nagl_in_global_toolkit_registry(self):
        assert "NAGL" in GLOBAL_TOOLKIT_REGISTRY.__repr__()


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
    @pytest.mark.parametrize(
        "nagl_model",
        [pathlib.Path(file).name for file in list_available_nagl_models()],
    )
    def test_assign_partial_charges_basic(self, molecule_function, nagl_model):
        molecule = molecule_function()

        molecule.assign_partial_charges(
            partial_charge_method="am1bccelf10",
            toolkit_registry=OpenEyeToolkitWrapper(),
        )

        openeye_charges = molecule.partial_charges

        molecule = molecule_function()

        molecule.assign_partial_charges(
            partial_charge_method=nagl_model,
            toolkit_registry=NAGLToolkitWrapper(),
        )

        assert molecule.partial_charges is not None

        nagl_charges = molecule.partial_charges.m_as(unit.elementary_charge)
        assert nagl_charges.dtype == float

        numpy.testing.assert_allclose(
            openeye_charges.m_as(unit.elementary_charge),
            nagl_charges,
            atol=0.07,
        )

    def test_atom_order_dependence(self):
        """Regression test against atom order dependence."""
        forward = create_ethanol()
        reverse = create_reversed_ethanol()

        for molecule in [forward, reverse]:
            molecule.assign_partial_charges(
                partial_charge_method=_DEFAULT_MODEL,
                toolkit_registry=NAGLToolkitWrapper(),
            )

        numpy.testing.assert_allclose(
            forward.partial_charges,
            reverse.partial_charges[::-1],
            atol=1e-7,
        )

    def test_conformer_argument(self):
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=2)

        with pytest.warns(
            UserWarning,
            match="optional argument.*use_conformers",
        ):
            molecule.assign_partial_charges(
                partial_charge_method=_DEFAULT_MODEL,
                toolkit_registry=NAGLToolkitWrapper(),
                use_conformers=[molecule.conformers[-1]],
            )

        with pytest.warns(
            UserWarning,
            match="optional argument.*strict_n_conformers",
        ):
            molecule.assign_partial_charges(
                partial_charge_method=_DEFAULT_MODEL,
                toolkit_registry=NAGLToolkitWrapper(),
                strict_n_conformers=1,
            )

    def test_unsupported_charge_method(self):
        with pytest.raises(
            BadFileSuffixError,
            match="NAGLToolkitWrapper does not recognize file path extension on filename='hartree_fock'",
        ):
            create_ethanol().assign_partial_charges(
                partial_charge_method="hartree_fock",
                toolkit_registry=NAGLToolkitWrapper(),
            )

    def test_unsupported_molecule_element(self):
        si = Molecule.from_smiles("[Si+4]")
        with pytest.raises(ValueError, match="Molecule contains forbidden element 14"):
            si.assign_partial_charges(
                partial_charge_method=_DEFAULT_MODEL,
                toolkit_registry=NAGLToolkitWrapper(),
            )

    def test_unsupported_molecule_bond(self):
        mol = Molecule.from_smiles("CCC=[Cl+]")
        err = re.escape("Molecule contains forbidden SMARTS pattern [#17:1]#,:,=[*:2]")
        with pytest.raises(ValueError, match=err):
            mol.assign_partial_charges(
                partial_charge_method=_DEFAULT_MODEL,
                toolkit_registry=NAGLToolkitWrapper(),
            )


@pytest.mark.skipif(
    has_package("openff.nagl"),
    reason="Test requires that OpenFF NAGL is not installed",
)
class TestNoNAGLToolkitWrapper:
    def test_nagl_toolkit_wrapper_unavailable(self):
        assert not NAGLToolkitWrapper.is_available()

    def test_init_unavailable(self):
        with pytest.raises(
            ToolkitUnavailableException,
            match="OpenFF NAGL is not available",
        ):
            NAGLToolkitWrapper()
