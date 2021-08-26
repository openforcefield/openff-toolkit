"""
Tests for chemiformatics toolkit wrapper helper classes
"""
from typing import Dict

import numpy as np
import pytest
from simtk import unit

from openff.toolkit.tests.create_molecules import create_acetate, create_ethanol
from openff.toolkit.tests.utils import (
    requires_ambertools,
    requires_openeye,
    requires_rdkit,
)
from openff.toolkit.topology.molecule import Molecule
from openff.toolkit.utils import get_data_file_path
from openff.toolkit.utils.exceptions import (
    ChargeMethodUnavailableError,
    IncorrectNumConformersError,
    IncorrectNumConformersWarning,
    InvalidToolkitError,
    ToolkitUnavailableException,
)
from openff.toolkit.utils.toolkits import (
    GLOBAL_TOOLKIT_REGISTRY,
    AmberToolsToolkitWrapper,
    BuiltInToolkitWrapper,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
    ToolkitRegistry,
    ToolkitWrapper,
)

# =============================================================================================
# FIXTURES
# =============================================================================================


def get_mini_drug_bank(toolkit_class, xfail_mols=None):
    """Read the mini drug bank sdf file with the toolkit and return the molecules"""

    # This is a work around a weird error where even though the test is skipped due to a missing toolkit
    #  we still try and read the file with the toolkit
    if toolkit_class.is_available():
        toolkit = toolkit_class()
        molecules = Molecule.from_file(
            get_data_file_path("molecules/MiniDrugBank.sdf"),
            "sdf",
            toolkit_registry=toolkit,
            allow_undefined_stereo=True,
        )
    else:
        molecules = []

    if xfail_mols is None:
        return molecules

    for i, mol in enumerate(molecules):
        if mol.name in xfail_mols:
            marker = pytest.mark.xfail(reason=xfail_mols[mol.name])
            molecules[i] = pytest.param(mol, marks=marker)

    return molecules


openeye_inchi_stereochemistry_lost = [
    "DrugBank_2799",
    "DrugBank_5414",
    "DrugBank_5415",
    "DrugBank_5418",
    "DrugBank_2955",
    "DrugBank_2987",
    "DrugBank_5555",
    "DrugBank_472",
    "DrugBank_5737",
    "DrugBank_3332",
    "DrugBank_3461",
    "DrugBank_794",
    "DrugBank_3502",
    "DrugBank_6026",
    "DrugBank_3622",
    "DrugBank_977",
    "DrugBank_3693",
    "DrugBank_3726",
    "DrugBank_3739",
    "DrugBank_6222",
    "DrugBank_6232",
    "DrugBank_3844",
    "DrugBank_6295",
    "DrugBank_6304",
    "DrugBank_6305",
    "DrugBank_3930",
    "DrugBank_6329",
    "DrugBank_6353",
    "DrugBank_6355",
    "DrugBank_6401",
    "DrugBank_4161",
    "DrugBank_4162",
    "DrugBank_6509",
    "DrugBank_6531",
    "DrugBank_1570",
    "DrugBank_4249",
    "DrugBank_1634",
    "DrugBank_1659",
    "DrugBank_6647",
    "DrugBank_1700",
    "DrugBank_1721",
    "DrugBank_1742",
    "DrugBank_1802",
    "DrugBank_6775",
    "DrugBank_1849",
    "DrugBank_1864",
    "DrugBank_6875",
    "DrugBank_1897",
    "DrugBank_4593",
    "DrugBank_1962",
    "DrugBank_4662",
    "DrugBank_7049",
    "DrugBank_4702",
    "DrugBank_2095",
    "DrugBank_4778",
    "DrugBank_2141",
    "DrugBank_2148",
    "DrugBank_2178",
    "DrugBank_4865",
    "DrugBank_2208",
    "DrugBank_2210",
    "DrugBank_2276",
    "DrugBank_4959",
    "DrugBank_4964",
    "DrugBank_5043",
    "DrugBank_2429",
    "DrugBank_5076",
    "DrugBank_2465",
    "DrugBank_2519",
    "DrugBank_2538",
    "DrugBank_5158",
    "DrugBank_5176",
    "DrugBank_2592",
]

openeye_inchi_isomorphic_fails = ["DrugBank_1661", "DrugBank_4346", "DrugBank_2467"]

rdkit_inchi_stereochemistry_lost = [
    "DrugBank_5414",
    "DrugBank_2955",
    "DrugBank_5737",
    "DrugBank_3332",
    "DrugBank_3461",
    "DrugBank_6026",
    "DrugBank_3622",
    "DrugBank_3726",
    "DrugBank_6222",
    "DrugBank_3844",
    "DrugBank_6304",
    "DrugBank_6305",
    "DrugBank_6329",
    "DrugBank_6509",
    "DrugBank_6647",
    "DrugBank_1897",
    "DrugBank_4778",
    "DrugBank_2148",
    "DrugBank_2178",
    "DrugBank_2538",
    "DrugBank_2592",
    "DrugBank_4249",
    "DrugBank_5076",
    "DrugBank_5418",
    "DrugBank_3930",
    "DrugBank_1634",
    "DrugBank_1962",
    "DrugBank_5043",
    "DrugBank_2519",
    "DrugBank_7124",
    "DrugBank_6865",
]

rdkit_inchi_roundtrip_mangled = ["DrugBank_2684"]

openeye_iupac_bad_stereo = [
    "DrugBank_977",
    "DrugBank_1634",
    "DrugBank_1700",
    "DrugBank_1962",
    "DrugBank_2148",
    "DrugBank_2178",
    "DrugBank_2186",
    "DrugBank_2208",
    "DrugBank_2519",
    "DrugBank_2538",
    "DrugBank_2592",
    "DrugBank_2651",
    "DrugBank_2987",
    "DrugBank_3332",
    "DrugBank_3502",
    "DrugBank_3622",
    "DrugBank_3726",
    "DrugBank_3844",
    "DrugBank_3930",
    "DrugBank_4161",
    "DrugBank_4162",
    "DrugBank_4778",
    "DrugBank_4593",
    "DrugBank_4959",
    "DrugBank_5043",
    "DrugBank_5076",
    "DrugBank_5176",
    "DrugBank_5418",
    "DrugBank_5737",
    "DrugBank_5902",
    "DrugBank_6295",
    "DrugBank_6304",
    "DrugBank_6305",
    "DrugBank_6329",
    "DrugBank_6355",
    "DrugBank_6401",
    "DrugBank_6509",
    "DrugBank_6531",
    "DrugBank_6647",
    "DrugBank_390",
    "DrugBank_810",
    "DrugBank_4316",
    "DrugBank_4346",
    "DrugBank_7124",
    "DrugBank_2799",
    "DrugBank_4662",
    "DrugBank_4865",
    "DrugBank_2465",
]


@pytest.fixture()
def formic_acid_molecule() -> Molecule:

    formic_acid = Molecule()
    formic_acid.add_atom(8, 0, False)  # O1
    formic_acid.add_atom(6, 0, False)  # C1
    formic_acid.add_atom(8, 0, False)  # O2
    formic_acid.add_atom(1, 0, False)  # H1
    formic_acid.add_atom(1, 0, False)  # H2
    formic_acid.add_bond(0, 1, 2, False)  # O1 - C1
    formic_acid.add_bond(1, 2, 1, False)  # C1 - O2
    formic_acid.add_bond(1, 3, 1, False)  # C1 - H1
    formic_acid.add_bond(2, 4, 1, False)  # O2 - H2

    return formic_acid


@pytest.fixture()
def formic_acid_conformers() -> Dict[str, unit.Quantity]:

    return {
        "cis": np.array(
            [
                [-0.95927322, -0.91789997, 0.36333418],
                [-0.34727824, 0.12828046, 0.22784603],
                [0.82766682, 0.26871252, -0.42284882],
                [-0.67153811, 1.10376000, 0.61921501],
                [1.15035689, -0.58282924, -0.78766006],
            ]
        )
        * unit.angstrom,
        "trans": np.array(
            [
                [-0.95927322, -0.91789997, 0.36333418],
                [-0.34727824, 0.12828046, 0.22784603],
                [0.82766682, 0.26871252, -0.42284882],
                [-0.67153811, 1.10376000, 0.61921501],
                [1.14532626, 1.19679034, -0.41266876],
            ]
        )
        * unit.angstrom,
    }


# =============================================================================================
# TESTS
# =============================================================================================


class TestBuiltInToolkitWrapper:
    """Test the BuiltInToolkitWrapper"""

    @pytest.mark.parametrize("partial_charge_method", ["zeros", "formal_charge"])
    def test_assign_partial_charges_neutral(self, partial_charge_method):
        """Test BuiltInToolkitWrapper assign_partial_charges()"""

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[BuiltInToolkitWrapper])
        molecule = create_ethanol()
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
            # The normalize_partial_charges kwarg won't cause a measurable effect here, this just
            # tests that the kwarg remains accepted.
            normalize_partial_charges=False,
        )
        charge_sum = sum(molecule.partial_charges, 0.0 * unit.elementary_charge)
        assert 1.0e-10 > abs(charge_sum / unit.elementary_charge)

    @pytest.mark.parametrize("partial_charge_method", ["formal_charge"])
    def test_assign_partial_charges_net_charge(self, partial_charge_method):
        """
        Test BuiltInToolkitWrapper assign_partial_charges(). Only formal_charge is tested, since zeros will not
        sum up to the proper number
        """

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[BuiltInToolkitWrapper])
        molecule = create_acetate()
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
        )
        charge_sum = sum(molecule.partial_charges, 0.0 * unit.elementary_charge)
        assert 1.0e-10 > abs(charge_sum.value_in_unit(unit.elementary_charge) + 1.0)

    def test_assign_partial_charges_bad_charge_method(self):
        """Test BuiltInToolkitWrapper assign_partial_charges() for a nonexistent charge method"""

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[BuiltInToolkitWrapper])
        molecule = create_ethanol()

        # For now, the Molecule API passes raise_exception_types=[] to ToolkitRegistry.call,
        # which loses track of what exception type
        # was thrown inside them, so we just check for a ValueError here
        with pytest.raises(
            ValueError, match="is not supported by the Built-in toolkit"
        ) as excinfo:
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method="NotARealChargeMethod",
            )

        # ToolkitWrappers raise a specific exception class, so we test that here
        with pytest.raises(
            ChargeMethodUnavailableError,
            match="is not supported by the Built-in toolkit",
        ) as excinfo:
            BITKW = BuiltInToolkitWrapper()
            BITKW.assign_partial_charges(
                molecule=molecule, partial_charge_method="NotARealChargeMethod"
            )

    def test_assign_partial_charges_wrong_n_confs(self):
        """
        Test BuiltInToolkitWrapper assign_partial_charges() when requesting to use an incorrect number of
        conformers
        """

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[BuiltInToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=1)
        with pytest.warns(
            IncorrectNumConformersWarning,
            match="has 1 conformers, but charge method 'zeros' expects exactly 0.",
        ):
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method="zeros",
                use_conformers=molecule.conformers,
                strict_n_conformers=False,
            )

        # Specify strict_n_conformers=True, but not use_conformers, so a recommended number of
        # conformers will be generated internally
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method="zeros",
            strict_n_conformers=True,
        )

        # For now, the Molecule API passes raise_exception_types=[] to ToolkitRegistry.call,
        # which loses track of what exception type
        # was thrown inside them, so we just check for a ValueError here
        with pytest.raises(
            ValueError,
            match=f"has 1 conformers, but charge method 'zeros' " f"expects exactly 0.",
        ):
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method="zeros",
                use_conformers=molecule.conformers,
                strict_n_conformers=True,
            )

        # Test calling the ToolkitWrapper _directly_, passing in the incorrect number of
        # confs, and specify strict_n_conformers, which should produce an IncorrectNumConformersError
        with pytest.raises(
            IncorrectNumConformersError,
            match=f"has 1 conformers, but charge method 'zeros' " f"expects exactly 0.",
        ):
            BITKW = BuiltInToolkitWrapper()
            BITKW.assign_partial_charges(
                molecule=molecule,
                partial_charge_method="zeros",
                use_conformers=molecule.conformers,
                strict_n_conformers=True,
            )


class TestToolkitWrapper:
    """Test the ToolkitWrapper class"""

    def test_check_n_conformers(self):
        """Ensure that _check_n_conformers is working properly"""
        tkw = ToolkitWrapper()
        mol = create_ethanol()

        ## Test molecule with no conformers
        # Check with no min or max should pass
        tkw._check_n_conformers(mol, "nocharge")
        # Check with min=1 should warn
        with pytest.warns(
            IncorrectNumConformersWarning,
            match="has 0 conformers, but charge method 'nocharge' expects at least 1",
        ):
            tkw._check_n_conformers(mol, "nocharge", min_confs=1)
        # Check with min=1 and strict_n_conformers should raise an error
        with pytest.raises(
            IncorrectNumConformersError,
            match="has 0 conformers, but charge method 'nocharge' expects at least 1",
        ):
            tkw._check_n_conformers(
                mol, "nocharge", min_confs=1, strict_n_conformers=True
            )
        # Check with min=1, max=1 and strict_n_conformers should raise an error
        with pytest.raises(
            IncorrectNumConformersError,
            match="has 0 conformers, but charge method 'nocharge' expects exactly 1",
        ):
            tkw._check_n_conformers(
                mol, "nocharge", min_confs=1, max_confs=1, strict_n_conformers=True
            )
        # Check with min=1, max=2 and strict_n_conformers should raise an error
        with pytest.raises(
            IncorrectNumConformersError,
            match="has 0 conformers, but charge method 'nocharge' expects between 1 and 2",
        ):
            tkw._check_n_conformers(
                mol, "nocharge", min_confs=1, max_confs=2, strict_n_conformers=True
            )
        # Check with max=1 should pass
        tkw._check_n_conformers(mol, "nocharge", max_confs=1, strict_n_conformers=True)

        ## Test molecule with conformers
        # Add some conformers
        mol.generate_conformers(n_conformers=1)
        for _ in range(9):
            mol.add_conformer(mol.conformers[0])

        # Check with no min or max should pass
        tkw._check_n_conformers(mol, "nocharge")

        ## min_confs checks
        # Check with min=1 should be fine
        tkw._check_n_conformers(mol, "nocharge", min_confs=1)
        # Check with min=10 should be fine
        tkw._check_n_conformers(mol, "nocharge", min_confs=10)
        # Check with min=11 should warn
        with pytest.warns(
            IncorrectNumConformersWarning,
            match="has 10 conformers, but charge method 'nocharge' expects at least 11",
        ):
            tkw._check_n_conformers(mol, "nocharge", min_confs=11)
        # Check with min=11 and strict_n_conformers should raise an error
        with pytest.raises(
            IncorrectNumConformersError,
            match="has 10 conformers, but charge method 'nocharge' expects at least 11",
        ):
            tkw._check_n_conformers(
                mol, "nocharge", min_confs=11, strict_n_conformers=True
            )

        ## max_confs checks
        # Check with max=1 and strict_n_conformers should raise an error
        with pytest.raises(
            IncorrectNumConformersError,
            match="has 10 conformers, but charge method 'nocharge' expects at most 1",
        ):
            tkw._check_n_conformers(
                mol, "nocharge", max_confs=1, strict_n_conformers=True
            )
        # Check with max=10 and strict_n_conformers should be OK
        tkw._check_n_conformers(mol, "nocharge", max_confs=10, strict_n_conformers=True)
        # Check with max=11 and strict_n_conformers should be OK
        tkw._check_n_conformers(mol, "nocharge", max_confs=11, strict_n_conformers=True)

        ## min_confs and max_confs checks
        # Check with max=10 and min=10 and strict_n_conformers should be OK
        tkw._check_n_conformers(
            mol, "nocharge", min_confs=10, max_confs=10, strict_n_conformers=True
        )
        # Check with max=10 and min=9 and strict_n_conformers should be OK
        tkw._check_n_conformers(
            mol, "nocharge", min_confs=9, max_confs=10, strict_n_conformers=True
        )
        # Check with max=11 and min=10 and strict_n_conformers should be OK
        tkw._check_n_conformers(
            mol, "nocharge", min_confs=10, max_confs=11, strict_n_conformers=True
        )
        # Check with max=11 and min=9 and strict_n_conformers should be OK
        tkw._check_n_conformers(
            mol, "nocharge", min_confs=9, max_confs=11, strict_n_conformers=True
        )
        # Check with min=9 and max=9 and strict_n_conformers should raise an error
        with pytest.raises(
            IncorrectNumConformersError,
            match="has 10 conformers, but charge method 'nocharge' expects exactly 9",
        ):
            tkw._check_n_conformers(
                mol, "nocharge", min_confs=9, max_confs=9, strict_n_conformers=True
            )
        # Check with min=1 and max=9 and strict_n_conformers should raise an error
        with pytest.raises(
            IncorrectNumConformersError,
            match="has 10 conformers, but charge method 'nocharge' expects between 1 and 9",
        ):
            tkw._check_n_conformers(
                mol, "nocharge", min_confs=1, max_confs=9, strict_n_conformers=True
            )
        # Check with min=11 and max=12 and strict_n_conformers should raise an error
        with pytest.raises(
            IncorrectNumConformersError,
            match="has 10 conformers, but charge method 'nocharge' expects between 11 and 12",
        ):
            tkw._check_n_conformers(
                mol, "nocharge", min_confs=11, max_confs=12, strict_n_conformers=True
            )


class TestToolkitRegistry:
    """Test the ToolkitRegistry class"""

    def test_register_empty_toolkit(self):
        """Ensure the default ToolkitRegistry init returns an empty registry"""
        empty_registry = ToolkitRegistry()

        assert empty_registry.registered_toolkits == []
        assert empty_registry.registered_toolkit_versions == {}

    @requires_openeye
    @requires_rdkit
    def test_register_imported_toolkit_wrappers(self):
        """Test that imported toolkits are registered, and in the expected order"""
        # Ensure a specified order is respected
        default_registry = ToolkitRegistry(
            toolkit_precedence=[
                OpenEyeToolkitWrapper,
                RDKitToolkitWrapper,
                AmberToolsToolkitWrapper,
                BuiltInToolkitWrapper,
            ],
            _register_imported_toolkit_wrappers=True,
        )

        assert len(default_registry.registered_toolkits) == 4

        expected_toolkits = [
            OpenEyeToolkitWrapper,
            RDKitToolkitWrapper,
            AmberToolsToolkitWrapper,
            BuiltInToolkitWrapper,
        ]
        for found, expected in zip(
            default_registry.registered_toolkits, expected_toolkits
        ):
            assert isinstance(found, expected)

        # Test forcing a non-default order
        non_default_registry = ToolkitRegistry(
            toolkit_precedence=[BuiltInToolkitWrapper, RDKitToolkitWrapper],
            _register_imported_toolkit_wrappers=True,
        )

        assert len(non_default_registry.registered_toolkits) == 2

        expected_toolkits = [BuiltInToolkitWrapper, RDKitToolkitWrapper]
        for found, expected in zip(
            non_default_registry.registered_toolkits, expected_toolkits
        ):
            assert isinstance(found, expected)

    @requires_rdkit
    def test_add_bad_toolkit(self):
        registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])
        with pytest.raises(InvalidToolkitError):
            registry.add_toolkit("rdkit as a string")

    @requires_rdkit
    @pytest.mark.skipif(
        OpenEyeToolkitWrapper.is_available(),
        reason="Skipping while OpenEye is available",
    )
    def test_register_unavailable_toolkit(self):
        registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])
        with pytest.raises(ToolkitUnavailableException):
            registry.register_toolkit(
                toolkit_wrapper=OpenEyeToolkitWrapper, exception_if_unavailable=True
            )

    @pytest.mark.skipif(
        RDKitToolkitWrapper.is_available(),
        reason="Skipping while The RDKit is available",
    )
    def test_requires_toolkit_exception(self):
        """Test that ToolkitUnavailableException, not LicenseError, is raised
        when RDKitToolkitWrapper is unavailable"""
        registry = ToolkitRegistry()
        with pytest.raises(ToolkitUnavailableException):
            registry.register_toolkit(
                toolkit_wrapper=RDKitToolkitWrapper, exception_if_unavailable=True
            )

    @requires_openeye
    def test_register_openeye(self):
        """Test creation of toolkit registry with OpenEye toolkit"""
        # Test registration of OpenEyeToolkitWrapper
        toolkit_precedence = [OpenEyeToolkitWrapper]
        registry = ToolkitRegistry(
            toolkit_precedence=toolkit_precedence,
        )

        assert set(type(c) for c in registry.registered_toolkits) == set(
            [OpenEyeToolkitWrapper]
        )

        # Test ToolkitRegistry.resolve()
        assert (
            registry.resolve("to_smiles") == registry.registered_toolkits[0].to_smiles
        )

        # Test ToolkitRegistry.call()
        smiles = "[H]C([H])([H])C([H])([H])[H]"
        molecule = registry.call("from_smiles", smiles)
        smiles2 = registry.call("to_smiles", molecule)
        assert smiles == smiles2

    @requires_rdkit
    def test_register_rdkit(self):
        """Test creation of toolkit registry with RDKit toolkit"""
        # Test registration of RDKitToolkitWrapper
        toolkit_precedence = [RDKitToolkitWrapper]
        registry = ToolkitRegistry(
            toolkit_precedence=toolkit_precedence,
        )

        assert set([type(c) for c in registry.registered_toolkits]) == set(
            [RDKitToolkitWrapper]
        )

        # Test ToolkitRegistry.resolve()
        assert (
            registry.resolve("to_smiles") == registry.registered_toolkits[0].to_smiles
        )

        # Test ToolkitRegistry.call()
        smiles = "[H][C]([H])([H])[C]([H])([H])[H]"
        molecule = registry.call("from_smiles", smiles)
        smiles2 = registry.call("to_smiles", molecule)
        assert smiles == smiles2

    @requires_ambertools
    def test_register_ambertools(self):
        """Test creation of toolkit registry with AmberToolsToolkitWrapper"""
        # Test registration of AmberToolsToolkitWrapper
        toolkit_precedence = [AmberToolsToolkitWrapper]
        registry = ToolkitRegistry(
            toolkit_precedence=toolkit_precedence,
        )

        assert set([type(c) for c in registry.registered_toolkits]) == set(
            [AmberToolsToolkitWrapper]
        )

        # Test ToolkitRegistry.resolve()
        registry.resolve("assign_partial_charges")
        assert (
            registry.resolve("assign_partial_charges")
            == registry.registered_toolkits[0].assign_partial_charges
        )

        # Test ToolkitRegistry.call()
        molecule = RDKitToolkitWrapper().from_file(
            file_path=get_data_file_path("molecules/ethanol.sdf"), file_format="SDF"
        )[0]
        registry.call("assign_partial_charges", molecule)
        charges_from_registry = molecule.partial_charges
        AmberToolsToolkitWrapper().assign_partial_charges(molecule)
        charges_from_toolkit = molecule.partial_charges

        assert np.allclose(charges_from_registry, charges_from_toolkit)

    @requires_ambertools
    def test_register_rdkit_and_ambertools(self):
        """Test creation of toolkit registry with RDKitToolkitWrapper and
        AmberToolsToolkitWrapper and test ToolkitRegistry.resolve()"""
        toolkit_precedence = [RDKitToolkitWrapper, AmberToolsToolkitWrapper]
        registry = ToolkitRegistry(
            toolkit_precedence=toolkit_precedence,
        )

        assert set([type(c) for c in registry.registered_toolkits]) == set(
            [RDKitToolkitWrapper, AmberToolsToolkitWrapper]
        )

        # Resolve to a method that is supported by AmberToolsToolkitWrapper
        # but _not_ RDKitToolkitWrapper. Note that this may change as more
        # functionality is added to to toolkit wrappers
        assert (
            registry.resolve("assign_fractional_bond_orders")
            == registry.registered_toolkits[1].assign_fractional_bond_orders
        )
        # Resolve a method supported by both to the highest-priority wrapper
        assert (
            registry.resolve("from_smiles")
            == registry.registered_toolkits[0].from_smiles
        )

        # Test ToolkitRegistry.call() for each toolkit
        smiles = "[H][C]([H])([H])[C]([H])([H])[H]"
        molecule = registry.call("from_smiles", smiles)
        smiles2 = registry.call("to_smiles", molecule)

        # Round-tripping SMILES is not 100% reliable, so just ensure it returned something
        assert isinstance(smiles2, str)

        # This method is available in AmberToolsToolkitWrapper, but not RDKitToolkitWrapper
        registry.call("assign_partial_charges", molecule)

    @requires_ambertools
    def test_deregister_toolkit(self):
        """Test removing an instantiated toolkit from the registry"""
        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )

        assert any(
            [
                isinstance(tk, AmberToolsToolkitWrapper)
                for tk in toolkit_registry._toolkits
            ]
        )
        assert any(
            [isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits]
        )

        toolkit_registry.deregister_toolkit(toolkit_registry._toolkits[-1])
        assert any(
            [
                isinstance(tk, AmberToolsToolkitWrapper)
                for tk in toolkit_registry._toolkits
            ]
        )
        assert not any(
            [isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits]
        )

        toolkit_registry.deregister_toolkit(toolkit_registry._toolkits[-1])
        assert not any(
            [
                isinstance(tk, AmberToolsToolkitWrapper)
                for tk in toolkit_registry._toolkits
            ]
        )
        assert not any(
            [isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits]
        )

    @requires_ambertools
    def test_deregister_toolkit_by_class(self):
        """Test removing a toolkit from the registry by matching class types"""
        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )

        assert any(
            [
                isinstance(tk, AmberToolsToolkitWrapper)
                for tk in toolkit_registry._toolkits
            ]
        )
        assert any(
            [isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits]
        )

        toolkit_registry.deregister_toolkit(RDKitToolkitWrapper)
        assert any(
            [
                isinstance(tk, AmberToolsToolkitWrapper)
                for tk in toolkit_registry._toolkits
            ]
        )
        assert not any(
            [isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits]
        )

        toolkit_registry.deregister_toolkit(AmberToolsToolkitWrapper)
        assert not any(
            [
                isinstance(tk, AmberToolsToolkitWrapper)
                for tk in toolkit_registry._toolkits
            ]
        )
        assert not any(
            [isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits]
        )

    @requires_ambertools
    def test_deregister_toolkit_bad_inputs(self):
        """Test bad inputs to deregister_toolkit"""
        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper]
        )

        with pytest.raises(InvalidToolkitError):
            toolkit_registry.deregister_toolkit("rdkit as a string")

        # Attempt to deregister a toolkit that is not registered
        with pytest.raises(ToolkitUnavailableException):
            toolkit_registry.deregister_toolkit(RDKitToolkitWrapper)

    def deregister_from_global_registry(self):
        # TODO: Update this, or move into a separate TestClass, pending GLOBAL_TOOLKIT_REGISTRY rewor
        # See issue #493
        # Whatever the first tookit it, de-register it and verify it's de-registered

        # Keep a copy of the original registry since this is a "global" variable accessible to other modules
        from copy import deepcopy

        global_registry_copy = deepcopy(GLOBAL_TOOLKIT_REGISTRY)
        first_toolkit = type(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits[0])
        num_toolkits = len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits)

        GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(first_toolkit)

        assert first_toolkit not in [
            type(tk) for tk in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits
        ]
        assert len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits) == num_toolkits - 1

        GLOBAL_TOOLKIT_REGISTRY = deepcopy(global_registry_copy)

    def test_register_builtintoolkit(self):
        """Test creation of toolkit registry with Built-in toolkit"""
        # Test registration of BuiltInToolkitWrapper
        toolkit_precedence = [BuiltInToolkitWrapper]
        registry = ToolkitRegistry(
            toolkit_precedence=toolkit_precedence,
        )
        # registry.register_toolkit(BuiltInToolkitWrapper)
        assert set([type(c) for c in registry.registered_toolkits]) == set(
            [BuiltInToolkitWrapper]
        )

        # Test ToolkitRegistry.resolve()
        assert (
            registry.resolve("assign_partial_charges")
            == registry.registered_toolkits[0].assign_partial_charges
        )

    @requires_rdkit
    @requires_openeye
    def test_toolkit_versions(self):
        """Test behavior of ToolkitRegistry.registered_toolkit_versions"""
        toolkit_precedence = [
            OpenEyeToolkitWrapper,
            RDKitToolkitWrapper,
            AmberToolsToolkitWrapper,
            BuiltInToolkitWrapper,
        ]
        all_toolkits = ToolkitRegistry(toolkit_precedence=toolkit_precedence)
        versions = all_toolkits.registered_toolkit_versions

        import openeye
        import rdkit

        assert versions["OpenEye Toolkit"] == openeye.__version__
        assert versions["The RDKit"] == rdkit.__version__
        assert versions["AmberTools"].startswith(
            "2"
        )  # TODO: Safer way of checking AmberTools version
        assert versions["Built-in Toolkit"] is None

        toolkit_precedence = [
            RDKitToolkitWrapper,
            AmberToolsToolkitWrapper,
            BuiltInToolkitWrapper,
        ]
        no_openeye = ToolkitRegistry(toolkit_precedence=toolkit_precedence)

        assert "OpenEye Toolkit" not in no_openeye.registered_toolkit_versions.keys()

    @requires_ambertools
    def test_call_raise_first_error(self):
        """Test to ensure proper behavior of raise_first_error kwarg to ToolkitRegistry.call"""
        toolkit_precedence = [
            BuiltInToolkitWrapper,
            RDKitToolkitWrapper,
            AmberToolsToolkitWrapper,
        ]
        registry = ToolkitRegistry(
            toolkit_precedence=toolkit_precedence,
        )
        mol = registry.call("from_smiles", "C")
        # Specify that the ToolkitRegistry should raise the first ChargeMethodUnavailableError it encounters
        with pytest.raises(
            ChargeMethodUnavailableError,
            match='"notarealchargemethod"" is not supported by the Built-in toolkit.',
        ):
            registry.call(
                "assign_partial_charges",
                molecule=mol,
                partial_charge_method="NotARealChargeMethod",
                raise_exception_types=[ChargeMethodUnavailableError],
            )
        # Specify that the ToolkitRegistry should collect all the errors it encounters and
        # ensure it raises a single ValueError when no ToolkitWrappers succeed
        with pytest.raises(
            ValueError,
            match="partial_charge_method 'notarealchargemethod' is not available from AmberToolsToolkitWrapper",
        ):
            registry.call(
                "assign_partial_charges",
                molecule=mol,
                partial_charge_method="NotARealChargeMethod",
                raise_exception_types=[],
            )


@requires_openeye
def test_license_check(monkeypatch):
    def MockIsLicensed():
        return False

    from openeye import oeiupac

    assert oeiupac.OEIUPACIsLicensed()

    from openff.toolkit.utils.toolkits import OpenEyeToolkitWrapper

    assert OpenEyeToolkitWrapper.is_available()

    # Mock OEIUPACIsLicensed to return False ...
    monkeypatch.setattr(oeiupac, "OEIUPACIsLicensed", MockIsLicensed)

    # ... ensure that the oeiupac module reflects this
    assert not oeiupac.OEIUPACIsLicensed()

    # ... and ensure that the toolkit wrapper is **still** available
    assert OpenEyeToolkitWrapper()._check_licenses()
    assert OpenEyeToolkitWrapper().is_available()

    from openff.toolkit.utils.toolkits import requires_openeye_module

    @requires_openeye_module("oeszybki")
    def func_using_extraneous_openeye_module():
        pass

    with pytest.raises(Exception, match="currently use oeszybki"):
        func_using_extraneous_openeye_module()

    @requires_openeye_module("oeiupac")
    def func_using_unlicsensed_openeye_module():
        pass

    with pytest.raises(AssertionError):
        func_using_unlicsensed_openeye_module()
