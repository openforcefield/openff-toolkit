"""
Tests for cheminformatics toolkit wrappers

"""

import logging
import os
import pathlib
from tempfile import NamedTemporaryFile
from typing import Dict

import numpy as np
import pytest
from numpy.testing import assert_almost_equal
from openff.units import unit

from openff.toolkit.tests.create_molecules import (
    create_acetaldehyde,
    create_acetate,
    create_cyclic_n3h3,
    create_cyclohexane,
    create_ethanol,
    create_reversed_ethanol,
    topology_with_metadata,
)
from openff.toolkit.tests.utils import (
    requires_ambertools,
    requires_openeye,
    requires_rdkit,
)
from openff.toolkit.topology.molecule import Atom, Molecule
from openff.toolkit.utils import get_data_file_path
from openff.toolkit.utils.exceptions import (
    ChargeMethodUnavailableError,
    ConformerGenerationError,
    IncorrectNumConformersError,
    IncorrectNumConformersWarning,
    InvalidIUPACNameError,
    InvalidToolkitError,
    NotAttachedToMoleculeError,
    RadicalsNotSupportedError,
    ToolkitUnavailableException,
    UndefinedStereochemistryError,
)
from openff.toolkit.utils.toolkits import (
    GLOBAL_TOOLKIT_REGISTRY,
    AmberToolsToolkitWrapper,
    BuiltInToolkitWrapper,
    GAFFAtomTypeWarning,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
    ToolkitRegistry,
    ToolkitWrapper,
)


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
        "cis": unit.Quantity(
            np.array(
                [
                    [-0.95927322, -0.91789997, 0.36333418],
                    [-0.34727824, 0.12828046, 0.22784603],
                    [0.82766682, 0.26871252, -0.42284882],
                    [-0.67153811, 1.10376000, 0.61921501],
                    [1.15035689, -0.58282924, -0.78766006],
                ]
            ),
            unit.angstrom,
        ),
        "trans": unit.Quantity(
            np.array(
                [
                    [-0.95927322, -0.91789997, 0.36333418],
                    [-0.34727824, 0.12828046, 0.22784603],
                    [0.82766682, 0.26871252, -0.42284882],
                    [-0.67153811, 1.10376000, 0.61921501],
                    [1.14532626, 1.19679034, -0.41266876],
                ]
            ),
            unit.angstrom,
        ),
    }


@requires_openeye
class TestOpenEyeToolkitWrapper:
    """Test the OpenEyeToolkitWrapper"""

    # TODO: Make separate smiles_add_H and smiles_explicit_H tests

    def test_smiles(self):
        """Test OpenEyeToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()

        # This differs from RDKit's SMILES due to different canonicalization schemes

        smiles = "[H]C([H])([H])C([H])([H])[H]"
        molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)
        # When creating an OFFMol from SMILES, partial charges should be initialized to None
        assert molecule.partial_charges is None
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2

    def test_smiles_missing_stereochemistry(self):
        """Test OpenEyeToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()

        unspec_chiral_smiles = r"C\C(F)=C(/F)CC(C)(Cl)Br"
        spec_chiral_smiles = r"C\C(F)=C(/F)C[C@@](C)(Cl)Br"
        unspec_db_smiles = r"CC(F)=C(F)C[C@@](C)(Cl)Br"
        spec_db_smiles = r"C\C(F)=C(/F)C[C@@](C)(Cl)Br"

        for title, smiles, raises_exception in [
            ("unspec_chiral_smiles", unspec_chiral_smiles, True),
            ("spec_chiral_smiles", spec_chiral_smiles, False),
            ("unspec_db_smiles", unspec_db_smiles, True),
            ("spec_db_smiles", spec_db_smiles, False),
        ]:
            if raises_exception:
                with pytest.raises(UndefinedStereochemistryError):
                    Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)
                Molecule.from_smiles(
                    smiles,
                    toolkit_registry=toolkit_wrapper,
                    allow_undefined_stereo=True,
                )
            else:
                Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)

    def test_openeye_from_smiles_radical(self):
        """Test that parsing an SMILES with a radical raises RadicalsNotSupportedError."""
        with pytest.raises(RadicalsNotSupportedError):
            OpenEyeToolkitWrapper().from_smiles("[CH3]")

    def test_openeye_from_smiles_transition_metal_radical(self):
        """Test that parsing an SMILES with a transition metal radical works."""
        OpenEyeToolkitWrapper().from_smiles("[Zn+2]")

    # TODO: test_smiles_round_trip

    def test_smiles_add_H(self):
        """Test OpenEyeToolkitWrapper for adding explicit hydrogens"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        # This differs from RDKit's SMILES due to different canonicalization schemes
        input_smiles = "CC"
        expected_output_smiles = "[H]C([H])([H])C([H])([H])[H]"
        molecule = Molecule.from_smiles(input_smiles, toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert expected_output_smiles == smiles2

    def test_smiles_charged(self):
        """Test OpenEyeToolkitWrapper functions for reading/writing charged SMILES"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        # This differs from RDKit's expected output due to different canonicalization schemes
        smiles = "[H]C([H])([H])[N+]([H])([H])[H]"
        molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2

    def test_to_from_openeye_core_props_filled(self):
        """Test OpenEyeToolkitWrapper to_openeye() and from_openeye()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()

        # Replacing with a simple molecule with stereochemistry
        input_smiles = r"C\C(F)=C(/F)C[C@@](C)(Cl)Br"
        expected_output_smiles = (
            r"[H]C([H])([H])/C(=C(/C([H])([H])[C@@](C([H])([H])[H])(Cl)Br)\F)/F"
        )
        molecule = Molecule.from_smiles(input_smiles, toolkit_registry=toolkit_wrapper)
        assert (
            molecule.to_smiles(toolkit_registry=toolkit_wrapper)
            == expected_output_smiles
        )

        # Populate core molecule property fields
        molecule.name = "Alice"
        partial_charges = unit.Quantity(
            np.array(
                [
                    -0.9,
                    -0.8,
                    -0.7,
                    -0.6,
                    -0.5,
                    -0.4,
                    -0.3,
                    -0.2,
                    -0.1,
                    0.0,
                    0.1,
                    0.2,
                    0.3,
                    0.4,
                    0.5,
                    0.6,
                    0.7,
                    0.8,
                ]
            ),
            unit.elementary_charge,
        )
        molecule.partial_charges = partial_charges
        coords = unit.Quantity(
            np.array(
                [
                    ["0.0", "1.0", "2.0"],
                    ["3.0", "4.0", "5.0"],
                    ["6.0", "7.0", "8.0"],
                    ["9.0", "10.0", "11.0"],
                    ["12.0", "13.0", "14.0"],
                    ["15.0", "16.0", "17.0"],
                    ["18.0", "19.0", "20.0"],
                    ["21.0", "22.0", "23.0"],
                    ["24.0", "25.0", "26.0"],
                    ["27.0", "28.0", "29.0"],
                    ["30.0", "31.0", "32.0"],
                    ["33.0", "34.0", "35.0"],
                    ["36.0", "37.0", "38.0"],
                    ["39.0", "40.0", "41.0"],
                    ["42.0", "43.0", "44.0"],
                    ["45.0", "46.0", "47.0"],
                    ["48.0", "49.0", "50.0"],
                    ["51.0", "52.0", "53.0"],
                ]
            ),
            unit.angstrom,
        )
        molecule.add_conformer(coords)
        # Populate core atom property fields
        molecule.atoms[2].name = "Bob"
        # Ensure one atom has its stereochemistry specified
        central_carbon_stereo_specified = False
        for atom in molecule.atoms:
            if (atom.atomic_number == 6) and atom.stereochemistry == "S":
                central_carbon_stereo_specified = True
        assert central_carbon_stereo_specified

        # Populate bond core property fields
        fractional_bond_orders = [float(val) for val in range(1, 19)]
        for fbo, bond in zip(fractional_bond_orders, molecule.bonds):
            bond.fractional_bond_order = fbo

        # Do a first conversion to/from oemol
        oemol = molecule.to_openeye()
        molecule2 = Molecule.from_openeye(oemol)

        # Test that properties survived first conversion
        # assert molecule.to_dict() == molecule2.to_dict()
        assert molecule.name == molecule2.name
        # NOTE: This expects the same indexing scheme in the original and new molecule

        central_carbon_stereo_specified = False
        for atom in molecule2.atoms:
            if (atom.atomic_number == 6) and atom.stereochemistry == "S":
                central_carbon_stereo_specified = True
        assert central_carbon_stereo_specified
        for atom1, atom2 in zip(molecule.atoms, molecule2.atoms):
            # OpenEye wrapper always adds hierarchy metadata (residue name + num) info, so account for that
            atom1_dict = atom1.to_dict()
            atom1_dict["metadata"].update(
                {
                    "residue_name": "UNL",
                    "residue_number": 1,
                    "insertion_code": " ",
                    "chain_id": " ",
                }
            )
            assert atom1_dict == atom2.to_dict()
        for bond1, bond2 in zip(molecule.bonds, molecule2.bonds):
            assert bond1.to_dict() == bond2.to_dict()
        assert (molecule.conformers[0] == molecule2.conformers[0]).all()
        for pc1, pc2 in zip(molecule._partial_charges, molecule2._partial_charges):
            pc1_ul = pc1.m_as(unit.elementary_charge)
            pc2_ul = pc2.m_as(unit.elementary_charge)
            assert_almost_equal(pc1_ul, pc2_ul, decimal=6)
        assert (
            molecule2.to_smiles(toolkit_registry=toolkit_wrapper)
            == expected_output_smiles
        )

    def test_to_from_openeye_core_props_unset(self):
        """Test OpenEyeToolkitWrapper to_openeye() and from_openeye() when given empty core property fields"""
        toolkit_wrapper = OpenEyeToolkitWrapper()

        # Using a simple molecule with tetrahedral and bond stereochemistry
        input_smiles = r"C\C(F)=C(/F)C[C@](C)(Cl)Br"

        expected_output_smiles = (
            r"[H]C([H])([H])/C(=C(/C([H])([H])[C@](C([H])([H])[H])(Cl)Br)\F)/F"
        )
        molecule = Molecule.from_smiles(input_smiles, toolkit_registry=toolkit_wrapper)
        assert (
            molecule.to_smiles(toolkit_registry=toolkit_wrapper)
            == expected_output_smiles
        )

        # Ensure one atom has its stereochemistry specified
        central_carbon_stereo_specified = False
        for atom in molecule.atoms:
            if (atom.atomic_number == 6) and atom.stereochemistry == "R":
                central_carbon_stereo_specified = True
        assert central_carbon_stereo_specified

        # Do a first conversion to/from oemol
        oemol = molecule.to_openeye()
        molecule2 = Molecule.from_openeye(oemol)

        # Test that properties survived first conversion
        assert molecule.name == molecule2.name
        # NOTE: This expects the same indexing scheme in the original and new molecule

        central_carbon_stereo_specified = False
        for atom in molecule2.atoms:
            if (atom.atomic_number == 6) and atom.stereochemistry == "R":
                central_carbon_stereo_specified = True
        assert central_carbon_stereo_specified
        for atom1, atom2 in zip(molecule.atoms, molecule2.atoms):
            # OpenEye wrapper always adds hierarchy metadata (residue name + num) info, so account for that
            atom1_dict = atom1.to_dict()
            atom1_dict["metadata"].update(
                {
                    "residue_name": "UNL",
                    "residue_number": 1,
                    "insertion_code": " ",
                    "chain_id": " ",
                }
            )
            assert atom1_dict == atom2.to_dict()
        for bond1, bond2 in zip(molecule.bonds, molecule2.bonds):
            assert bond1.to_dict() == bond2.to_dict()
        # The molecule was initialized from SMILES, so mol.conformers arrays should be None for both
        assert molecule.conformers is None
        assert molecule2.conformers is None
        # The molecule was initialized from SMILES, so mol.partial_charges arrays should be None for both
        assert molecule.partial_charges is None
        assert molecule2.partial_charges is None

        assert (
            molecule2.to_smiles(toolkit_registry=toolkit_wrapper)
            == expected_output_smiles
        )

    def test_to_from_openeye_none_partial_charges(self):
        """Test to ensure that to_openeye and from_openeye correctly handle None partial charges"""
        import math

        # Create ethanol, which has partial charges defined with float values
        ethanol = create_ethanol()
        assert ethanol.partial_charges is not None
        # Convert to OEMol, which should populate the partial charges on
        # the OEAtoms with the same partial charges
        oemol = ethanol.to_openeye()
        for oeatom in oemol.GetAtoms():
            assert not math.isnan(oeatom.GetPartialCharge())
        # Change the first OEAtom's partial charge to nan, and ensure that it comes
        # back to OFFMol with only the first atom as nan
        for oeatom in oemol.GetAtoms():
            oeatom.SetPartialCharge(float("nan"))
            break
        eth_from_oe = Molecule.from_openeye(oemol)
        assert math.isnan(eth_from_oe.partial_charges[0].m_as(unit.elementary_charge))
        for pc in eth_from_oe.partial_charges[1:]:
            assert not math.isnan(pc.m_as(unit.elementary_charge))
        # Then, set all the OEMol's partial charges to nan, and ensure that
        # from_openeye produces an OFFMol with partial_charges = None
        for oeatom in oemol.GetAtoms():
            oeatom.SetPartialCharge(float("nan"))
        eth_from_oe = Molecule.from_openeye(oemol)
        assert eth_from_oe.partial_charges is None

        # Send the OFFMol with partial_charges = None back to OEMol, and
        # ensure that all its charges are nan
        oemol2 = eth_from_oe.to_openeye()
        for oeatom in oemol2.GetAtoms():
            assert math.isnan(oeatom.GetPartialCharge())

    def test_to_from_openeye_hierarchy_metadata(self):
        """
        Test roundtripping to/from ``OpenEyeToolkitWrapper`` for molecules with PDB hierarchy metadata
        """
        from openeye import oechem

        for molecule in topology_with_metadata().molecules:
            oemol = molecule.to_openeye()
            roundtrip_mol = Molecule.from_openeye(oemol)

            # Check OEMol
            for orig_atom, oe_atom in zip(molecule.atoms, oemol.GetAtoms()):
                if "residue_name" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["residue_name"]
                        == oechem.OEAtomGetResidue(oe_atom).GetName()
                    )

                if "residue_number" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["residue_number"]
                        == oechem.OEAtomGetResidue(oe_atom).GetResidueNumber()
                    )

                if "insertion_code" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["insertion_code"]
                        == oechem.OEAtomGetResidue(oe_atom).GetInsertCode()
                    )

                if "chain_id" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["chain_id"]
                        == oechem.OEAtomGetResidue(oe_atom).GetChainID()
                    )

            # Check roundtripped OFFMol
            for orig_atom, roundtrip_atom in zip(molecule.atoms, roundtrip_mol.atoms):
                # If ANY atom in the OEMol has any hierarchy metadata set, then the ENTIRE molecule is considered
                # to have metadata. Anything that wasn't set defaults to ("UNK", 1, " ").
                if oechem.OEHasResidues(oemol):
                    if "residue_name" in orig_atom.metadata:
                        assert (
                            orig_atom.metadata["residue_name"]
                            == roundtrip_atom.metadata["residue_name"]
                        )

                    if "residue_number" in orig_atom.metadata:
                        assert (
                            orig_atom.metadata["residue_number"]
                            == roundtrip_atom.metadata["residue_number"]
                        )

                    if "insertion_code" in orig_atom.metadata:
                        assert (
                            orig_atom.metadata["insertion_code"]
                            == roundtrip_atom.metadata["insertion_code"]
                        )

                    if "chain_id" in orig_atom.metadata:
                        assert (
                            orig_atom.metadata["chain_id"]
                            == roundtrip_atom.metadata["chain_id"]
                        )

                else:
                    assert "residue_name" not in roundtrip_atom.metadata
                    assert "residue_number" not in roundtrip_atom.metadata
                    assert "insertion_code" not in roundtrip_atom.metadata
                    assert "chain_id" not in roundtrip_atom.metadata

    def test_from_openeye_mutable_input(self):
        """
        Test ``OpenEyeToolkitWrapper.from_openeye`` does not mutate the input molecule.
        """
        from openeye import oechem

        oe_molecule = oechem.OEMol()
        oechem.OESmilesToMol(oe_molecule, "C")

        assert oechem.OEHasImplicitHydrogens(oe_molecule)
        Molecule.from_openeye(oe_molecule)
        assert oechem.OEHasImplicitHydrogens(oe_molecule)

    def test_from_openeye_radical(self):
        """Test that parsing an oemol with a radical raises RadicalsNotSupportedError."""
        from openeye import oechem

        smiles = "[H][C]([H])[H]"
        oemol = oechem.OEGraphMol()
        oechem.OESmilesToMol(oemol, smiles)

        with pytest.raises(RadicalsNotSupportedError):
            OpenEyeToolkitWrapper().from_openeye(oemol)

    def test_from_openeye_transition_metal_radical(self):
        """Test that parsing an oemol with a transition metal radical works."""
        from openeye import oechem

        smiles = "[Zn+2]"
        oemol = oechem.OEGraphMol()
        oechem.OESmilesToMol(oemol, smiles)

        OpenEyeToolkitWrapper().from_openeye(oemol)

    def test_from_openeye_implicit_hydrogen(self):
        """
        Test OpenEyeToolkitWrapper for loading a molecule with implicit
        hydrogens (correct behavior is to add them explicitly)
        """
        from openeye import oechem

        smiles_impl = "C#C"
        oemol_impl = oechem.OEMol()
        oechem.OESmilesToMol(oemol_impl, smiles_impl)
        molecule_from_impl = Molecule.from_openeye(oemol_impl)

        assert molecule_from_impl.n_atoms == 4

        smiles_expl = "HC#CH"
        oemol_expl = oechem.OEMol()
        oechem.OESmilesToMol(oemol_expl, smiles_expl)
        molecule_from_expl = Molecule.from_openeye(oemol_expl)
        assert molecule_from_expl.to_smiles() == molecule_from_impl.to_smiles()

    def test_openeye_from_smiles_hydrogens_are_explicit(self):
        """
        Test to ensure that OpenEyeToolkitWrapper.from_smiles has the proper behavior with
        respect to its hydrogens_are_explicit kwarg
        """
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles_impl = "C#C"
        with pytest.raises(
            ValueError,
            match="but OpenEye Toolkit interpreted SMILES 'C#C' as having implicit hydrogen",
        ):
            offmol = Molecule.from_smiles(
                smiles_impl,
                toolkit_registry=toolkit_wrapper,
                hydrogens_are_explicit=True,
            )
        offmol = Molecule.from_smiles(
            smiles_impl, toolkit_registry=toolkit_wrapper, hydrogens_are_explicit=False
        )
        assert offmol.n_atoms == 4

        smiles_expl = "HC#CH"
        offmol = Molecule.from_smiles(
            smiles_expl, toolkit_registry=toolkit_wrapper, hydrogens_are_explicit=True
        )
        assert offmol.n_atoms == 4
        # It's debatable whether this next function should pass. Strictly speaking, the hydrogens in this SMILES
        # _are_ explicit, so allowing "hydrogens_are_explicit=False" through here is allowing a contradiction.
        # We might rethink the name of this kwarg.

        offmol = Molecule.from_smiles(
            smiles_expl, toolkit_registry=toolkit_wrapper, hydrogens_are_explicit=False
        )
        assert offmol.n_atoms == 4

    @pytest.mark.parametrize(
        "smiles, expected_map", [("[Cl:1][H]", {0: 1}), ("[Cl:1][H:2]", {0: 1, 1: 2})]
    )
    def test_from_openeye_atom_map(self, smiles, expected_map):
        """
        Test OpenEyeToolkitWrapper for loading a molecule with implicit
        hydrogens (correct behavior is to add them explicitly)
        """
        from openeye import oechem

        oemol = oechem.OEMol()
        oechem.OESmilesToMol(oemol, smiles)

        off_molecule = Molecule.from_openeye(oemol)
        assert off_molecule.properties["atom_map"] == expected_map

    @pytest.mark.parametrize("molecule", get_mini_drug_bank(OpenEyeToolkitWrapper))
    def test_to_inchi(self, molecule):
        """Test, but do not validate, conversion to standard and non-standard InChI"""

        toolkit = OpenEyeToolkitWrapper()
        molecule.to_inchi(toolkit_registry=toolkit)
        molecule.to_inchi(True, toolkit_registry=toolkit)

    @pytest.mark.parametrize("molecule", get_mini_drug_bank(OpenEyeToolkitWrapper))
    def test_to_inchikey(self, molecule):
        """Test, but do not validate, the conversion to standard and non-standard InChIKey"""

        toolkit = OpenEyeToolkitWrapper()
        molecule.to_inchikey(toolkit_registry=toolkit)
        molecule.to_inchikey(True, toolkit_registry=toolkit)

    def test_from_bad_inchi(self):
        """Test building a molecule from a bad InChI string"""

        toolkit = OpenEyeToolkitWrapper()
        inchi = "InChI=1S/ksbfksfksfksbfks"
        with pytest.raises(RuntimeError):
            Molecule.from_inchi(inchi, toolkit_registry=toolkit)

    @pytest.mark.parametrize(
        "molecule",
        get_mini_drug_bank(
            OpenEyeToolkitWrapper,
            xfail_mols={
                "DrugBank_3046": "Molecule is corrupted and is interpreted as a radical",
                "DrugBank_3655": "Molecule is corrupted and is interpreted as a radical",
                "DrugBank_1594": "Molecule is corrupted and is interpreted as a radical",
                "DrugBank_4346": "Molecule is corrupted and is interpreted as a radical",
                "DrugBank_6947": "Molecule is corrupted and is interpreted as a radical",
            },
        ),
    )
    def test_non_standard_inchi_round_trip(self, molecule):
        """Test if a molecule can survive an InChi round trip test in some cases the standard InChI
        will not enough to ensure information is preserved so we test the non-standard inchi here.
        """

        from openff.toolkit.utils.toolkits import UndefinedStereochemistryError

        toolkit = OpenEyeToolkitWrapper()
        inchi = molecule.to_inchi(fixed_hydrogens=True, toolkit_registry=toolkit)
        # make a copy of the molecule from the inchi string
        if molecule.name in openeye_inchi_stereochemistry_lost:
            # some molecules lose sterorchemsitry so they are skipped
            # if we fail here the molecule may of been fixed
            with pytest.raises(UndefinedStereochemistryError):
                mol2 = molecule.from_inchi(inchi, toolkit_registry=toolkit)

        else:
            mol2 = molecule.from_inchi(inchi, toolkit_registry=toolkit)
            # compare the full molecule excluding the properties dictionary
            # turn of the bond order matching as this could move in the aromatic rings
            if molecule.name in openeye_inchi_isomorphic_fails:
                # Some molecules graphs change during the round trip testing
                # we test quite strict isomorphism here
                with pytest.raises(AssertionError):
                    assert molecule.is_isomorphic_with(
                        mol2, bond_order_matching=False, toolkit_registry=toolkit
                    )
            else:
                assert molecule.is_isomorphic_with(
                    mol2, bond_order_matching=False, toolkit_registry=toolkit
                )

    @pytest.mark.parametrize(
        "molecule",
        get_mini_drug_bank(
            OpenEyeToolkitWrapper,
            xfail_mols={
                "DrugBank_2397": 'OpenEye cannot generate a correct IUPAC name and raises a "Warning: Incorrect name:" or simply return "BLAH".',
                "DrugBank_2543": 'OpenEye cannot generate a correct IUPAC name and raises a "Warning: Incorrect name:" or simply return "BLAH".',
                "DrugBank_2642": 'OpenEye cannot generate a correct IUPAC name and raises a "Warning: Incorrect name:" or simply return "BLAH".',
                "DrugBank_1212": "the roundtrip generates molecules with very different IUPAC/SMILES!",
                "DrugBank_2210": "the roundtrip generates molecules with very different IUPAC/SMILES!",
                "DrugBank_4584": "the roundtrip generates molecules with very different IUPAC/SMILES!",
                "DrugBank_390": 'raises warning "Unable to make OFFMol from OEMol: OEMol has unspecified stereochemistry."',
                "DrugBank_810": 'raises warning "Unable to make OFFMol from OEMol: OEMol has unspecified stereochemistry."',
                "DrugBank_4316": 'raises warning "Unable to make OFFMol from OEMol: OEMol has unspecified stereochemistry."',
                "DrugBank_7124": 'raises warning "Unable to make OFFMol from OEMol: OEMol has unspecified stereochemistry."',
                "DrugBank_3739": 'raises warning "Failed to parse name:"',
                "DrugBank_4346": 'raises warning "Failed to parse name:"',
                "DrugBank_5415": 'raises warning "Failed to parse name:"',
                "DrugBank_1661": "fails roundtrip test",
                "DrugBank_6353": "fails roundtrip test",
                "DrugBank_2799": "from_iupac fails to read what to_iupac returns",
                "DrugBank_4865": "from_iupac fails to read what to_iupac returns",
                "DrugBank_2465": "from_iupac fails to read what to_iupac returns",
            },
        ),
    )
    def test_iupac_round_trip(self, molecule):
        """Test round-trips with IUPAC names"""
        undefined_stereo = molecule.name in openeye_iupac_bad_stereo

        iupac = molecule.to_iupac()

        if undefined_stereo:
            with pytest.raises(UndefinedStereochemistryError):
                Molecule.from_iupac(iupac)

        molecule_copy = Molecule.from_iupac(
            iupac, allow_undefined_stereo=undefined_stereo
        )
        if not undefined_stereo:
            assert molecule.is_isomorphic_with(
                molecule_copy, atom_stereochemistry_matching=not undefined_stereo
            )

    def test_from_iupac_failure(self):
        """Test that invalid IUPAC names are handled properly"""
        toolkit = OpenEyeToolkitWrapper()

        with pytest.raises(InvalidIUPACNameError):
            toolkit.from_iupac(".BETA.-PINENE")

    def test_atom_names_round_trip(self):
        """Test that atom names are correctly preserved in a round trip"""
        # Create a molecule with unique atom names
        molecule = Molecule.from_smiles("c1ccccc1")
        molecule.generate_unique_atom_names()
        # Convert it to an OEMol
        oemol = molecule.to_openeye()
        # Compare atom names
        for oeatom, atom in zip(oemol.GetAtoms(), molecule.atoms):
            assert oeatom.GetName() == atom.name
        # Change the OEMol atom names
        for index, oeatom in enumerate(oemol.GetAtoms()):
            oeatom.SetName(f"X{index}")
        # Round trip back to molecule
        molecule2 = molecule.from_openeye(oemol)
        # Compare atom names
        for oeatom, atom2 in zip(oemol.GetAtoms(), molecule2.atoms):
            assert oeatom.GetName() == atom2.name

    def test_from_pathlib_path(self):
        ethanol = create_ethanol()
        with NamedTemporaryFile(suffix=".sdf") as outfile:
            filename = str(outfile.name)
            ethanol.to_file(filename, file_format="sdf")

            toolkit = OpenEyeToolkitWrapper()
            toolkit.from_file(pathlib.Path(filename), file_format="sdf")

    def test_write_multiconformer_pdb(self):
        """
        Make sure OpenEye can write multi conformer PDB files.
        """
        from io import StringIO

        toolkit = OpenEyeToolkitWrapper()
        # load up a multiconformer sdf file and condense down the conformers
        molecules = Molecule.from_file(
            get_data_file_path("molecules/butane_multi.sdf"), toolkit_registry=toolkit
        )
        butane = molecules.pop(0)
        for mol in molecules:
            butane.add_conformer(mol.conformers[0])
        assert butane.n_conformers == 7
        sio = StringIO()
        butane.to_file(sio, "pdb", toolkit_registry=toolkit)
        # we need to make sure each conformer is wrote to the file
        pdb = sio.getvalue()
        assert pdb.count("END") == 7

    def test_write_pdb_preserving_atom_order(self):
        """
        Make sure OpenEye does not rearrange hydrogens when writing PDBs
        (reference: https://github.com/openforcefield/openff-toolkit/issues/475).
        """
        from io import StringIO

        toolkit = OpenEyeToolkitWrapper()
        water = Molecule()
        water.add_atom(1, 0, False)
        water.add_atom(8, 0, False)
        water.add_atom(1, 0, False)
        water.add_bond(0, 1, 1, False)
        water.add_bond(1, 2, 1, False)
        water.add_conformer(
            unit.Quantity(
                np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
                unit.angstrom,
            )
        )
        sio = StringIO()
        water.to_file(sio, "pdb", toolkit_registry=toolkit)
        water_from_pdb = sio.getvalue()
        water_from_pdb_split = water_from_pdb.split("\n")
        # Check serial number
        assert water_from_pdb_split[0].split()[1].rstrip() == "1"
        assert water_from_pdb_split[1].split()[1].rstrip() == "2"
        assert water_from_pdb_split[2].split()[1].rstrip() == "3"
        # Check atom name
        assert water_from_pdb_split[0].split()[2].rstrip() == "H"
        assert water_from_pdb_split[1].split()[2].rstrip() == "O"
        assert water_from_pdb_split[2].split()[2].rstrip() == "H"

    def test_get_sdf_coordinates(self):
        """Test OpenEyeToolkitWrapper for importing a single set of coordinates from a sdf file"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path("molecules/toluene.sdf")
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule.conformers) == 1
        assert molecule.conformers[0].shape == (15, 3)

    def test_load_multiconformer_sdf_as_separate_molecules(self):
        """
        Test OpenEyeToolkitWrapper for reading a "multiconformer" SDF, which the OFF
        Toolkit should treat as separate molecules
        """
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path("molecules/methane_multiconformer.sdf")
        molecules = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecules) == 2
        assert len(molecules[0].conformers) == 1
        assert len(molecules[1].conformers) == 1
        assert molecules[0].conformers[0].shape == (5, 3)

    def test_load_multiconformer_sdf_as_separate_molecules_properties(self):
        """
        Test OpenEyeToolkitWrapper for reading a "multiconformer" SDF, which the OFF
        Toolkit should treat as separate molecules, and it should load their SD properties
        and partial charges separately
        """
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path("molecules/methane_multiconformer_properties.sdf")
        molecules = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecules) == 2
        assert len(molecules[0].conformers) == 1
        assert len(molecules[1].conformers) == 1
        assert molecules[0].conformers[0].shape == (5, 3)
        # The first molecule in the SDF has the following properties and charges:
        assert molecules[0].properties["test_property_key"] == "test_property_value"
        np.testing.assert_allclose(
            molecules[0].partial_charges.m_as(unit.elementary_charge),
            [-0.108680, 0.027170, 0.027170, 0.027170, 0.027170],
        )
        # The second molecule in the SDF has the following properties and charges:
        assert molecules[1].properties["test_property_key"] == "test_property_value2"
        assert (
            molecules[1].properties["another_test_property_key"]
            == "another_test_property_value"
        )
        np.testing.assert_allclose(
            molecules[1].partial_charges.m_as(unit.elementary_charge),
            [0.027170, 0.027170, 0.027170, 0.027170, -0.108680],
        )

    def test_file_extension_case(self):
        """
        Test round-trips of some file extensions when called directly from the toolkit wrappers,
        including lower- and uppercase file extensions. Note that this test does not ensure
        accuracy, it only tests that reading/writing without raising an exception.
        """
        mols_in = OpenEyeToolkitWrapper().from_file(
            file_path=get_data_file_path("molecules/ethanol.sdf"), file_format="sdf"
        )

        assert len(mols_in) > 0

        mols_in = OpenEyeToolkitWrapper().from_file(
            file_path=get_data_file_path("molecules/ethanol.sdf"), file_format="SDF"
        )

        assert len(mols_in) > 0

    def test_write_sdf_charges(self):
        """Test OpenEyeToolkitWrapper for writing partial charges to a sdf file"""
        from io import StringIO

        toolkit_wrapper = OpenEyeToolkitWrapper()
        ethanol = create_ethanol()
        sio = StringIO()
        ethanol.to_file(sio, "SDF", toolkit_registry=toolkit_wrapper)
        sdf_text = sio.getvalue()
        # The output lines of interest here will look like
        # > <atom.dprop.PartialCharge>
        # -0.400000 -0.300000 -0.200000 -0.100000 0.000010 0.100000 0.200000 0.300000 0.400000
        # Parse the SDF text, grabbing the numeric line above
        sdf_split = sdf_text.split("\n")
        charge_line_found = False
        for line in sdf_split:
            if charge_line_found:
                charges = [float(i) for i in line.split()]
                break
            if "> <atom.dprop.PartialCharge>" in line:
                charge_line_found = True

        # Make sure that a charge line was ever found
        assert charge_line_found

        # Make sure that the charges found were correct
        assert_almost_equal(
            charges, [-0.4, -0.3, -0.2, -0.1, 0.00001, 0.1, 0.2, 0.3, 0.4]
        )

    def test_write_sdf_no_charges(self):
        """Test OpenEyeToolkitWrapper for writing an SDF file without charges"""
        from io import StringIO

        toolkit_wrapper = OpenEyeToolkitWrapper()
        ethanol = create_ethanol()
        ethanol.partial_charges = None
        sio = StringIO()
        ethanol.to_file(sio, "SDF", toolkit_registry=toolkit_wrapper)
        sdf_text = sio.getvalue()
        # In our current configuration, if the OFFMol doesn't have partial charges, we DO NOT want a partial charge
        # block to be written. For reference, it's possible to indicate that a partial charge is not known by writing
        # out "n/a" (or another placeholder) in the partial charge block atoms without charges.
        assert "<atom.dprop.PartialCharge>" not in sdf_text

    def test_sdf_properties_roundtrip(self):
        """Test OpenEyeToolkitWrapper for performing a round trip of a molecule with defined partial charges
        and entries in the properties dict to and from a sdf file"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        ethanol = create_ethanol()
        ethanol.properties["test_property"] = "test_value"
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix=".sdf") as iofile:
            ethanol.to_file(
                iofile.name, file_format="SDF", toolkit_registry=toolkit_wrapper
            )
            ethanol2 = Molecule.from_file(
                iofile.name, file_format="SDF", toolkit_registry=toolkit_wrapper
            )
        np.testing.assert_allclose(
            ethanol.partial_charges.m_as(unit.elementary_charge),
            ethanol2.partial_charges.m_as(unit.elementary_charge),
        )
        assert ethanol2.properties["test_property"] == "test_value"

        # Now test with no properties or charges
        ethanol = create_ethanol()
        ethanol.partial_charges = None
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix=".sdf") as iofile:
            ethanol.to_file(
                iofile.name, file_format="SDF", toolkit_registry=toolkit_wrapper
            )
            ethanol2 = Molecule.from_file(
                iofile.name, file_format="SDF", toolkit_registry=toolkit_wrapper
            )
        assert ethanol2.partial_charges is None
        assert ethanol2.properties == {}

    def test_write_multiconformer_mol_as_sdf(self):
        """
        Test OpenEyeToolkitWrapper for writing a multiconformer molecule to SDF. The OFF toolkit should only
        save the first conformer.
        """
        from io import StringIO

        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path("molecules/ethanol.sdf")
        ethanol = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        ethanol.partial_charges = unit.Quantity(
            np.array([-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]),
            unit.elementary_charge,
        )
        ethanol.properties["test_prop"] = "test_value"
        new_conf = ethanol.conformers[0] + (
            np.ones(ethanol.conformers[0].shape) * unit.angstrom
        )
        ethanol.add_conformer(new_conf)
        sio = StringIO()
        ethanol.to_file(sio, "sdf", toolkit_registry=toolkit_wrapper)
        data = sio.getvalue()
        # In SD format, each molecule ends with "$$$$"
        assert data.count("$$$$") == 1
        # A basic SDF for ethanol would be 27 lines, though the properties add three more
        assert len(data.split("\n")) == 30
        assert "test_prop" in data
        assert "<atom.dprop.PartialCharge>" in data
        # Ensure the first conformer's first atom's X coordinate is in the file
        assert str(ethanol.conformers[0][0][0].to(unit.angstrom))[:5] in data
        # Ensure the SECOND conformer's first atom's X coordinate is NOT in the file
        assert str(ethanol.conformers[1][0][0].to(unit.angstrom))[:5] not in data

    def test_get_mol2_coordinates(self):
        """Test OpenEyeToolkitWrapper for importing a single set of molecule coordinates"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path("molecules/toluene.mol2")
        molecule1 = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule1.conformers) == 1
        assert molecule1.conformers[0].shape == (15, 3)
        assert_almost_equal(
            molecule1.conformers[0][5][1].m_as(unit.angstrom), 22.98, decimal=2
        )

        # Test loading from file-like object
        with open(filename, "r") as infile:
            molecule2 = Molecule(
                infile, file_format="MOL2", toolkit_registry=toolkit_wrapper
            )
        assert molecule1.is_isomorphic_with(molecule2)
        assert len(molecule2.conformers) == 1
        assert molecule2.conformers[0].shape == (15, 3)
        assert_almost_equal(
            molecule2.conformers[0][5][1].m_as(unit.angstrom), 22.98, decimal=2
        )

        # Test loading from gzipped mol2
        import gzip

        with gzip.GzipFile(filename + ".gz", "r") as infile:
            molecule3 = Molecule(
                infile, file_format="MOL2", toolkit_registry=toolkit_wrapper
            )
        assert molecule1.is_isomorphic_with(molecule3)
        assert len(molecule3.conformers) == 1
        assert molecule3.conformers[0].shape == (15, 3)
        assert_almost_equal(
            molecule3.conformers[0][5][1].m_as(unit.angstrom), 22.98, decimal=2
        )

    def test_get_mol2_charges(self):
        """Test OpenEyeToolkitWrapper for importing a mol2 file specifying partial charges"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path("molecules/toluene_charged.mol2")
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule.conformers) == 1
        assert molecule.conformers[0].shape == (15, 3)
        target_charges = unit.Quantity(
            np.array(
                [
                    -0.1342,
                    -0.1271,
                    -0.1271,
                    -0.1310,
                    -0.1310,
                    -0.0765,
                    -0.0541,
                    0.1314,
                    0.1286,
                    0.1286,
                    0.1303,
                    0.1303,
                    0.0440,
                    0.0440,
                    0.0440,
                ]
            ),
            unit.elementary_charge,
        )
        for pc1, pc2 in zip(molecule._partial_charges, target_charges):
            pc1_ul = pc1.m_as(unit.elementary_charge)
            pc2_ul = pc2.m_as(unit.elementary_charge)
            assert_almost_equal(pc1_ul, pc2_ul, decimal=4)

    def test_mol2_charges_roundtrip(self):
        """Test OpenEyeToolkitWrapper for performing a round trip of a molecule with partial charge to and from
        a mol2 file"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        ethanol = create_ethanol()
        # we increase the magnitude of the partial charges here, since mol2 is only
        # written to 4 digits of precision, and the default middle charge for our test ethanol is 1e-5
        ethanol.partial_charges *= 100
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix=".mol2") as iofile:
            ethanol.to_file(
                iofile.name, file_format="mol2", toolkit_registry=toolkit_wrapper
            )
            ethanol2 = Molecule.from_file(
                iofile.name, file_format="mol2", toolkit_registry=toolkit_wrapper
            )
        np.testing.assert_allclose(
            ethanol.partial_charges.m_as(unit.elementary_charge),
            ethanol2.partial_charges.m_as(unit.elementary_charge),
        )

        # Now test with no properties or charges
        ethanol = create_ethanol()
        ethanol.partial_charges = None
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix=".mol2") as iofile:
            ethanol.to_file(
                iofile.name, file_format="mol2", toolkit_registry=toolkit_wrapper
            )
            ethanol2 = Molecule.from_file(
                iofile.name, file_format="mol2", toolkit_registry=toolkit_wrapper
            )
        assert ethanol2.partial_charges is None
        assert ethanol2.properties == {}

    def test_get_mol2_gaff_atom_types(self):
        """Test that a warning is raised OpenEyeToolkitWrapper when it detects GAFF atom types in a mol2 file."""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        mol2_file_path = get_data_file_path("molecules/AlkEthOH_test_filt1_ff.mol2")
        with pytest.warns(GAFFAtomTypeWarning, match="SYBYL"):
            Molecule.from_file(mol2_file_path, toolkit_registry=toolkit_wrapper)

    def test_generate_conformers(self):
        """Test OpenEyeToolkitWrapper generate_conformers()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = "[H]C([H])([H])C([H])([H])[H]"
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers()
        assert molecule.n_conformers != 0
        assert not (molecule.conformers[0] == (0.0 * unit.angstrom)).all()

    def test_generate_multiple_conformers(self):
        """Test OpenEyeToolkitWrapper generate_conformers() for generating multiple conformers"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = "CCCCCCCCCN"
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(
            rms_cutoff=1 * unit.angstrom,
            n_conformers=100,
            toolkit_registry=toolkit_wrapper,
        )
        assert molecule.n_conformers > 1
        assert not (molecule.conformers[0] == (0.0 * unit.angstrom)).all()

        # Ensure rms_cutoff kwarg is working
        molecule2 = toolkit_wrapper.from_smiles(smiles)
        molecule2.generate_conformers(
            rms_cutoff=0.1 * unit.angstrom,
            n_conformers=100,
            toolkit_registry=toolkit_wrapper,
        )
        assert molecule2.n_conformers > molecule.n_conformers

        # Ensure n_conformers kwarg is working
        molecule2 = toolkit_wrapper.from_smiles(smiles)
        molecule2.generate_conformers(
            rms_cutoff=0.1 * unit.angstrom,
            n_conformers=10,
            toolkit_registry=toolkit_wrapper,
        )
        assert molecule2.n_conformers == 10

    def test_generate_conformers_failure(self):
        toolkit = OpenEyeToolkitWrapper()

        molecule = Molecule.from_smiles("F[U](F)(F)(F)(F)F")

        with pytest.raises(ConformerGenerationError, match="Omega conf.*fail"):
            toolkit.generate_conformers(molecule, n_conformers=1)

    def test_apply_elf_conformer_selection(self):
        """Test applying the ELF10 method."""

        toolkit = OpenEyeToolkitWrapper()

        molecule = Molecule.from_file(
            get_data_file_path(os.path.join("molecules", "z_3_hydroxy_propenal.sdf")),
            "SDF",
        )

        # Test that the simple case of no conformers does not yield an exception.
        toolkit.apply_elf_conformer_selection(molecule)

        initial_conformers = [
            # Add a conformer with an internal H-bond.
            unit.Quantity(
                np.array(
                    [
                        [0.5477, 0.3297, -0.0621],
                        [-0.1168, -0.7881, 0.2329],
                        [-1.4803, -0.8771, 0.1667],
                        [-0.2158, 1.5206, -0.4772],
                        [-1.4382, 1.5111, -0.5580],
                        [1.6274, 0.3962, -0.0089],
                        [0.3388, -1.7170, 0.5467],
                        [-1.8612, -0.0347, -0.1160],
                        [0.3747, 2.4222, -0.7115],
                    ]
                ),
                unit.angstrom,
            ),
            # Add a conformer without an internal H-bond.
            unit.Quantity(
                np.array(
                    [
                        [0.5477, 0.3297, -0.0621],
                        [-0.1168, -0.7881, 0.2329],
                        [-1.4803, -0.8771, 0.1667],
                        [-0.2158, 1.5206, -0.4772],
                        [0.3353, 2.5772, -0.7614],
                        [1.6274, 0.3962, -0.0089],
                        [0.3388, -1.7170, 0.5467],
                        [-1.7743, -1.7634, 0.4166],
                        [-1.3122, 1.4082, -0.5180],
                    ]
                ),
                unit.angstrom,
            ),
        ]

        molecule._conformers = [*initial_conformers]

        # Apply ELF10
        toolkit.apply_elf_conformer_selection(molecule)
        elf10_conformers = molecule.conformers

        assert len(elf10_conformers) == 1

        assert np.allclose(
            elf10_conformers[0].m_as(unit.angstrom),
            initial_conformers[1].m_as(unit.angstrom),
        )

    def test_assign_partial_charges_am1bcc(self):
        """Test OpenEyeToolkitWrapper assign_partial_charges() with am1bcc"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()
        molecule.assign_partial_charges(
            partial_charge_method="am1bcc", toolkit_registry=toolkit_registry
        )
        charge_sum = np.sum(molecule.partial_charges)
        abs_charge_sum = np.sum(abs(molecule.partial_charges))

        assert abs(charge_sum) < 1e-10 * unit.elementary_charge
        assert abs_charge_sum > 0.25 * unit.elementary_charge

    def test_assign_partial_charges_am1bcc_no_normalization(self):
        """Test OpenEyeToolkitWrapper assign_partial_charges() with am1bcc, with
        normalize_partial_charges=False"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        # Use a cyclic N3H3 molecule, since the threefold symmetry makes it likely to expose rounding
        # errors. (can we find a cyclic molecule with exactly 3 atoms?)
        molecule = create_cyclic_n3h3()
        molecule.assign_partial_charges(
            partial_charge_method="am1bcc",
            toolkit_registry=toolkit_registry,
            normalize_partial_charges=False,
        )
        charge_sum = np.sum(molecule.partial_charges)
        abs_charge_sum = np.sum(abs(molecule.partial_charges))
        # Rounding error should be on the order of 1e-3

        assert 1e-7 > abs(charge_sum.m_as(unit.elementary_charge)) > 1e-8
        assert abs_charge_sum > 0.25 * unit.elementary_charge

    def test_assign_partial_charges_am1bcc_net_charge(self):
        """Test OpenEyeToolkitWrapper assign_partial_charges() on a molecule with a net +1 charge"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_acetate()
        molecule.assign_partial_charges(
            partial_charge_method="am1bcc", toolkit_registry=toolkit_registry
        )
        charge_sum = np.sum(molecule.partial_charges)
        assert 1e-10 > abs(
            (charge_sum - molecule.total_charge).m_as(unit.elementary_charge)
        )

    def test_assign_partial_charges_am1bcc_wrong_n_confs(self):
        """
        Test OpenEyeToolkitWrapper assign_partial_charges() with am1bcc when requesting to use an incorrect number of
        conformers. This test is a bit shorter than that for AmberToolsToolkitWrapper because OETK uses the
        ELF10 multiconformer method of AM1BCC, which doesn't have a maximum number of conformers.
        """

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(
            n_conformers=2,
            rms_cutoff=0.1 * unit.angstrom,
            toolkit_registry=toolkit_registry,
        )

        # Try again, with strict_n_confs as true, but not including use_confs, so the
        # recommended number of confs will be generated
        molecule.assign_partial_charges(
            partial_charge_method="am1bcc",
            toolkit_registry=toolkit_registry,
            strict_n_conformers=True,
        )

    @pytest.mark.parametrize(
        "partial_charge_method", ["am1bcc", "am1elf10", "am1-mulliken", "gasteiger"]
    )
    def test_assign_partial_charges_neutral(self, partial_charge_method):
        """Test OpenEyeToolkitWrapper assign_partial_charges()"""

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
        )
        charge_sum = np.sum(molecule.partial_charges)
        assert 1.0e-10 > abs(charge_sum.m_as(unit.elementary_charge))

        # Atoms 6 and 7 are hydrogens on the central C. If we don't symmetrize charges they'll have slight differences
        assert 1.0e-10 > abs(
            molecule.partial_charges[6] - molecule.partial_charges[7]
        ).m_as(unit.elementary_charge)

    @pytest.mark.parametrize("partial_charge_method", ["am1bcc", "am1-mulliken"])
    def test_assign_partial_charges_conformer_dependence(self, partial_charge_method):
        """Test OpenEyeToolkitWrapper assign_partial_charges()'s use_conformers kwarg
        to ensure charges are really conformer dependent. Skip Gasteiger because it isn't
        conformer dependent."""
        import copy

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=1)
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
            use_conformers=molecule.conformers,
        )
        pcs1 = copy.deepcopy(molecule.partial_charges)
        molecule._conformers[0][0][0] += 0.2 * unit.angstrom
        molecule._conformers[0][1][1] -= 0.2 * unit.angstrom
        molecule._conformers[0][2][1] += 0.2 * unit.angstrom
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
            use_conformers=molecule.conformers,
        )
        for pc1, pc2 in zip(pcs1, molecule.partial_charges):
            assert abs(pc1 - pc2) > 1.0e-5 * unit.elementary_charge

    @pytest.mark.parametrize(
        "partial_charge_method", ["am1bcc", "am1elf10", "am1-mulliken", "gasteiger"]
    )
    def test_assign_partial_charges_net_charge(self, partial_charge_method):
        """
        Test OpenEyeToolkitWrapper assign_partial_charges() on a molecule with net charge.
        """
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_acetate()
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
        )
        charge_sum = np.sum(molecule.partial_charges)
        assert -1.0e-10 < abs(charge_sum.m_as(unit.elementary_charge) + 1.0)

    def test_assign_partial_charges_bad_charge_method(self):
        """Test OpenEyeToolkitWrapper assign_partial_charges() for a nonexistent charge method"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()

        # Molecule.assign_partial_charges calls the ToolkitRegistry with raise_exception_types = [],
        # which means it will only ever return ValueError
        with pytest.raises(
            ValueError, match="is not available from OpenEyeToolkitWrapper"
        ):
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method="NotARealChargeMethod",
            )

        # ToolkitWrappers raise a specific exception class, so we test that here
        with pytest.raises(
            ChargeMethodUnavailableError,
            match="is not available from OpenEyeToolkitWrapper",
        ):
            OETKW = OpenEyeToolkitWrapper()
            OETKW.assign_partial_charges(
                molecule=molecule, partial_charge_method="NotARealChargeMethod"
            )

    @pytest.mark.parametrize(
        "partial_charge_method,expected_n_confs",
        [("am1bcc", 1), ("am1-mulliken", 1), ("gasteiger", 0)],
    )
    def test_assign_partial_charges_wrong_n_confs(
        self, partial_charge_method, expected_n_confs
    ):
        """
        Test OpenEyeToolkitWrapper assign_partial_charges() when requesting to use an incorrect number of
        conformers
        """
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=2, rms_cutoff=0.01 * unit.angstrom)

        # Try passing in the incorrect number of confs, but without specifying strict_n_conformers,
        # which should produce a warning
        with pytest.warns(
            IncorrectNumConformersWarning,
            match=f"has 2 conformers, but charge method '{partial_charge_method}' "
            f"expects exactly {expected_n_confs}.",
        ):
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method=partial_charge_method,
                use_conformers=molecule.conformers,
                strict_n_conformers=False,
            )

        # Try again, with strict_n_confs as true, but not including use_confs, so the
        # recommended number of confs will be generated
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
            strict_n_conformers=True,
        )

        # Test calling the ToolkitWrapper _indirectly_, though a ToolkitRegistry,
        # which should aggregate any exceptions and bundle all of the messages
        # in a failed task together in a single ValueError.
        with pytest.raises(
            ValueError,
            match=f"has 2 conformers, but charge method '{partial_charge_method}' "
            f"expects exactly {expected_n_confs}.",
        ):
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method=partial_charge_method,
                use_conformers=molecule.conformers,
                strict_n_conformers=True,
            )

        # Test calling the ToolkitWrapper _directly_, passing in the incorrect number of
        # confs, and specify strict_n_conformers, which should produce an IncorrectNumConformersError
        with pytest.raises(
            IncorrectNumConformersError,
            match=f"has 2 conformers, but charge method '{partial_charge_method}' "
            f"expects exactly {expected_n_confs}.",
        ):
            OETKW = OpenEyeToolkitWrapper()
            OETKW.assign_partial_charges(
                molecule=molecule,
                partial_charge_method=partial_charge_method,
                use_conformers=molecule.conformers,
                strict_n_conformers=True,
            )

    def test_assign_partial_charges_failure(self):
        """Test OpenEyeToolkitWrapper assign_partial_charges() on a molecule it cannot assign charges to"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = "[Li+1]"
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)

        # For now, I'm just testing AM1-BCC (will test more when the SMIRNOFF spec for other charges is finalized)
        with pytest.raises(Exception) as excinfo:
            molecule.assign_partial_charges(
                partial_charge_method="am1-bcc", toolkit_registry=toolkit_wrapper
            )
            assert "Unable to assign charges" in str(excinfo)
            assert "OE Error: " in str(excinfo)

    def test_assign_partial_charges_trans_cooh_am1bcc(self):
        """Test OpenEyeToolkitWrapper for computing partial charges for problematic molecules, as exemplified by
        Issue 346 (https://github.com/openforcefield/openff-toolkit/issues/346)"""

        lysine = Molecule.from_smiles("C(CC[NH3+])C[C@@H](C(=O)O)N")
        toolkit_wrapper = OpenEyeToolkitWrapper()
        lysine.generate_conformers(toolkit_registry=toolkit_wrapper)
        lysine.assign_partial_charges(
            partial_charge_method="am1bcc", toolkit_registry=toolkit_wrapper
        )

    @pytest.mark.parametrize(
        "bond_order_model",
        ["am1-wiberg", "am1-wiberg-elf10", "pm3-wiberg", "pm3-wiberg-elf10"],
    )
    @pytest.mark.parametrize(
        "smiles",
        [
            "[H]C([H])([H])C([H])([H])[H]",
            "[H]C([H])([H])[N+]([H])([H])[H]",
            r"C\C(F)=C(/F)C[C@@](C)(Cl)Br",
        ],
    )
    def test_assign_fractional_bond_orders(self, bond_order_model, smiles):
        """Test OpenEyeToolkitWrapper assign_fractional_bond_orders()"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.assign_fractional_bond_orders(
            toolkit_registry=toolkit_wrapper, bond_order_model=bond_order_model
        )
        # TODO: Add test for equivalent Wiberg orders for equivalent bonds

        # Sanity check single bonds.
        assert all(
            0.75 < bond.fractional_bond_order < 1.25
            for bond in molecule.bonds
            if bond.bond_order == 1
        )
        # Sanity check double bonds.
        assert all(
            1.75 < bond.fractional_bond_order < 2.25
            for bond in molecule.bonds
            if bond.bond_order == 2
        )

    def test_assign_fractional_bond_orders_multi_conf(
        self, formic_acid_molecule, formic_acid_conformers
    ):
        """Test that the OpenEyeToolkitWrapper assign_fractional_bond_orders()
        function correctly averages over all conformers."""

        toolkit_wrapper = OpenEyeToolkitWrapper()

        # Compute the WBO from a single conformer.
        formic_acid_molecule.assign_fractional_bond_orders(
            toolkit_registry=toolkit_wrapper,
            bond_order_model="am1-wiberg",
            use_conformers=[formic_acid_conformers["cis"]],
        )
        cis_bond_orders = [
            bond.fractional_bond_order for bond in formic_acid_molecule.bonds
        ]
        formic_acid_molecule.assign_fractional_bond_orders(
            toolkit_registry=toolkit_wrapper,
            bond_order_model="am1-wiberg",
            use_conformers=[formic_acid_conformers["trans"]],
        )
        trans_bond_orders = [
            bond.fractional_bond_order for bond in formic_acid_molecule.bonds
        ]

        # Use the method to average the conformers.
        formic_acid_molecule.assign_fractional_bond_orders(
            toolkit_registry=toolkit_wrapper,
            bond_order_model="am1-wiberg",
            use_conformers=[
                formic_acid_conformers["cis"],
                formic_acid_conformers["trans"],
            ],
        )
        avg_bond_orders = [
            bond.fractional_bond_order for bond in formic_acid_molecule.bonds
        ]

        # The average should be distinct from the WBO from either conformer.
        assert not np.allclose(cis_bond_orders, avg_bond_orders)
        assert not np.allclose(trans_bond_orders, avg_bond_orders)

        assert np.allclose(
            np.mean([trans_bond_orders, cis_bond_orders], axis=0), avg_bond_orders
        )

    def test_assign_fractional_bond_orders_conformer_dependence(self):
        """
        Test that OpenEyeToolkitWrapper assign_fractional_bond_orders() provides different results when using
        different conformers
        """
        toolkit_wrapper = OpenEyeToolkitWrapper()
        # Get the WBOs using one conformer
        molecule = create_ethanol()
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        molecule.assign_fractional_bond_orders(
            toolkit_registry=toolkit_wrapper,
            use_conformers=molecule.conformers,
            bond_order_model="am1-wiberg",
        )

        # Do the same again, but change the conformer to yield a different result
        molecule_diff_coords = create_ethanol()
        molecule_diff_coords.generate_conformers(toolkit_registry=toolkit_wrapper)
        molecule_diff_coords._conformers[0][0][0] = (
            molecule_diff_coords._conformers[0][0][0] + 1.0 * unit.angstrom
        )
        molecule_diff_coords._conformers[0][1][0] = (
            molecule_diff_coords._conformers[0][1][0] - 1.0 * unit.angstrom
        )
        molecule_diff_coords._conformers[0][2][0] = (
            molecule_diff_coords._conformers[0][2][0] + 1.0 * unit.angstrom
        )
        molecule_diff_coords.assign_fractional_bond_orders(
            toolkit_registry=toolkit_wrapper,
            use_conformers=molecule_diff_coords.conformers,
            bond_order_model="am1-wiberg",
        )

        for bond1, bond2 in zip(molecule.bonds, molecule_diff_coords.bonds):
            assert abs(bond1.fractional_bond_order - bond2.fractional_bond_order) > 1e-3

    @pytest.mark.parametrize(
        "bond_order_model",
        ["am1-wiberg", "am1-wiberg-elf10", "pm3-wiberg", "pm3-wiberg-elf10"],
    )
    def test_assign_fractional_bond_orders_neutral_charge_mol(self, bond_order_model):
        """Test OpenEyeToolkitWrapper assign_fractional_bond_orders() for neutral and
        charged molecule"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        # Reading neutral molecule from file
        filename = get_data_file_path("molecules/CID20742535_neutral.sdf")
        molecule1 = Molecule.from_file(filename)
        # Reading negative molecule from file
        filename = get_data_file_path("molecules/CID20742535_anion.sdf")
        molecule2 = Molecule.from_file(filename)

        # Checking that only one additional bond is present in the neutral molecule
        assert len(molecule1.bonds) == len(molecule2.bonds) + 1

        molecule1.assign_fractional_bond_orders(
            toolkit_registry=toolkit_wrapper,
            bond_order_model=bond_order_model,
            use_conformers=molecule1.conformers,
        )

        for i in molecule1.bonds:
            if i.is_aromatic:
                # Checking aromatic bonds
                assert 1.05 < i.fractional_bond_order < 1.65
            elif i.atom1.atomic_number == 1 or i.atom2.atomic_number == 1:
                # Checking bond order of C-H or O-H bonds are around 1
                assert 0.85 < i.fractional_bond_order < 1.05
            elif i.atom1.atomic_number == 8 or i.atom2.atomic_number == 8:
                # Checking C-O single bond
                wbo_C_O_neutral = i.fractional_bond_order
                assert 1.0 < wbo_C_O_neutral < 1.5
            else:
                # Should be C-C single bond
                assert (i.atom1_index == 4 and i.atom2_index == 6) or (
                    i.atom1_index == 6 and i.atom2_index == 4
                )
                wbo_C_C_neutral = i.fractional_bond_order
                assert 1.0 < wbo_C_C_neutral < 1.3

        molecule2.assign_fractional_bond_orders(
            toolkit_registry=toolkit_wrapper,
            bond_order_model=bond_order_model,
            use_conformers=molecule2.conformers,
        )
        for i in molecule2.bonds:
            if i.is_aromatic:
                # Checking aromatic bonds
                assert 1.05 < i.fractional_bond_order < 1.65
            elif i.atom1.atomic_number == 1 or i.atom2.atomic_number == 1:
                # Checking bond order of C-H or O-H bonds are around 1
                assert 0.85 < i.fractional_bond_order < 1.05
            elif i.atom1.atomic_number == 8 or i.atom2.atomic_number == 8:
                # Checking C-O single bond
                wbo_C_O_anion = i.fractional_bond_order
                assert 1.3 < wbo_C_O_anion < 1.8
            else:
                # Should be C-C single bond
                assert (i.atom1_index == 4 and i.atom2_index == 6) or (
                    i.atom1_index == 6 and i.atom2_index == 4
                )
                wbo_C_C_anion = i.fractional_bond_order
                assert 1.0 < wbo_C_C_anion < 1.3

        # Wiberg bond order of C-C single bond is higher in the anion
        assert wbo_C_C_anion > wbo_C_C_neutral
        # Wiberg bond order of C-O bond is higher in the anion
        assert wbo_C_O_anion > wbo_C_O_neutral

    def test_assign_fractional_bond_orders_invalid_method(self):
        """
        Test that OpenEyeToolkitWrapper assign_fractional_bond_orders() raises the
        correct error if an invalid charge model is provided
        """
        toolkit_wrapper = OpenEyeToolkitWrapper()

        molecule = toolkit_wrapper.from_smiles("C")

        expected_error = (
            "Bond order model 'not a real bond order model' is not supported by "
            "OpenEyeToolkitWrapper. Supported models are "
            r"\['am1-wiberg', 'am1-wiberg-elf10', 'pm3-wiberg', 'pm3-wiberg-elf10'\]"
        )

        with pytest.raises(ValueError, match=expected_error):
            molecule.assign_fractional_bond_orders(
                toolkit_registry=toolkit_wrapper,
                bond_order_model="not a real bond order model",
            )

    @pytest.mark.slow
    def test_max_substructure_matches_can_handle_large_molecule(self):
        """Test OpenEyeToolkitWrapper substructure search handles more than the default of MaxMatches = 1024
        See https://github.com/openforcefield/openff-toolkit/pull/509 .
        """

        tk = OpenEyeToolkitWrapper()
        smiles = "C" * 600
        molecule = tk.from_smiles(smiles)
        query = "[C:1]~[C:2]"
        ret = molecule.chemical_environment_matches(query, toolkit_registry=tk)
        assert len(ret) == 1198
        assert len(ret[0]) == 2

    def test_find_rotatable_bonds(self):
        """Test finding rotatable bonds while ignoring some groups"""

        # test a simple molecule
        ethanol = create_ethanol()
        bonds = ethanol.find_rotatable_bonds()
        assert len(bonds) == 2
        for bond in bonds:
            assert ethanol.atoms[bond.atom1_index].atomic_number != 1
            assert ethanol.atoms[bond.atom2_index].atomic_number != 1

        # now ignore the C-O bond, forwards
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups="[#6:1]-[#8:2]")
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 6

        # now ignore the O-C bond, backwards
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups="[#8:1]-[#6:2]")
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 6

        # now ignore the C-C bond
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups="[#6:1]-[#6:2]")
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 8

        # ignore a list of searches, forward
        bonds = ethanol.find_rotatable_bonds(
            ignore_functional_groups=["[#6:1]-[#8:2]", "[#6:1]-[#6:2]"]
        )
        assert bonds == []

        # ignore a list of searches, backwards
        bonds = ethanol.find_rotatable_bonds(
            ignore_functional_groups=["[#6:1]-[#6:2]", "[#8:1]-[#6:2]"]
        )
        assert bonds == []

        # test  molecules that should have no rotatable bonds
        cyclohexane = create_cyclohexane()
        bonds = cyclohexane.find_rotatable_bonds()
        assert bonds == []

        methane = Molecule.from_smiles("C")
        bonds = methane.find_rotatable_bonds()
        assert bonds == []

        ethene = Molecule.from_smiles("C=C")
        bonds = ethene.find_rotatable_bonds()
        assert bonds == []

        terminal_forwards = "[*]~[*:1]-[X2H1,X3H2,X4H3:2]-[#1]"
        terminal_backwards = "[#1]-[X2H1,X3H2,X4H3:1]-[*:2]~[*]"
        # test removing terminal rotors
        toluene = Molecule.from_file(get_data_file_path("molecules/toluene.sdf"))
        bonds = toluene.find_rotatable_bonds()
        assert len(bonds) == 1
        assert toluene.atoms[bonds[0].atom1_index].atomic_number == 6
        assert toluene.atoms[bonds[0].atom2_index].atomic_number == 6

        # find terminal bonds forward
        bonds = toluene.find_rotatable_bonds(ignore_functional_groups=terminal_forwards)
        assert bonds == []

        # find terminal bonds backwards
        bonds = toluene.find_rotatable_bonds(
            ignore_functional_groups=terminal_backwards
        )
        assert bonds == []

        # TODO: Check partial charge invariants (total charge, charge equivalence)

        # TODO: Add test for aromaticity
        # TODO: Add test and molecule functionality for isotopes

    @pytest.mark.parametrize(
        ("smiles", "n_atom_rings", "n_bond_rings"),
        [
            ("c1ccc2ccccc2c1", 10, 11),
            ("c1ccc(cc1)c2ccccc2", 12, 12),
            ("Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C", 30, 30),
        ],
    )
    def test_is_in_ring(self, smiles, n_atom_rings, n_bond_rings):
        """Test Atom.is_in_ring and Bond.is_in_ring"""
        mol = Molecule.from_smiles(smiles)

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper()])

        atoms_in_ring = [
            atom
            for atom in mol.atoms
            if atom.is_in_ring(toolkit_registry=toolkit_registry)
        ]

        bonds_in_ring = [
            bond
            for bond in mol.bonds
            if bond.is_in_ring(toolkit_registry=toolkit_registry)
        ]

        assert len(atoms_in_ring) == n_atom_rings
        assert len(bonds_in_ring) == n_bond_rings

    def test_central_biphenyl_bond(self):
        """Test that `Bond.is_in_ring` is False for the central bond in a phenyl"""
        # Use a mapped smiles to ensure atom order while looking up central bond
        # Generated via Molecule.from_smiles("c1ccc(cc1)c2ccccc2").to_smiles(mapped=True)

        biphenyl = Molecule.from_mapped_smiles(
            "[H:13][c:1]1[c:2]([c:3]([c:4]([c:5]([c:6]1[H:17])[H:16])[c:7]2[c:8]([c:9]([c:10]"
            "([c:11]([c:12]2[H:22])[H:21])[H:20])[H:19])[H:18])[H:15])[H:14]"
        )

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])

        this_bond = biphenyl.get_bond_between(3, 6)
        assert this_bond.is_in_ring(toolkit_registry=toolkit_registry) is False

        ring_bonds = [
            b for b in biphenyl.bonds if b.is_in_ring(toolkit_registry=toolkit_registry)
        ]
        assert len(ring_bonds) == 12

    def test_unattached_is_in_ring(self):
        toolkit = OpenEyeToolkitWrapper()
        dummy_atom = Atom(1, 0, False)

        with pytest.raises(NotAttachedToMoleculeError, match="Atom"):
            toolkit.atom_is_in_ring(dummy_atom)

        # The Bond constructor checks to see that the atoms are in a molecule,
        # so not so straightforward to make one from the constructor
        dummy_bond = Molecule.from_smiles("O").bonds[0]
        dummy_bond._molecule = None

        with pytest.raises(NotAttachedToMoleculeError, match="Bond"):
            toolkit.bond_is_in_ring(dummy_bond)

    def test_find_matches_unique(self):
        """Test the expected behavior of the `unique` argument in find_matches"""
        smirks = "[C:1]~[C:2]~[C:3]"
        tk = OpenEyeToolkitWrapper()

        mol = Molecule.from_smiles("CCC")

        assert len(tk.find_smarts_matches(mol, smirks, unique=True)) == 1
        assert len(tk.find_smarts_matches(mol, smirks, unique=False)) == 2


@requires_rdkit
class TestRDKitToolkitWrapper:
    """Test the RDKitToolkitWrapper"""

    def test_smiles(self):
        """Test RDKitToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        # This differs from OE's expected output due to different canonicalization schemes
        smiles = "[H][C]([H])([H])[C]([H])([H])[H]"
        molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)
        # When making a molecule from SMILES, partial charges should be initialized to None
        assert molecule.partial_charges is None
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        # print(smiles, smiles2)
        assert smiles == smiles2

    @pytest.mark.parametrize(
        "smiles,exception_regex",
        [
            (r"C\C(F)=C(/F)CC(C)(Cl)Br", "Undefined chiral centers"),
            (r"C\C(F)=C(/F)C[C@@](C)(Cl)Br", None),
            (r"CC(F)=C(F)C[C@@](C)(Cl)Br", "Bonds with undefined stereochemistry"),
        ],
    )
    def test_smiles_missing_stereochemistry(self, smiles, exception_regex):
        """Test RDKitToolkitWrapper to_smiles() and from_smiles() when given ambiguous stereochemistry"""
        toolkit_wrapper = RDKitToolkitWrapper()

        if exception_regex is not None:
            with pytest.raises(UndefinedStereochemistryError, match=exception_regex):
                Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)
            Molecule.from_smiles(
                smiles, toolkit_registry=toolkit_wrapper, allow_undefined_stereo=True
            )
        else:
            Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)

    # TODO: test_smiles_round_trip

    def test_smiles_add_H(self):
        """Test RDKitToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        input_smiles = "CC"
        # This differs from OE's expected output due to different canonicalization schemes
        expected_output_smiles = "[H][C]([H])([H])[C]([H])([H])[H]"
        molecule = Molecule.from_smiles(input_smiles, toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles2 == expected_output_smiles

    def test_rdkit_from_smiles_radical(self):
        """Test that parsing an SMILES with a radical raises RadicalsNotSupportedError."""
        with pytest.raises(RadicalsNotSupportedError):
            RDKitToolkitWrapper().from_smiles("[CH3]")

    def test_rdkit_from_smiles_transition_metal_radical(self):
        """Test that parsing an SMILES with a transition metal radical works."""
        RDKitToolkitWrapper().from_smiles("[Zn+2]")

    def test_rdkit_from_smiles_hydrogens_are_explicit(self):
        """
        Test to ensure that RDKitToolkitWrapper.from_smiles has the proper behavior with
        respect to its hydrogens_are_explicit kwarg
        """
        toolkit_wrapper = RDKitToolkitWrapper()
        smiles_impl = "C#C"
        with pytest.raises(
            ValueError,
            match="but RDKit toolkit interpreted SMILES 'C#C' as having implicit hydrogen",
        ):
            offmol = Molecule.from_smiles(
                smiles_impl,
                toolkit_registry=toolkit_wrapper,
                hydrogens_are_explicit=True,
            )
        offmol = Molecule.from_smiles(
            smiles_impl, toolkit_registry=toolkit_wrapper, hydrogens_are_explicit=False
        )
        assert offmol.n_atoms == 4

        smiles_expl = "[H][C]#[C][H]"
        offmol = Molecule.from_smiles(
            smiles_expl, toolkit_registry=toolkit_wrapper, hydrogens_are_explicit=True
        )
        assert offmol.n_atoms == 4
        # It's debatable whether this next function should pass. Strictly speaking, the hydrogens in this SMILES
        # _are_ explicit, so allowing "hydrogens_are_explicit=False" through here is allowing a contradiction.
        # We might rethink the name of this kwarg.

        offmol = Molecule.from_smiles(
            smiles_expl, toolkit_registry=toolkit_wrapper, hydrogens_are_explicit=False
        )
        assert offmol.n_atoms == 4

    @pytest.mark.parametrize("molecule", get_mini_drug_bank(RDKitToolkitWrapper))
    def test_to_inchi(self, molecule):
        """Test, but do not validate, conversion to standard and non-standard InChI"""

        toolkit = RDKitToolkitWrapper()
        molecule.to_inchi(toolkit_registry=toolkit)
        molecule.to_inchi(fixed_hydrogens=True, toolkit_registry=toolkit)

    @pytest.mark.parametrize("molecule", get_mini_drug_bank(RDKitToolkitWrapper))
    def test_to_inchikey(self, molecule):
        """Test, but do not validate, the conversion to standard and non-standard InChIKey"""

        toolkit = RDKitToolkitWrapper()
        molecule.to_inchikey(toolkit_registry=toolkit)
        molecule.to_inchikey(fixed_hydrogens=True, toolkit_registry=toolkit)

    def test_from_bad_inchi(self):
        """Test building a molecule from a bad InChI string"""

        toolkit = RDKitToolkitWrapper()
        inchi = "InChI=1S/ksbfksfksfksbfks"
        with pytest.raises(RuntimeError):
            Molecule.from_inchi(inchi, toolkit_registry=toolkit)

    inchi_data = [
        {
            "molecule": create_ethanol(),
            "standard_inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            "fixed_hydrogen_inchi": "InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3",
        },
        {
            "molecule": create_reversed_ethanol(),
            "standard_inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            "fixed_hydrogen_inchi": "InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3",
        },
        {
            "molecule": create_acetaldehyde(),
            "standard_inchi": "InChI=1S/C2H4O/c1-2-3/h2H,1H3",
            "fixed_hydrogen_inchi": "InChI=1/C2H4O/c1-2-3/h2H,1H3",
        },
        {
            "molecule": create_cyclohexane(),
            "standard_inchi": "InChI=1S/C6H12/c1-2-4-6-5-3-1/h1-6H2",
            "fixed_hydrogen_inchi": "InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2",
        },
    ]

    @pytest.mark.parametrize("data", inchi_data)
    def test_from_inchi(self, data):
        """Test building a molecule from standard and non-standard InChI strings."""

        toolkit = RDKitToolkitWrapper()

        ref_mol = data["molecule"]
        # make a molecule from inchi
        inchi_mol = Molecule.from_inchi(
            data["standard_inchi"], toolkit_registry=toolkit
        )
        assert inchi_mol.to_inchi(toolkit_registry=toolkit) == data["standard_inchi"]

        def compare_mols(ref_mol, inchi_mol):
            assert ref_mol.n_atoms == inchi_mol.n_atoms
            assert ref_mol.n_bonds == inchi_mol.n_bonds
            assert ref_mol.n_angles == inchi_mol.n_angles
            assert ref_mol.n_propers == inchi_mol.n_propers
            assert ref_mol.is_isomorphic_with(inchi_mol) is True

        compare_mols(ref_mol, inchi_mol)

        # now make the molecule from the non-standard inchi and compare
        nonstandard_inchi_mol = Molecule.from_inchi(data["fixed_hydrogen_inchi"])
        assert (
            nonstandard_inchi_mol.to_inchi(
                fixed_hydrogens=True, toolkit_registry=toolkit
            )
            == data["fixed_hydrogen_inchi"]
        )

        compare_mols(ref_mol, nonstandard_inchi_mol)

    @pytest.mark.parametrize("molecule", get_mini_drug_bank(RDKitToolkitWrapper))
    def test_non_standard_inchi_round_trip(self, molecule):
        """Test if a molecule can survive an InChi round trip test in some cases the standard InChI
        will not be enough to ensure information is preserved so we test the non-standard inchi here.
        """

        from openff.toolkit.utils.toolkits import UndefinedStereochemistryError

        toolkit = RDKitToolkitWrapper()
        inchi = molecule.to_inchi(fixed_hydrogens=True, toolkit_registry=toolkit)
        # make a copy of the molecule from the inchi string
        if molecule.name in rdkit_inchi_stereochemistry_lost:
            # some molecules lose stereochemsitry so they are skipped
            # if we fail here the molecule may of been fixed
            with pytest.raises(UndefinedStereochemistryError):
                mol2 = molecule.from_inchi(inchi, toolkit_registry=toolkit)

        else:
            print(molecule.name)
            mol2 = molecule.from_inchi(inchi, toolkit_registry=toolkit)

            # compare the full molecule excluding the properties dictionary
            # turn of the bond order matching as this could move in the aromatic rings
            assert molecule.is_isomorphic_with(
                mol2, bond_order_matching=False, toolkit_registry=toolkit
            )

    def test_smiles_charged(self):
        """Test RDKitWrapper functions for reading/writing charged SMILES"""
        toolkit_wrapper = RDKitToolkitWrapper()
        # This differs from OE's expected output due to different canonicalization schemes
        smiles = "[H][C]([H])([H])[N+]([H])([H])[H]"
        molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2

    def test_to_from_rdkit_core_props_filled(self):
        """Test RDKitToolkitWrapper to_rdkit() and from_rdkit() when given populated core property fields"""
        toolkit_wrapper = RDKitToolkitWrapper()

        # Replacing with a simple molecule with stereochemistry
        input_smiles = r"C\C(F)=C(/F)C[C@@](C)(Cl)Br"
        expected_output_smiles = r"[H][C]([H])([H])/[C]([F])=[C](\[F])[C]([H])([H])[C@@]([Cl])([Br])[C]([H])([H])[H]"
        molecule = Molecule.from_smiles(input_smiles, toolkit_registry=toolkit_wrapper)
        assert (
            molecule.to_smiles(toolkit_registry=toolkit_wrapper)
            == expected_output_smiles
        )

        # Populate core molecule property fields
        molecule.name = "Alice"
        partial_charges = unit.Quantity(
            np.array(
                [
                    -0.9,
                    -0.8,
                    -0.7,
                    -0.6,
                    -0.5,
                    -0.4,
                    -0.3,
                    -0.2,
                    -0.1,
                    0.0,
                    0.1,
                    0.2,
                    0.3,
                    0.4,
                    0.5,
                    0.6,
                    0.7,
                    0.8,
                ]
            ),
            unit.elementary_charge,
        )
        molecule.partial_charges = partial_charges
        coords = unit.Quantity(
            np.array(
                [
                    ["0.0", "1.0", "2.0"],
                    ["3.0", "4.0", "5.0"],
                    ["6.0", "7.0", "8.0"],
                    ["9.0", "10.0", "11.0"],
                    ["12.0", "13.0", "14.0"],
                    ["15.0", "16.0", "17.0"],
                    ["18.0", "19.0", "20.0"],
                    ["21.0", "22.0", "23.0"],
                    ["24.0", "25.0", "26.0"],
                    ["27.0", "28.0", "29.0"],
                    ["30.0", "31.0", "32.0"],
                    ["33.0", "34.0", "35.0"],
                    ["36.0", "37.0", "38.0"],
                    ["39.0", "40.0", "41.0"],
                    ["42.0", "43.0", "44.0"],
                    ["45.0", "46.0", "47.0"],
                    ["48.0", "49.0", "50.0"],
                    ["51.0", "52.0", "53.0"],
                ]
            ),
            unit.angstrom,
        )
        molecule.add_conformer(coords)
        # Populate core atom property fields
        molecule.atoms[2].name = "Bob"
        # Ensure one atom has its stereochemistry specified
        central_carbon_stereo_specified = False
        for atom in molecule.atoms:
            if (atom.atomic_number == 6) and atom.stereochemistry == "S":
                central_carbon_stereo_specified = True
        assert central_carbon_stereo_specified

        # Populate bond core property fields
        fractional_bond_orders = [float(val) for val in range(18)]
        for fbo, bond in zip(fractional_bond_orders, molecule.bonds):
            bond.fractional_bond_order = fbo

        # Do a first conversion to/from rdmol
        rdmol = molecule.to_rdkit()
        molecule2 = Molecule.from_rdkit(rdmol)

        # Test that properties survived first conversion
        # assert molecule.to_dict() == molecule2.to_dict()
        assert molecule.name == molecule2.name
        # NOTE: This expects the same indexing scheme in the original and new molecule

        central_carbon_stereo_specified = False
        for atom in molecule2.atoms:
            if (atom.atomic_number == 6) and atom.stereochemistry == "S":
                central_carbon_stereo_specified = True
        assert central_carbon_stereo_specified
        for atom1, atom2 in zip(molecule.atoms, molecule2.atoms):
            assert atom1.to_dict() == atom2.to_dict()
        for bond1, bond2 in zip(molecule.bonds, molecule2.bonds):
            assert bond1.to_dict() == bond2.to_dict()
        assert (molecule.conformers[0] == molecule2.conformers[0]).all()
        for pc1, pc2 in zip(molecule._partial_charges, molecule2._partial_charges):
            pc1_ul = pc1.m_as(unit.elementary_charge)
            pc2_ul = pc2.m_as(unit.elementary_charge)
            assert_almost_equal(pc1_ul, pc2_ul, decimal=6)
        assert (
            molecule2.to_smiles(toolkit_registry=toolkit_wrapper)
            == expected_output_smiles
        )
        # TODO: This should be its own test

    def test_to_from_rdkit_core_props_unset(self):
        """Test RDKitToolkitWrapper to_rdkit() and from_rdkit() when given empty core property fields"""
        toolkit_wrapper = RDKitToolkitWrapper()

        # Replacing with a simple molecule with stereochemistry
        input_smiles = r"C\C(F)=C(/F)C[C@](C)(Cl)Br"
        expected_output_smiles = r"[H][C]([H])([H])/[C]([F])=[C](\[F])[C]([H])([H])[C@]([Cl])([Br])[C]([H])([H])[H]"
        molecule = Molecule.from_smiles(input_smiles, toolkit_registry=toolkit_wrapper)
        assert (
            molecule.to_smiles(toolkit_registry=toolkit_wrapper)
            == expected_output_smiles
        )

        # Ensure one atom has its stereochemistry specified
        central_carbon_stereo_specified = False
        for atom in molecule.atoms:
            if (atom.atomic_number == 6) and atom.stereochemistry == "R":
                central_carbon_stereo_specified = True
        assert central_carbon_stereo_specified

        # Do a first conversion to/from rdmol
        rdmol = molecule.to_rdkit()
        molecule2 = Molecule.from_rdkit(rdmol)

        # Test that properties survived first conversion
        assert molecule.name == molecule2.name
        # NOTE: This expects the same indexing scheme in the original and new molecule

        central_carbon_stereo_specified = False
        for atom in molecule2.atoms:
            if (atom.atomic_number == 6) and atom.stereochemistry == "R":
                central_carbon_stereo_specified = True
        assert central_carbon_stereo_specified
        for atom1, atom2 in zip(molecule.atoms, molecule2.atoms):
            assert atom1.to_dict() == atom2.to_dict()
        for bond1, bond2 in zip(molecule.bonds, molecule2.bonds):
            assert bond1.to_dict() == bond2.to_dict()
        # The molecule was initialized from SMILES, so mol.conformers arrays should be None for both
        assert molecule.conformers is None
        assert molecule2.conformers is None
        # The molecule was initialized from SMILES, so mol.partial_charges arrays should be None for both
        assert molecule.partial_charges is None
        assert molecule2.partial_charges is None

        assert (
            molecule2.to_smiles(toolkit_registry=toolkit_wrapper)
            == expected_output_smiles
        )

    def test_to_from_rdkit_hierarchy_metadata(self):
        """
        Test roundtripping to/from ``OpenEyeToolkitWrapper`` for molecules with PDB hierarchy metadata
        """
        for molecule in topology_with_metadata().molecules:
            rdmol = molecule.to_rdkit()
            roundtrip_mol = Molecule.from_rdkit(rdmol)

            # Check RDMol
            for orig_atom, rd_atom in zip(molecule.atoms, rdmol.GetAtoms()):
                atom_has_any_metadata = (
                    ("residue_name" in orig_atom.metadata)
                    or ("residue_number" in orig_atom.metadata)
                    or ("insertion_code" in orig_atom.metadata)
                    or ("chain_id" in orig_atom.metadata)
                )

                if not (atom_has_any_metadata):
                    assert rd_atom.GetPDBResidueInfo() is None
                    continue

                if "residue_name" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["residue_name"]
                        == rd_atom.GetPDBResidueInfo().GetResidueName()
                    )
                else:
                    assert rd_atom.GetPDBResidueInfo().GetResidueName() == ""

                if "residue_number" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["residue_number"]
                        == rd_atom.GetPDBResidueInfo().GetResidueNumber()
                    )
                else:
                    assert rd_atom.GetPDBResidueInfo().GetResidueNumber() == 0

                if "insertion_code" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["insertion_code"]
                        == rd_atom.GetPDBResidueInfo().GetInsertionCode()
                    )
                else:
                    assert rd_atom.GetPDBResidueInfo().GetInsertionCode() == " "

                if "chain_id" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["chain_id"]
                        == rd_atom.GetPDBResidueInfo().GetChainId()
                    )
                else:
                    assert rd_atom.GetPDBResidueInfo().GetChainId() == ""

            # Check roundtripped OFFMol
            for orig_atom, roundtrip_atom in zip(molecule.atoms, roundtrip_mol.atoms):
                atom_has_any_metadata = (
                    ("residue_name" in orig_atom.metadata)
                    or ("residue_number" in orig_atom.metadata)
                    or ("insertion_code" in orig_atom.metadata)
                    or ("chain_id" in orig_atom.metadata)
                )
                if not (atom_has_any_metadata):
                    assert roundtrip_atom.metadata == {}
                    continue

                if "residue_name" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["residue_name"]
                        == roundtrip_atom.metadata["residue_name"]
                    )
                else:
                    assert roundtrip_atom.metadata["residue_name"] == ""

                if "residue_number" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["residue_number"]
                        == roundtrip_atom.metadata["residue_number"]
                    )
                else:
                    assert roundtrip_atom.metadata["residue_number"] == 0

                if "insertion_code" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["insertion_code"]
                        == roundtrip_atom.metadata["insertion_code"]
                    )
                else:
                    assert roundtrip_atom.metadata["insertion_code"] == 0

                if "chain_id" in orig_atom.metadata:
                    assert (
                        orig_atom.metadata["chain_id"]
                        == roundtrip_atom.metadata["chain_id"]
                    )
                else:
                    assert roundtrip_atom.metadata["chain_id"] == ""

    def test_from_rdkit_implicit_hydrogens(self):
        """
        Test that hydrogens are inferred from hydrogen-less RDKit molecules,
        unless the option is turned off.
        """
        from rdkit import Chem

        rdmol = Chem.MolFromSmiles("CC")
        offmol = Molecule.from_rdkit(rdmol)

        assert any([a.atomic_number == 1 for a in offmol.atoms])

        offmol_no_h = Molecule.from_rdkit(rdmol, hydrogens_are_explicit=True)
        assert not any([a.atomic_number == 1 for a in offmol_no_h.atoms])

    def test_from_rdkit_radical(self):
        """Test that parsing an rdmol with a radical raises RadicalsNotSupportedError."""
        from rdkit import Chem

        rdmol = Chem.MolFromSmiles("[CH3]")

        with pytest.raises(RadicalsNotSupportedError):
            RDKitToolkitWrapper().from_rdkit(rdmol)

    def test_from_rdkit_transition_metal_radical(self):
        """Test that parsing an rdmol with a transition metal radical works."""
        from rdkit import Chem

        rdmol = Chem.MolFromSmiles("[Zn+2]")

        RDKitToolkitWrapper().from_rdkit(rdmol)

    @pytest.mark.parametrize(
        "smiles, expected_map", [("[Cl:1][Cl]", {0: 1}), ("[Cl:1][Cl:2]", {0: 1, 1: 2})]
    )
    def test_from_rdkit_atom_map(self, smiles, expected_map):
        """
        Test OpenEyeToolkitWrapper for loading a molecule with implicit
        hydrogens (correct behavior is to add them explicitly)
        """
        from rdkit import Chem

        off_molecule = Molecule.from_rdkit(Chem.MolFromSmiles(smiles))
        assert off_molecule.properties["atom_map"] == expected_map

    def test_from_pathlib_path(self):
        ethanol = create_ethanol()
        with NamedTemporaryFile(suffix=".sdf") as outfile:
            filename = str(outfile.name)
            ethanol.to_file(filename, file_format="sdf")

            toolkit = RDKitToolkitWrapper()
            toolkit.from_file(pathlib.Path(filename), file_format="sdf")

    def test_file_extension_case(self):
        """
        Test round-trips of some file extensions when called directly from the toolkit wrappers,
        including lower- and uppercase file extensions. Note that this test does not ensure
        accuracy, it only tests that reading/writing without raising an exception.
        """
        mols_in = RDKitToolkitWrapper().from_file(
            file_path=get_data_file_path("molecules/ethanol.sdf"), file_format="sdf"
        )

        assert len(mols_in) > 0

        mols_in = RDKitToolkitWrapper().from_file(
            file_path=get_data_file_path("molecules/ethanol.sdf"), file_format="SDF"
        )

        assert len(mols_in) > 0

    def test_get_sdf_coordinates(self):
        """Test RDKitToolkitWrapper for importing a single set of coordinates from a sdf file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path("molecules/toluene.sdf")
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule.conformers) == 1
        assert molecule.conformers[0].shape == (15, 3)
        assert_almost_equal(
            molecule.conformers[0][5][1].m_as(unit.angstrom), 22.9800, decimal=4
        )

    def test_read_sdf_charges(self):
        """Test RDKitToolkitWrapper for importing a charges from a sdf file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path("molecules/ethanol_partial_charges.sdf")
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert molecule.partial_charges is not None
        assert molecule.partial_charges[0] == -0.4 * unit.elementary_charge
        assert molecule.partial_charges[-1] == 0.4 * unit.elementary_charge

    def test_write_sdf_charges(self):
        """Test RDKitToolkitWrapper for writing partial charges to a sdf file"""
        from io import StringIO

        toolkit_wrapper = RDKitToolkitWrapper()
        ethanol = create_ethanol()
        sio = StringIO()
        ethanol.to_file(sio, "SDF", toolkit_registry=toolkit_wrapper)
        sdf_text = sio.getvalue()
        # The output lines of interest here will look like
        # >  <atom.dprop.PartialCharge>  (1)
        # -0.40000000000000002 -0.29999999999999999 -0.20000000000000001 -0.10000000000000001 0.01 0.10000000000000001 0.20000000000000001 0.29999999999999999 0.40000000000000002

        # Parse the SDF text, grabbing the numeric line above
        sdf_split = sdf_text.split("\n")
        charge_line_found = False
        for line in sdf_split:
            if charge_line_found:
                charges = [float(i) for i in line.split()]
                break
            if ">  <atom.dprop.PartialCharge>" in line:
                charge_line_found = True

        # Make sure that a charge line was ever found
        assert charge_line_found

        # Make sure that the charges found were correct
        assert_almost_equal(
            charges, [-0.4, -0.3, -0.2, -0.1, 0.00001, 0.1, 0.2, 0.3, 0.4]
        )

    def test_sdf_properties_roundtrip(self):
        """Test RDKitToolkitWrapper for performing a round trip of a molecule with defined partial charges
        and entries in the properties dict to and from a sdf file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        ethanol = create_ethanol()
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix=".sdf") as iofile:
            ethanol.to_file(
                iofile.name, file_format="SDF", toolkit_registry=toolkit_wrapper
            )
            ethanol2 = Molecule.from_file(
                iofile.name, file_format="SDF", toolkit_registry=toolkit_wrapper
            )
        assert (ethanol.partial_charges == ethanol2.partial_charges).all()

        # Now test with no properties or charges
        ethanol = create_ethanol()
        ethanol.partial_charges = None
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix=".sdf") as iofile:
            ethanol.to_file(
                iofile.name, file_format="SDF", toolkit_registry=toolkit_wrapper
            )
            ethanol2 = Molecule.from_file(
                iofile.name, file_format="SDF", toolkit_registry=toolkit_wrapper
            )
        assert ethanol2.partial_charges is None
        assert ethanol2.properties == {}

    def test_write_sdf_no_charges(self):
        """Test RDKitToolkitWrapper for writing an SDF file with no charges"""
        from io import StringIO

        toolkit_wrapper = RDKitToolkitWrapper()
        ethanol = create_ethanol()
        ethanol.partial_charges = None
        sio = StringIO()
        ethanol.to_file(sio, "SDF", toolkit_registry=toolkit_wrapper)
        sdf_text = sio.getvalue()
        # In our current configuration, if the OFFMol doesn't have partial charges, we DO NOT want a partial charge
        # block to be written. For reference, it's possible to indicate that a partial charge is not known by writing
        # out "n/a" (or another placeholder) in the partial charge block atoms without charges.
        assert ">  <atom.dprop.PartialCharge>" not in sdf_text

    def test_read_ethene_sdf(self):
        """
        Test that RDKitToolkitWrapper can load an ethene molecule without complaining about bond stereo.
        See https://github.com/openforcefield/openff-toolkit/issues/785
        """
        ethene_file_path = get_data_file_path("molecules/ethene_rdkit.sdf")
        toolkit_wrapper = RDKitToolkitWrapper()
        toolkit_wrapper.from_file(ethene_file_path, file_format="sdf")

    def test_load_multiconformer_sdf_as_separate_molecules(self):
        """
        Test RDKitToolkitWrapper for reading a "multiconformer" SDF, which the OFF
        Toolkit should treat as separate molecules
        """
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path("molecules/methane_multiconformer.sdf")
        molecules = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecules) == 2
        assert len(molecules[0].conformers) == 1
        assert len(molecules[1].conformers) == 1
        assert molecules[0].conformers[0].shape == (5, 3)

    def test_load_multiconformer_sdf_as_separate_molecules_properties(self):
        """
        Test RDKitToolkitWrapper for reading a "multiconformer" SDF, which the OFF
        Toolkit should treat as separate molecules
        """
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path("molecules/methane_multiconformer_properties.sdf")
        molecules = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecules) == 2
        assert len(molecules[0].conformers) == 1
        assert len(molecules[1].conformers) == 1
        assert molecules[0].conformers[0].shape == (5, 3)
        # The first molecule in the SDF has the following properties and charges:
        assert molecules[0].properties["test_property_key"] == "test_property_value"
        np.testing.assert_allclose(
            molecules[0].partial_charges.m_as(unit.elementary_charge),
            [-0.108680, 0.027170, 0.027170, 0.027170, 0.027170],
        )
        # The second molecule in the SDF has the following properties and charges:
        assert molecules[1].properties["test_property_key"] == "test_property_value2"
        assert (
            molecules[1].properties["another_test_property_key"]
            == "another_test_property_value"
        )
        np.testing.assert_allclose(
            molecules[1].partial_charges.m_as(unit.elementary_charge),
            [0.027170, 0.027170, 0.027170, 0.027170, -0.108680],
        )

    def test_write_multiconformer_mol_as_sdf(self):
        """
        Test RDKitToolkitWrapper for writing a multiconformer molecule to SDF. The OFF toolkit should only
        save the first conformer
        """
        from io import StringIO

        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path("molecules/ethanol.sdf")
        ethanol = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        ethanol.partial_charges = unit.Quantity(
            np.array([-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]),
            unit.elementary_charge,
        )
        ethanol.properties["test_prop"] = "test_value"
        new_conf = ethanol.conformers[0] + (
            np.ones(ethanol.conformers[0].shape) * unit.angstrom
        )
        ethanol.add_conformer(new_conf)
        sio = StringIO()
        ethanol.to_file(sio, "sdf", toolkit_registry=toolkit_wrapper)
        data = sio.getvalue()
        # In SD format, each molecule ends with "$$$$"
        assert data.count("$$$$") == 1
        # A basic SDF for ethanol would be 27 lines, though the properties add three more
        assert len(data.split("\n")) == 30
        assert "test_prop" in data
        assert "<atom.dprop.PartialCharge>" in data
        # Ensure the first conformer's first atom's X coordinate is in the file
        assert str(ethanol.conformers[0][0][0].m_as(unit.angstrom))[:5] in data
        # Ensure the SECOND conformer's first atom's X coordinate is NOT in the file
        assert str(ethanol.conformers[1][0][0].m_as(unit.angstrom))[:5] not in data

    def test_write_multiconformer_pdb(self):
        """
        Make sure RDKit can write multi conformer PDB files.
        """
        from io import StringIO

        toolkit = RDKitToolkitWrapper()
        # load up a multiconformer pdb file and condense down the conformers
        molecules = Molecule.from_file(
            get_data_file_path("molecules/butane_multi.sdf"), toolkit_registry=toolkit
        )
        butane = molecules.pop(0)
        for mol in molecules:
            butane.add_conformer(mol.conformers[0])
        assert butane.n_conformers == 7
        sio = StringIO()
        butane.to_file(sio, "pdb", toolkit_registry=toolkit)
        # we need to make sure each conformer is wrote to the file
        pdb = sio.getvalue()
        for i in range(1, 8):
            assert f"MODEL        {i}" in pdb

    # Unskip this when we implement PDB-reading support for RDKitToolkitWrapper
    @pytest.mark.skip
    def test_get_pdb_coordinates(self):
        """Test RDKitToolkitWrapper for importing a single set of coordinates from a pdb file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path("molecules/toluene.pdb")
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule.conformers) == 1
        assert molecule.conformers[0].shape == (15, 3)

    # Unskip this when we implement PDB-reading support for RDKitToolkitWrapper
    @pytest.mark.skip
    def test_load_aromatic_pdb(self):
        """Test OpenEyeToolkitWrapper for importing molecule conformers"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path("molecules/toluene.pdb")
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule.conformers) == 1
        assert molecule.conformers[0].shape == (15, 3)

    def test_generate_conformers(self):
        """Test RDKitToolkitWrapper generate_conformers()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        smiles = "[H]C([H])([H])C([H])([H])[H]"
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers()
        # TODO: Make this test more robust

    def test_generate_multiple_conformers(self):
        """Test RDKitToolkitWrapper generate_conformers() for generating multiple conformers"""
        toolkit_wrapper = RDKitToolkitWrapper()
        smiles = "CCCCCCCCCN"
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(
            rms_cutoff=1 * unit.angstrom,
            n_conformers=100,
            toolkit_registry=toolkit_wrapper,
        )
        assert molecule.n_conformers > 1
        assert not (molecule.conformers[0] == (0.0 * unit.angstrom)).all()

        # Ensure rms_cutoff kwarg is working
        molecule2 = toolkit_wrapper.from_smiles(smiles)
        molecule2.generate_conformers(
            rms_cutoff=0.1 * unit.angstrom,
            n_conformers=100,
            toolkit_registry=toolkit_wrapper,
        )
        assert molecule2.n_conformers > molecule.n_conformers

        # Ensure n_conformers kwarg is working
        molecule2 = toolkit_wrapper.from_smiles(smiles)
        molecule2.generate_conformers(
            rms_cutoff=0.1 * unit.angstrom,
            n_conformers=10,
            toolkit_registry=toolkit_wrapper,
        )
        assert molecule2.n_conformers == 10

    def test_generate_conformers_failure(self):
        toolkit = RDKitToolkitWrapper()

        molecule = Molecule.from_smiles("F[U](F)(F)(F)(F)F")

        with pytest.raises(ConformerGenerationError, match="RDKit conf.*fail"):
            toolkit.generate_conformers(molecule, n_conformers=1)

    def test_generate_conformers_large_molecule(self):
        """Ensure that we don't get error caused by this molecule being too big for conf gen.  See issue #882 / OpenMM #3550."""
        ql8 = Molecule.from_file(get_data_file_path("molecules/QL8.sdf"))

        ql8.generate_conformers(
            n_conformers=1,
            toolkit_registry=RDKitToolkitWrapper(),
        )

    @pytest.mark.parametrize("partial_charge_method", ["mmff94", "gasteiger"])
    def test_assign_partial_charges_neutral(self, partial_charge_method):
        """Test RDKitToolkitWrapper assign_partial_charges()"""

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])
        # TODO: create_ethanol should be replaced by a function scope fixture.
        molecule = create_ethanol()
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
            # the partial charges returned by this method already have no rounding errors that require normalization
            # so this is just testing that the kwarg remains accessible
            normalize_partial_charges=False,
        )
        charge_sum = np.sum(molecule.partial_charges)
        assert 1.0e-10 > abs(charge_sum.m_as(unit.elementary_charge))

    @pytest.mark.parametrize("partial_charge_method", ["mmff94", "gasteiger"])
    def test_assign_partial_charges_net_charge(self, partial_charge_method):
        """
        Test RDKitToolkitWrapper assign_partial_charges() on a molecule with net charge.
        """

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])
        # TODO: create_acetate should be replaced by a function scope fixture.
        molecule = create_acetate()
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
        )
        charge_sum = np.sum(molecule.partial_charges)
        assert 1.0e-10 > abs(charge_sum.m_as(unit.elementary_charge) + 1.0)

    def test_assign_partial_charges_bad_charge_method(self):
        """Test RDKitToolkitWrapper assign_partial_charges() for a nonexistent charge method"""

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])
        molecule = create_ethanol()

        # Molecule.assign_partial_charges calls the ToolkitRegistry with raise_exception_types = [],
        # which means it will only ever return ValueError
        with pytest.raises(
            ValueError, match="is not available from RDKitToolkitWrapper"
        ):
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method="NotARealChargeMethod",
            )

        # ToolkitWrappers raise a specific exception class, so we test that here
        with pytest.raises(
            ChargeMethodUnavailableError,
            match="is not available from RDKitToolkitWrapper",
        ):
            RDTKW = RDKitToolkitWrapper()
            RDTKW.assign_partial_charges(
                molecule=molecule, partial_charge_method="NotARealChargeMethod"
            )

    def test_elf_is_problematic_conformer_acid(
        self, formic_acid_molecule, formic_acid_conformers
    ):
        problematic, reason = RDKitToolkitWrapper._elf_is_problematic_conformer(
            formic_acid_molecule, formic_acid_conformers["cis"]
        )
        assert not problematic
        assert reason is None

        problematic, reason = RDKitToolkitWrapper._elf_is_problematic_conformer(
            formic_acid_molecule, formic_acid_conformers["trans"]
        )
        assert problematic
        assert reason is not None

    def test_elf_prune_problematic_conformers_acid(
        self, formic_acid_molecule, formic_acid_conformers
    ):
        formic_acid_molecule._conformers = [*formic_acid_conformers.values()]

        pruned_conformers = RDKitToolkitWrapper._elf_prune_problematic_conformers(
            formic_acid_molecule
        )

        assert len(pruned_conformers) == 1
        assert np.allclose(
            formic_acid_conformers["cis"].m_as(unit.angstrom),
            pruned_conformers[0].m_as(unit.angstrom),
        )

    def test_elf_compute_electrostatic_energy(self, formic_acid_molecule: Molecule):
        """Test the computation of the ELF electrostatic energy function."""

        # Set some partial charges and a dummy conformer with values which make
        # computing the expected energy by hand easier.
        formic_acid_molecule.partial_charges = (
            np.ones(formic_acid_molecule.n_atoms) * 1.0 * unit.elementary_charge
        )

        formic_acid_molecule.partial_charges[0] *= 2.0
        formic_acid_molecule.partial_charges[4] *= 3.0

        conformer = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
            ]
        )

        # Compute the conformers electrostatic energy.
        computed_energy = RDKitToolkitWrapper._elf_compute_electrostatic_energy(
            formic_acid_molecule, conformer * unit.angstrom
        )
        # q_O1 * q_H2 / d_O1,H2 + q_H1 * q_H2 / d_H1,H2
        expected_energy = 2.0 * 3.0 / np.sqrt(2.0) + 1.0 * 3.0 / 2.0

        assert np.isclose(computed_energy, expected_energy)

    def test_elf_compute_rms_matrix(self, formic_acid_molecule: Molecule):
        """Test the computation of the ELF conformer RMS matrix."""
        formic_acid_molecule.add_conformer(np.random.random((5, 3)) * unit.angstrom)
        formic_acid_molecule.add_conformer(np.random.random((5, 3)) * unit.angstrom)

        rms_matrix = RDKitToolkitWrapper._elf_compute_rms_matrix(formic_acid_molecule)

        assert rms_matrix.shape == (2, 2)

        assert np.isclose(rms_matrix[0, 0], 0.0)
        assert np.isclose(rms_matrix[1, 1], 0.0)

        assert np.isclose(rms_matrix[0, 1], rms_matrix[1, 0])
        assert not np.isclose(rms_matrix[0, 1], 0.0)

    def test_elf_compute_rms_matrix_symmetry(self):
        """Test the computation of the ELF conformer RMS matrix for matrices which
        contain symmetry."""

        # Create a molecule which can have two different automorphs.
        n_methyl_aniline: Molecule = Molecule.from_smiles("CNc1ccccc1")
        n_methyl_aniline.generate_conformers(n_conformers=1)

        # Add a second conformer with the benzene ring flipped 180
        original_conformer = n_methyl_aniline.conformers[0].m_as(unit.angstrom)

        ring_atoms = RDKitToolkitWrapper().find_smarts_matches(
            n_methyl_aniline,
            "[#6]-[#7](-[#6]1:[#6:1](-[#1:2]):[#6:3](-[#1:4]):[#6]:[#6:6](-[#1:5]):[#6:8](-[#1:7])1)",
        )[0]

        flipped_conformer = np.copy(original_conformer)

        for i in range(8):
            flipped_conformer[ring_atoms[i], :] = original_conformer[ring_atoms[7 - i]]

        n_methyl_aniline.add_conformer(flipped_conformer * unit.angstrom)

        # Compute the RMS matrix.
        rms_matrix = RDKitToolkitWrapper._elf_compute_rms_matrix(n_methyl_aniline)

        assert rms_matrix.shape == (2, 2)
        assert np.allclose(rms_matrix, 0.0, atol=1e-7)

    @pytest.mark.parametrize(
        "expected_conformer_map, rms_tolerance",
        [({0: 0, 1: 2}, 0.001 * unit.angstrom), ({0: 0}, 100.0 * unit.angstrom)],
    )
    def test_elf_select_diverse_conformers(
        self,
        formic_acid_molecule: Molecule,
        expected_conformer_map: Dict[int, int],
        rms_tolerance: unit.Quantity,
    ):
        """Test the greedy selection of 'diverse' ELF conformers."""

        formic_acid_molecule.add_conformer(np.random.random((5, 3)) * unit.angstrom)
        formic_acid_molecule.add_conformer(formic_acid_molecule.conformers[0] * 1.1)
        formic_acid_molecule.add_conformer(formic_acid_molecule.conformers[0] * 1.2)

        conformers = RDKitToolkitWrapper._elf_select_diverse_conformers(
            formic_acid_molecule, formic_acid_molecule.conformers, 2, rms_tolerance
        )

        assert len(conformers) == len(expected_conformer_map)

        for elf_index, original_index in expected_conformer_map.items():
            assert np.allclose(
                conformers[elf_index].m_as(unit.angstrom),
                formic_acid_molecule.conformers[original_index].m_as(unit.angstrom),
            )

    def test_apply_elf_conformer_selection(self):
        """Test applying the ELF10 method."""

        toolkit = RDKitToolkitWrapper()

        molecule = Molecule.from_file(
            get_data_file_path(os.path.join("molecules", "z_3_hydroxy_propenal.sdf")),
            "SDF",
        )

        # Test that the simple case of no conformers does not yield an exception.
        toolkit.apply_elf_conformer_selection(molecule)

        initial_conformers = [
            # Add a conformer with an internal H-bond.
            unit.Quantity(
                np.array(
                    [
                        [0.5477, 0.3297, -0.0621],
                        [-0.1168, -0.7881, 0.2329],
                        [-1.4803, -0.8771, 0.1667],
                        [-0.2158, 1.5206, -0.4772],
                        [-1.4382, 1.5111, -0.5580],
                        [1.6274, 0.3962, -0.0089],
                        [0.3388, -1.7170, 0.5467],
                        [-1.8612, -0.0347, -0.1160],
                        [0.3747, 2.4222, -0.7115],
                    ]
                ),
                unit.angstrom,
            ),
            # Add a conformer without an internal H-bond.
            unit.Quantity(
                np.array(
                    [
                        [0.5477, 0.3297, -0.0621],
                        [-0.1168, -0.7881, 0.2329],
                        [-1.4803, -0.8771, 0.1667],
                        [-0.2158, 1.5206, -0.4772],
                        [0.3353, 2.5772, -0.7614],
                        [1.6274, 0.3962, -0.0089],
                        [0.3388, -1.7170, 0.5467],
                        [-1.7743, -1.7634, 0.4166],
                        [-1.3122, 1.4082, -0.5180],
                    ]
                ),
                unit.angstrom,
            ),
        ]

        molecule._conformers = [*initial_conformers]

        # Apply ELF10
        toolkit.apply_elf_conformer_selection(molecule)
        elf10_conformers = molecule.conformers

        assert len(elf10_conformers) == 1

        assert np.allclose(
            elf10_conformers[0].m_as(unit.angstrom),
            initial_conformers[1].m_as(unit.angstrom),
        )

    def test_apply_elf_conformer_selection_acid(
        self, formic_acid_molecule, formic_acid_conformers, caplog
    ):
        """Test applying the ELF10 method."""

        toolkit = RDKitToolkitWrapper()

        # Add the conformers to the molecule and apply ELF.
        formic_acid_molecule._conformers = [
            formic_acid_conformers["trans"],
            formic_acid_conformers["cis"],
        ]

        # Only the CIS conformer should remain after pruning and a warning raised to
        # explain why the conformer was discarded.
        with caplog.at_level(logging.WARNING):
            toolkit.apply_elf_conformer_selection(formic_acid_molecule)

        assert formic_acid_molecule.n_conformers == 1
        assert "Discarding conformer 0" in caplog.text
        assert "Molecules which contain COOH functional groups in a" in caplog.text

        assert np.allclose(
            formic_acid_molecule.conformers[0].m_as(unit.angstrom),
            formic_acid_conformers["cis"].m_as(unit.angstrom),
        )

        # Check that an exception is raised if no conformers remain after removing the
        # trans conformer.
        formic_acid_molecule._conformers = [formic_acid_conformers["trans"]]

        with pytest.raises(ValueError) as error_info:
            toolkit.apply_elf_conformer_selection(formic_acid_molecule)

        assert (
            "There were no conformers to select from after discarding conformers"
            in str(error_info.value)
        )

    def test_find_rotatable_bonds(self):
        """Test finding rotatable bonds while ignoring some groups"""

        # test a simple molecule
        ethanol = create_ethanol()
        bonds = ethanol.find_rotatable_bonds()
        assert len(bonds) == 2
        for bond in bonds:
            assert ethanol.atoms[bond.atom1_index].atomic_number != 1
            assert ethanol.atoms[bond.atom2_index].atomic_number != 1

        # now ignore the C-O bond, forwards
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups="[#6:1]-[#8:2]")
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 6

        # now ignore the O-C bond, backwards
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups="[#8:1]-[#6:2]")
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 6

        # now ignore the C-C bond
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups="[#6:1]-[#6:2]")
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 8

        # ignore a list of searches, forward
        bonds = ethanol.find_rotatable_bonds(
            ignore_functional_groups=["[#6:1]-[#8:2]", "[#6:1]-[#6:2]"]
        )
        assert bonds == []

        # ignore a list of searches, backwards
        bonds = ethanol.find_rotatable_bonds(
            ignore_functional_groups=["[#6:1]-[#6:2]", "[#8:1]-[#6:2]"]
        )
        assert bonds == []

        # test  molecules that should have no rotatable bonds
        cyclohexane = create_cyclohexane()
        bonds = cyclohexane.find_rotatable_bonds()
        assert bonds == []

        methane = Molecule.from_smiles("C")
        bonds = methane.find_rotatable_bonds()
        assert bonds == []

        ethene = Molecule.from_smiles("C=C")
        bonds = ethene.find_rotatable_bonds()
        assert bonds == []

        terminal_forwards = "[*]~[*:1]-[X2H1,X3H2,X4H3:2]-[#1]"
        terminal_backwards = "[#1]-[X2H1,X3H2,X4H3:1]-[*:2]~[*]"
        # test removing terminal rotors
        toluene = Molecule.from_file(get_data_file_path("molecules/toluene.sdf"))
        bonds = toluene.find_rotatable_bonds()
        assert len(bonds) == 1
        assert toluene.atoms[bonds[0].atom1_index].atomic_number == 6
        assert toluene.atoms[bonds[0].atom2_index].atomic_number == 6

        # find terminal bonds forward
        bonds = toluene.find_rotatable_bonds(ignore_functional_groups=terminal_forwards)
        assert bonds == []

        # find terminal bonds backwards
        bonds = toluene.find_rotatable_bonds(
            ignore_functional_groups=terminal_backwards
        )
        assert bonds == []

    @pytest.mark.parametrize(
        ("smiles", "n_atom_rings", "n_bond_rings"),
        [
            ("c1ccc2ccccc2c1", 10, 11),
            ("c1ccc(cc1)c2ccccc2", 12, 12),
            ("Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C", 30, 30),
        ],
    )
    def test_is_in_ring(self, smiles, n_atom_rings, n_bond_rings):
        """Test Atom.is_in_ring and Bond.is_in_ring"""
        mol = Molecule.from_smiles(smiles)

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper()])

        atoms_in_ring = [
            atom
            for atom in mol.atoms
            if atom.is_in_ring(toolkit_registry=toolkit_registry)
        ]

        bonds_in_ring = [
            bond
            for bond in mol.bonds
            if bond.is_in_ring(toolkit_registry=toolkit_registry)
        ]

        assert len(atoms_in_ring) == n_atom_rings
        assert len(bonds_in_ring) == n_bond_rings

    def test_central_biphenyl_bond(self):
        """Test that `Bond.is_in_ring` is False for the central bond in a phenyl"""
        # Use a mapped smiles to ensure atom order while looking up central bond
        # Generated via Molecule.from_smiles("c1ccc(cc1)c2ccccc2").to_smiles(mapped=True)

        biphenyl = Molecule.from_mapped_smiles(
            "[H:13][c:1]1[c:2]([c:3]([c:4]([c:5]([c:6]1[H:17])[H:16])[c:7]2[c:8]([c:9]([c:10]"
            "([c:11]([c:12]2[H:22])[H:21])[H:20])[H:19])[H:18])[H:15])[H:14]"
        )

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper()])

        this_bond = biphenyl.get_bond_between(3, 6)
        assert this_bond.is_in_ring(toolkit_registry=toolkit_registry) is False

        ring_bonds = [
            b for b in biphenyl.bonds if b.is_in_ring(toolkit_registry=toolkit_registry)
        ]
        assert len(ring_bonds) == 12

    def test_unattached_is_in_ring(self):
        toolkit = RDKitToolkitWrapper()
        dummy_atom = Atom(1, 0, False)

        with pytest.raises(NotAttachedToMoleculeError, match="Atom"):
            toolkit.atom_is_in_ring(dummy_atom)

        # The Bond constructor checks to see that the atoms are in a molecule,
        # so not so straightforward to make one from the constructor
        dummy_bond = Molecule.from_smiles("O").bonds[0]
        dummy_bond._molecule = None

        with pytest.raises(NotAttachedToMoleculeError, match="Bond"):
            toolkit.bond_is_in_ring(dummy_bond)

    def test_find_matches_unique(self):
        """Test the expected behavior of the `unique` argument in find_matches"""
        smirks = "[C:1]~[C:2]~[C:3]"
        tk = RDKitToolkitWrapper()

        mol = Molecule.from_smiles("CCC")

        assert len(tk.find_smarts_matches(mol, smirks, unique=True)) == 1
        assert len(tk.find_smarts_matches(mol, smirks, unique=False)) == 2

    def test_to_rdkit_losing_aromaticity_(self):
        # test the example given in issue #513
        # <https://github.com/openforcefield/openff-toolkit/issues/513>
        smiles = (
            "[H]c1c(c(c(c(c1OC2=C(C(=C(N3C2=C(C(=C3[H])C#N)[H])[H])F)[H])OC([H])([H])C([H])([H])"
            "N4C(=C(C(=O)N(C4=O)[H])[H])[H])[H])F)[H]"
        )

        mol = Molecule.from_smiles(smiles)
        rdmol = mol.to_rdkit()

        # now make sure the aromaticity matches for each atom
        for offatom, rdatom in zip(mol.atoms, rdmol.GetAtoms()):
            assert offatom.is_aromatic is rdatom.GetIsAromatic()

    @pytest.mark.slow
    def test_max_substructure_matches_can_handle_large_molecule(self):
        """Test RDKitToolkitWrapper substructure search handles more than the default of maxMatches = 1000
        See https://github.com/openforcefield/openff-toolkit/pull/509 .
        """
        tk = RDKitToolkitWrapper()
        smiles = "C" * 3000
        molecule = tk.from_smiles(smiles)
        query = "[C:1]~[C:2]"
        ret = molecule.chemical_environment_matches(query, toolkit_registry=tk)
        assert len(ret) == 5998
        assert len(ret[0]) == 2

        # TODO: Add test for higher bonds orders
        # TODO: Add test for aromaticity
        # TODO: Add test and molecule functionality for isotopes
        # TODO: Add read tests for MOL/SDF, SMI
        # TODO: Add read tests fpr multi-SMI files
        # TODO: Add read tests for both files and file-like objects
        # TODO: Add read/write tests for gzipped files
        # TODO: Add write tests for all formats


@requires_ambertools
@requires_rdkit
class TestAmberToolsToolkitWrapper:
    """Test the AmberToolsToolkitWrapper"""

    def test_assign_partial_charges_am1bcc(self):
        """Test AmberToolsToolkitWrapper assign_partial_charges() with am1bcc"""
        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )
        molecule = create_ethanol()
        molecule.assign_partial_charges(
            partial_charge_method="am1bcc", toolkit_registry=toolkit_registry
        )
        charge_sum = np.sum(molecule.partial_charges)
        abs_charge_sum = np.sum(abs(molecule.partial_charges))

        assert abs(charge_sum) < 1e-10 * unit.elementary_charge
        assert abs_charge_sum > 0.25 * unit.elementary_charge

    def test_assign_partial_charges_am1bcc_no_normalization(self):
        """Test AmberToolsToolkitWrapper assign_partial_charges() with am1bcc, with
        normalize_partial_charges=False"""
        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )
        # Use a cyclic N3H3 molecule since this is likely to produce a rounding error
        # Antechamber outputs 6 digits after the decimal in charges.txt, so I (Jeff) don't know
        # why this N3H3 molecule ends up with an error of 1e-3, but it's the smallest reproducing
        # case of this that I could find.
        molecule = create_cyclic_n3h3()
        molecule.assign_partial_charges(
            partial_charge_method="am1bcc",
            toolkit_registry=toolkit_registry,
            normalize_partial_charges=False,
        )

        charge_sum = np.sum(molecule.partial_charges)
        abs_charge_sum = np.sum(abs(molecule.partial_charges))

        # Rounding error should be on the order of 1e-3
        assert 1e-2 > abs(charge_sum.m_as(unit.elementary_charge)) > 1e-4
        assert abs_charge_sum > 0.25 * unit.elementary_charge

    def test_assign_partial_charges_am1bcc_net_charge(self):
        """Test AmberToolsToolkitWrapper assign_partial_charges() on a molecule with a net -1 charge"""
        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )
        molecule = create_acetate()
        molecule.assign_partial_charges(
            partial_charge_method="am1bcc", toolkit_registry=toolkit_registry
        )
        charge_sum = np.sum(molecule.partial_charges)
        assert 1e-10 > abs(charge_sum.m_as(unit.elementary_charge) + 1)

    def test_assign_partial_charges_am1bcc_wrong_n_confs(self):
        """
        Test AmberToolsToolkitWrapper assign_partial_charges() with am1bcc when requesting to use an incorrect number
        of conformers
        """

        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=2, rms_cutoff=0.01 * unit.angstrom)

        # Try passing in the incorrect number of confs, but without specifying strict_n_conformers,
        # which should produce a warning
        with pytest.warns(
            IncorrectNumConformersWarning,
            match="has 2 conformers, but charge method 'am1bcc' expects exactly 1.",
        ):
            molecule.assign_partial_charges(
                partial_charge_method="am1bcc",
                toolkit_registry=toolkit_registry,
                use_conformers=molecule.conformers,
                strict_n_conformers=False,
            )

        # Try again, with strict_n_confs as true, but not including use_confs, so the
        # recommended number of confs will be generated
        molecule.assign_partial_charges(
            partial_charge_method="am1bcc",
            toolkit_registry=toolkit_registry,
            strict_n_conformers=True,
        )

        # Test calling the ToolkitWrapper _indirectly_, though the Molecule API,
        # which should raise the first error encountered
        with pytest.raises(
            ValueError,
            match="has 2 conformers, but charge method 'am1bcc' expects exactly 1.",
        ):
            molecule.assign_partial_charges(
                partial_charge_method="am1bcc",
                toolkit_registry=toolkit_registry,
                use_conformers=molecule.conformers,
                strict_n_conformers=True,
            )

        # Test calling the ToolkitWrapper _indirectly_, though a ToolkitRegistry,
        # specifying raise_exception_types=[]
        # which should aggregate any exceptions and bundle all of the messages
        # in a failed task together in a single ValueError.
        with pytest.raises(
            ValueError,
            match="has 2 conformers, but charge method 'am1bcc' expects exactly 1.",
        ):
            toolkit_registry.call(
                "assign_partial_charges",
                partial_charge_method="am1bcc",
                molecule=molecule,
                use_conformers=molecule.conformers,
                strict_n_conformers=True,
                raise_exception_types=[],
            )

        # Test calling the ToolkitWrapper _directly_, passing in the incorrect number of
        # confs, and specify strict_n_conformers, which should produce an IncorrectNumConformersError
        with pytest.raises(
            IncorrectNumConformersError,
            match="has 2 conformers, but charge method 'am1bcc' expects exactly 1.",
        ):
            ATTKW = AmberToolsToolkitWrapper()
            ATTKW.assign_partial_charges(
                partial_charge_method="am1bcc",
                molecule=molecule,
                use_conformers=molecule.conformers,
                strict_n_conformers=True,
            )

    @pytest.mark.parametrize(
        "partial_charge_method", ["am1bcc", "am1-mulliken", "gasteiger"]
    )
    def test_assign_partial_charges_neutral(self, partial_charge_method):
        """Test AmberToolsToolkitWrapper assign_partial_charges()"""

        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )
        molecule = create_ethanol()
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
        )
        charge_sum = np.sum(molecule.partial_charges)
        assert 1e-10 > charge_sum.m_as(unit.elementary_charge)

    @pytest.mark.xfail(strict=False)
    @pytest.mark.parametrize("partial_charge_method", ["am1bcc", "am1-mulliken"])
    def test_assign_partial_charges_conformer_dependence(self, partial_charge_method):
        """Test AmberToolsToolkitWrapper assign_partial_charges()'s use_conformers kwarg
        to ensure charges are really conformer dependent. Skip Gasteiger because it isn't
        conformer dependent."""
        import copy

        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=1)
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
            use_conformers=molecule.conformers,
        )
        pcs1 = copy.deepcopy(molecule.partial_charges)
        # This test case needs a pretty extreme coordinate change since ambertools only
        # stores partial charges to 1e-3
        molecule._conformers[0][0][0] += 3.0 * unit.angstrom
        molecule._conformers[0][1][1] += 3.0 * unit.angstrom
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
            use_conformers=molecule.conformers,
        )
        for pc1, pc2 in zip(pcs1, molecule.partial_charges):
            assert abs(pc1 - pc2) > 1.0e-3 * unit.elementary_charge

    @pytest.mark.parametrize(
        "partial_charge_method", ["am1bcc", "am1-mulliken", "gasteiger"]
    )
    def test_assign_partial_charges_net_charge(self, partial_charge_method):
        """
        Test AmberToolsToolkitWrapper assign_partial_charges().
        """

        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )
        molecule = create_acetate()
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
        )
        charge_sum = np.sum(molecule.partial_charges)
        assert 1e-10 > abs(charge_sum.m_as(unit.elementary_charge) + 1)

    def test_assign_partial_charges_bad_charge_method(self):
        """Test AmberToolsToolkitWrapper assign_partial_charges() for a nonexistent charge method"""

        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )
        molecule = create_ethanol()

        # For now, ToolkitRegistries lose track of what exception type
        # was thrown inside them, so we just check for a ValueError here
        with pytest.raises(
            ValueError, match="is not available from AmberToolsToolkitWrapper"
        ):
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method="NotARealChargeMethod",
            )

        # ToolkitWrappers raise a specific exception class, so we test that here
        with pytest.raises(
            ChargeMethodUnavailableError,
            match="is not available from AmberToolsToolkitWrapper",
        ):
            ATTKW = AmberToolsToolkitWrapper()
            ATTKW.assign_partial_charges(
                molecule=molecule, partial_charge_method="NotARealChargeMethod"
            )

    @pytest.mark.parametrize(
        "partial_charge_method, expected_n_confs, toolkit_wrappers",
        [
            ("am1bcc", 1, [AmberToolsToolkitWrapper, RDKitToolkitWrapper]),
            ("am1-mulliken", 1, [AmberToolsToolkitWrapper, RDKitToolkitWrapper]),
            ("gasteiger", 0, [AmberToolsToolkitWrapper]),
        ],
    )
    def test_assign_partial_charges_wrong_n_confs(
        self, partial_charge_method, expected_n_confs, toolkit_wrappers
    ):
        """
        Test AmberToolsToolkitWrapper assign_partial_charges() when requesting to use an incorrect number of
        conformers
        """

        toolkit_registry = ToolkitRegistry(toolkit_precedence=toolkit_wrappers)
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=2, rms_cutoff=0.01 * unit.angstrom)

        # Try passing in the incorrect number of confs, but without specifying strict_n_conformers,
        # which should produce a warning
        with pytest.warns(
            IncorrectNumConformersWarning,
            match=f"has 2 conformers, but charge method '{partial_charge_method}' "
            f"expects exactly {expected_n_confs}.",
        ):
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method=partial_charge_method,
                use_conformers=molecule.conformers,
                strict_n_conformers=False,
            )

        # Try again, with strict_n_confs as true, but not including use_confs, so the
        # recommended number of confs will be generated
        molecule.assign_partial_charges(
            toolkit_registry=toolkit_registry,
            partial_charge_method=partial_charge_method,
            strict_n_conformers=True,
        )

        # Test calling the ToolkitWrapper _indirectly_, though the Molecule API
        # which should aggregate any exceptions and bundle all of the messages
        # in a failed task together in a single ValueError.
        with pytest.raises(
            ValueError,
            match=f"has 2 conformers, but charge method '{partial_charge_method}' "
            f"expects exactly {expected_n_confs}.",
        ):
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method=partial_charge_method,
                use_conformers=molecule.conformers,
                strict_n_conformers=True,
            )

        # Test calling the ToolkitWrapper _directly_, passing in the incorrect number of
        # confs, and specify strict_n_conformers, which should produce an IncorrectNumConformersError
        with pytest.raises(
            IncorrectNumConformersError,
            match=f"has 2 conformers, but charge method '{partial_charge_method}' "
            f"expects exactly {expected_n_confs}.",
        ):
            ATTKW = AmberToolsToolkitWrapper()
            ATTKW.assign_partial_charges(
                molecule=molecule,
                partial_charge_method=partial_charge_method,
                use_conformers=molecule.conformers,
                strict_n_conformers=True,
            )

    @pytest.mark.parametrize("bond_order_model", ["am1-wiberg"])
    @pytest.mark.parametrize(
        "smiles",
        [
            "[H]C([H])([H])C([H])([H])[H]",
            "[H]C([H])([H])[N+]([H])([H])[H]",
            r"C\C(F)=C(/F)C[C@@](C)(Cl)Br",
        ],
    )
    def test_assign_fractional_bond_orders(self, bond_order_model, smiles):
        """Test AmbetToolsToolkitWrapper assign_fractional_bond_orders()"""

        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )

        molecule = toolkit_registry.call("from_smiles", smiles)
        molecule.assign_fractional_bond_orders(
            toolkit_registry=toolkit_registry, bond_order_model=bond_order_model
        )
        # TODO: Add test for equivalent Wiberg orders for equivalent bonds

        # Sanity check single bonds.
        assert all(
            0.75 < bond.fractional_bond_order < 1.25
            for bond in molecule.bonds
            if bond.bond_order == 1
        )
        # Sanity check double bonds.
        assert all(
            1.75 < bond.fractional_bond_order < 2.25
            for bond in molecule.bonds
            if bond.bond_order == 2
        )

    def test_assign_fractional_bond_orders_conformer_dependence(self):
        """
        Test that RDKitToolkitWrapper assign_fractional_bond_orders() provides different results when using
        different conformers
        """

        toolkit_wrapper = ToolkitRegistry(
            [RDKitToolkitWrapper, AmberToolsToolkitWrapper]
        )
        # Get the WBOs using one conformer
        molecule = create_ethanol()
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        molecule.assign_fractional_bond_orders(
            toolkit_registry=toolkit_wrapper,
            use_conformers=molecule.conformers,
            bond_order_model="am1-wiberg",
        )

        # Do the same again, but change the conformer to yield a different result
        molecule_diff_coords = create_ethanol()
        molecule_diff_coords.generate_conformers(toolkit_registry=toolkit_wrapper)
        molecule_diff_coords._conformers[0][0][0] = (
            molecule_diff_coords._conformers[0][0][0] + 1.0 * unit.angstrom
        )
        molecule_diff_coords._conformers[0][1][0] = (
            molecule_diff_coords._conformers[0][1][0] - 1.0 * unit.angstrom
        )
        molecule_diff_coords._conformers[0][2][0] = (
            molecule_diff_coords._conformers[0][2][0] + 1.0 * unit.angstrom
        )
        molecule_diff_coords.assign_fractional_bond_orders(
            toolkit_registry=toolkit_wrapper,
            use_conformers=molecule_diff_coords.conformers,
            bond_order_model="am1-wiberg",
        )

        for bond1, bond2 in zip(molecule.bonds, molecule_diff_coords.bonds):
            assert abs(bond1.fractional_bond_order - bond2.fractional_bond_order) > 1e-3

    @pytest.mark.parametrize("bond_order_model", ["am1-wiberg"])
    def test_assign_fractional_bond_orders_neutral_charge_mol(self, bond_order_model):
        """Test AmberToolsToolkitWrapper assign_fractional_bond_orders() for neutral and charged molecule.
        Also tests using existing conformers"""

        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )
        # Reading neutral molecule from file
        filename = get_data_file_path("molecules/CID20742535_neutral.sdf")
        molecule1 = Molecule.from_file(filename)
        # Reading negative molecule from file
        filename = get_data_file_path("molecules/CID20742535_anion.sdf")
        molecule2 = Molecule.from_file(filename)

        # Checking that only one additional bond is present in the neutral molecule
        assert len(molecule1.bonds) == len(molecule2.bonds) + 1

        molecule1.assign_fractional_bond_orders(
            toolkit_registry=toolkit_registry,
            bond_order_model=bond_order_model,
            use_conformers=molecule1.conformers,
        )

        for i in molecule1.bonds:
            if i.is_aromatic:
                # Checking aromatic bonds
                assert 1.05 < i.fractional_bond_order < 1.65
            elif i.atom1.atomic_number == 1 or i.atom2.atomic_number == 1:
                # Checking bond order of C-H or O-H bonds are around 1
                assert 0.85 < i.fractional_bond_order < 1.05
            elif i.atom1.atomic_number == 8 or i.atom2.atomic_number == 8:
                # Checking C-O single bond
                wbo_C_O_neutral = i.fractional_bond_order
                assert 1.0 < wbo_C_O_neutral < 1.5
            else:
                # Should be C-C single bond
                assert (i.atom1_index == 4 and i.atom2_index == 6) or (
                    i.atom1_index == 6 and i.atom2_index == 4
                )
                wbo_C_C_neutral = i.fractional_bond_order
                assert 1.0 < wbo_C_C_neutral < 1.3

        molecule2.assign_fractional_bond_orders(
            toolkit_registry=toolkit_registry,
            bond_order_model=bond_order_model,
            use_conformers=molecule2.conformers,
        )
        for i in molecule2.bonds:
            if i.is_aromatic:
                # Checking aromatic bonds
                assert 1.05 < i.fractional_bond_order < 1.65

            elif i.atom1.atomic_number == 1 or i.atom2.atomic_number == 1:
                # Checking bond order of C-H or O-H bonds are around 1
                assert 0.85 < i.fractional_bond_order < 1.05
            elif i.atom1.atomic_number == 8 or i.atom2.atomic_number == 8:
                # Checking C-O single bond
                wbo_C_O_anion = i.fractional_bond_order
                assert 1.3 < wbo_C_O_anion < 1.8
            else:
                # Should be C-C single bond
                assert (i.atom1_index == 4 and i.atom2_index == 6) or (
                    i.atom1_index == 6 and i.atom2_index == 4
                )
                wbo_C_C_anion = i.fractional_bond_order
                assert 1.0 < wbo_C_C_anion < 1.3

        # Wiberg bond order of C-C single bond is higher in the anion
        assert wbo_C_C_anion > wbo_C_C_neutral
        # Wiberg bond order of C-O bond is higher in the anion
        assert wbo_C_O_anion > wbo_C_O_neutral

    def test_assign_fractional_bond_orders_invalid_method(self):
        """
        Test that AmberToolsToolkitWrapper.assign_fractional_bond_orders() raises the
        correct error if an invalid charge model is provided
        """

        toolkit_registry = ToolkitRegistry(
            toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        )
        molecule = toolkit_registry.call("from_smiles", "C")

        expected_error = (
            "Bond order model 'not a real charge model' is not supported by "
            "AmberToolsToolkitWrapper. Supported models are ([[]'am1-wiberg'[]])"
        )
        with pytest.raises(ValueError, match=expected_error):
            molecule.assign_fractional_bond_orders(
                toolkit_registry=AmberToolsToolkitWrapper(),
                bond_order_model="not a real charge model",
            )

    @requires_openeye
    def test_assign_fractional_bond_orders_openeye_installed(self):
        """Test that assign_fractional_bond_orders produces the same result
        with and without OpenEye toolkits installed"""
        mol = Molecule.from_smiles("CCO")
        AmberToolsToolkitWrapper().assign_fractional_bond_orders(mol)
        with_oe = [b.fractional_bond_order for b in mol.bonds]
        GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)
        AmberToolsToolkitWrapper().assign_fractional_bond_orders(mol)
        without_oe = [b.fractional_bond_order for b in mol.bonds]
        GLOBAL_TOOLKIT_REGISTRY.register_toolkit(OpenEyeToolkitWrapper)

        assert with_oe == pytest.approx(without_oe, abs=1e-5)


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
        charge_sum = np.sum(molecule.partial_charges)
        assert 1.0e-10 > abs(charge_sum.m_as(unit.elementary_charge))

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
        charge_sum = np.sum(molecule.partial_charges)
        assert 1.0e-10 > abs(charge_sum.m_as(unit.elementary_charge) + 1.0)

    def test_assign_partial_charges_bad_charge_method(self):
        """Test BuiltInToolkitWrapper assign_partial_charges() for a nonexistent charge method"""

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[BuiltInToolkitWrapper])
        molecule = create_ethanol()

        # For now, the Molecule API passes raise_exception_types=[] to ToolkitRegistry.call,
        # which loses track of what exception type
        # was thrown inside them, so we just check for a ValueError here
        with pytest.raises(
            ValueError, match="is not supported by the Built-in toolkit"
        ):
            molecule.assign_partial_charges(
                toolkit_registry=toolkit_registry,
                partial_charge_method="NotARealChargeMethod",
            )

        # ToolkitWrappers raise a specific exception class, so we test that here
        with pytest.raises(
            ChargeMethodUnavailableError,
            match="is not supported by the Built-in toolkit",
        ):
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
            match="has 1 conformers, but charge method 'zeros' expects exactly 0.",
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
            match="has 1 conformers, but charge method 'zeros' expects exactly 0.",
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

        # Test molecule with no conformers
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

        # Test molecule with conformers
        # Add some conformers
        mol.generate_conformers(n_conformers=1)
        for _ in range(9):
            mol.add_conformer(mol.conformers[0])

        # Check with no min or max should pass
        tkw._check_n_conformers(mol, "nocharge")

        # min_confs checks
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

        # max_confs checks
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

        # min_confs and max_confs checks
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
        no_precedence = ToolkitRegistry(_register_imported_toolkit_wrappers=True)

        assert len(no_precedence.registered_toolkits) == 4

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
        assert (
            len(GLOBAL_TOOLKIT_REGISTRY.registered_toolkits) == num_toolkits - 1  # noqa
        )

        GLOBAL_TOOLKIT_REGISTRY = deepcopy(global_registry_copy)  # noqa

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
        resolved = registry.resolve("assign_partial_charges")
        assert resolved == registry.registered_toolkits[0].assign_partial_charges

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

    from openff.toolkit.utils.openeye_wrapper import requires_openeye_module

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
