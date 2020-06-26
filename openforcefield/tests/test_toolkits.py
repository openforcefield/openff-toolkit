#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for cheminformatics toolkit wrappers

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

from simtk import unit
import numpy as np
from numpy.testing import assert_almost_equal
from tempfile import NamedTemporaryFile

import pytest

from openforcefield.utils.toolkits import (OpenEyeToolkitWrapper, RDKitToolkitWrapper,
                                           AmberToolsToolkitWrapper, BuiltInToolkitWrapper,
                                           ToolkitRegistry, ToolkitWrapper,
                                           GAFFAtomTypeWarning, UndefinedStereochemistryError,
                                           ChargeMethodUnavailableError, IncorrectNumConformersError,
                                           IncorrectNumConformersWarning, InvalidToolkitError,
                                           ToolkitUnavailableException)

from openforcefield.utils import get_data_file_path
from openforcefield.topology.molecule import Molecule
from openforcefield.tests.test_forcefield import create_ethanol, create_cyclohexane, create_acetaldehyde, \
    create_reversed_ethanol, create_acetate

#=============================================================================================
# FIXTURES
#=============================================================================================


def get_mini_drug_bank(toolkit_class):
    """Read the mini drug bank sdf file with the toolkit and return the molecules"""

    # This is a work around a weird error where even though the test is skipped due to a missing toolkit
    #  we still try and read the file with the toolkit
    if toolkit_class.is_available():
        toolkit = toolkit_class()
        molecules = Molecule.from_file(get_data_file_path('molecules/MiniDrugBank.sdf'), 'sdf', toolkit_registry=toolkit,
                                       allow_undefined_stereo=True)
    else:
        molecules = []

    return molecules


openeye_inchi_stereochemistry_lost = ['DrugBank_2799', 'DrugBank_5414', 'DrugBank_5415', 'DrugBank_5418',
                                      'DrugBank_2955', 'DrugBank_2987', 'DrugBank_5555', 'DrugBank_472',
                                      'DrugBank_5737', 'DrugBank_3332', 'DrugBank_3461', 'DrugBank_794',
                                      'DrugBank_3502', 'DrugBank_6026', 'DrugBank_3622', 'DrugBank_977',
                                      'DrugBank_3693', 'DrugBank_3726', 'DrugBank_3739', 'DrugBank_6222',
                                      'DrugBank_6232', 'DrugBank_3844', 'DrugBank_6295', 'DrugBank_6304',
                                      'DrugBank_6305', 'DrugBank_3930', 'DrugBank_6329', 'DrugBank_6353',
                                      'DrugBank_6355', 'DrugBank_6401', 'DrugBank_4161', 'DrugBank_4162',
                                      'DrugBank_6509', 'DrugBank_6531', 'DrugBank_1570', 'DrugBank_4249',
                                      'DrugBank_1634', 'DrugBank_1659', 'DrugBank_6647', 'DrugBank_1700',
                                      'DrugBank_1721', 'DrugBank_1742', 'DrugBank_1802', 'DrugBank_6775',
                                      'DrugBank_1849', 'DrugBank_1864', 'DrugBank_6875', 'DrugBank_1897',
                                      'DrugBank_4593', 'DrugBank_1962', 'DrugBank_4662', 'DrugBank_7049',
                                      'DrugBank_4702', 'DrugBank_2095', 'DrugBank_4778', 'DrugBank_2141',
                                      'DrugBank_2148', 'DrugBank_2178', 'DrugBank_4865', 'DrugBank_2208',
                                      'DrugBank_2210', 'DrugBank_2276', 'DrugBank_4959', 'DrugBank_4964',
                                      'DrugBank_5043', 'DrugBank_2429', 'DrugBank_5076', 'DrugBank_2465',
                                      'DrugBank_2519', 'DrugBank_2538', 'DrugBank_5158', 'DrugBank_5176',
                                      'DrugBank_2592']

openeye_inchi_isomorphic_fails = ['DrugBank_1661', 'DrugBank_4346', 'DrugBank_2467']

rdkit_inchi_stereochemistry_lost = ['DrugBank_5414', 'DrugBank_2955', 'DrugBank_5737', 'DrugBank_3332', 'DrugBank_3461',
                                    'DrugBank_6026', 'DrugBank_3622', 'DrugBank_3726', 'DrugBank_6222', 'DrugBank_3844',
                                    'DrugBank_6304', 'DrugBank_6305', 'DrugBank_6329', 'DrugBank_6509', 'DrugBank_6647',
                                    'DrugBank_1897', 'DrugBank_4778', 'DrugBank_2148', 'DrugBank_2178', 'DrugBank_2538',
                                    'DrugBank_2592', 'DrugBank_4249', 'DrugBank_5076', 'DrugBank_5418', 'DrugBank_3930',
                                    'DrugBank_1634', 'DrugBank_1962', 'DrugBank_5043', 'DrugBank_2519']

rdkit_inchi_isomorphic_fails = ['DrugBank_178', 'DrugBank_246', 'DrugBank_5847', 'DrugBank_700', 'DrugBank_1564',
                                'DrugBank_1700', 'DrugBank_4662', 'DrugBank_2052', 'DrugBank_2077', 'DrugBank_2082',
                                'DrugBank_2210', 'DrugBank_2642']
#=============================================================================================
# TESTS
#=============================================================================================

class TestOpenEyeToolkitWrapper:
    """Test the OpenEyeToolkitWrapper"""

    # TODO: Make separate smiles_add_H and smiles_explicit_H tests

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_smiles(self):
        """Test OpenEyeToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()

        # This differs from RDKit's SMILES due to different canonicalization schemes

        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles,
                                        toolkit_registry=toolkit_wrapper)
        # When creating an OFFMol from SMILES, partial charges should be initialized to None
        assert molecule.partial_charges is None
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_smiles_missing_stereochemistry(self):
        """Test OpenEyeToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()

        unspec_chiral_smiles = r"C\C(F)=C(/F)CC(C)(Cl)Br"
        spec_chiral_smiles = r"C\C(F)=C(/F)C[C@@](C)(Cl)Br"
        unspec_db_smiles = r"CC(F)=C(F)C[C@@](C)(Cl)Br"
        spec_db_smiles = r"C\C(F)=C(/F)C[C@@](C)(Cl)Br"

        for title, smiles, raises_exception in [("unspec_chiral_smiles", unspec_chiral_smiles, True),
                                                ("spec_chiral_smiles", spec_chiral_smiles, False),
                                                ("unspec_db_smiles", unspec_db_smiles, True),
                                                ("spec_db_smiles", spec_db_smiles, False),
                                                ]:
            if raises_exception:
                with pytest.raises(UndefinedStereochemistryError) as context:
                    Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)
                Molecule.from_smiles(smiles,
                                     toolkit_registry=toolkit_wrapper,
                                     allow_undefined_stereo=True)
            else:
                Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)

    # TODO: test_smiles_round_trip


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_smiles_add_H(self):
        """Test OpenEyeToolkitWrapper for adding explicit hydrogens"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        # This differs from RDKit's SMILES due to different canonicalization schemes
        input_smiles = 'CC'
        expected_output_smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = Molecule.from_smiles(input_smiles,
                                        toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert expected_output_smiles == smiles2

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_smiles_charged(self):
        """Test OpenEyeToolkitWrapper functions for reading/writing charged SMILES"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        # This differs from RDKit's expected output due to different canonicalization schemes
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles,
                                        toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_to_from_openeye_core_props_filled(self):
        """Test OpenEyeToolkitWrapper to_openeye() and from_openeye()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()

        # Replacing with a simple molecule with stereochemistry
        input_smiles = r'C\C(F)=C(/F)C[C@@](C)(Cl)Br'
        expected_output_smiles = r'[H]C([H])([H])/C(=C(/C([H])([H])[C@@](C([H])([H])[H])(Cl)Br)\F)/F'
        molecule = Molecule.from_smiles(input_smiles, toolkit_registry=toolkit_wrapper)
        assert molecule.to_smiles(toolkit_registry=toolkit_wrapper) == expected_output_smiles

        # Populate core molecule property fields
        molecule.name = 'Alice'
        partial_charges = unit.Quantity(np.array([-.9, -.8, -.7, -.6,
                                                  -.5, -.4, -.3, -.2,
                                                  -.1, 0., .1, .2,
                                                  .3, .4, .5, .6,
                                                  .7, .8]), unit.elementary_charge)
        molecule.partial_charges = partial_charges
        coords = unit.Quantity(np.array([['0.0', '1.0', '2.0'], ['3.0', '4.0', '5.0'], ['6.0', '7.0', '8.0'],
                                         ['9.0', '10.0', '11.0'], ['12.0', '13.0', '14.0'],
                                         ['15.0', '16.0', '17.0'],
                                         ['18.0', '19.0', '20.0'], ['21.0', '22.0', '23.0'],
                                         ['24.0', '25.0', '26.0'],
                                         ['27.0', '28.0', '29.0'], ['30.0', '31.0', '32.0'],
                                         ['33.0', '34.0', '35.0'],
                                         ['36.0', '37.0', '38.0'], ['39.0', '40.0', '41.0'],
                                         ['42.0', '43.0', '44.0'],
                                         ['45.0', '46.0', '47.0'], ['48.0', '49.0', '50.0'],
                                         ['51.0', '52.0', '53.0']]),
                               unit.angstrom)
        molecule.add_conformer(coords)
        # Populate core atom property fields
        molecule.atoms[2].name = 'Bob'
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
            assert atom1.to_dict() == atom2.to_dict()
        for bond1, bond2 in zip(molecule.bonds, molecule2.bonds):
            assert bond1.to_dict() == bond2.to_dict()
        assert (molecule.conformers[0] == molecule2.conformers[0]).all()
        for pc1, pc2 in zip(molecule._partial_charges, molecule2._partial_charges):
            pc1_ul = pc1 / unit.elementary_charge
            pc2_ul = pc2 / unit.elementary_charge
            assert_almost_equal(pc1_ul, pc2_ul, decimal=6)
        assert molecule2.to_smiles(toolkit_registry=toolkit_wrapper) == expected_output_smiles


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_to_from_openeye_core_props_unset(self):
        """Test OpenEyeToolkitWrapper to_openeye() and from_openeye() when given empty core property fields"""
        toolkit_wrapper = OpenEyeToolkitWrapper()

        # Using a simple molecule with tetrahedral and bond stereochemistry
        input_smiles = r'C\C(F)=C(/F)C[C@](C)(Cl)Br'

        expected_output_smiles = r'[H]C([H])([H])/C(=C(/C([H])([H])[C@](C([H])([H])[H])(Cl)Br)\F)/F'
        molecule = Molecule.from_smiles(input_smiles, toolkit_registry=toolkit_wrapper)
        assert molecule.to_smiles(toolkit_registry=toolkit_wrapper) == expected_output_smiles

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
            assert atom1.to_dict() == atom2.to_dict()
        for bond1, bond2 in zip(molecule.bonds, molecule2.bonds):
            assert bond1.to_dict() == bond2.to_dict()
        # The molecule was initialized from SMILES, so mol.conformers arrays should be None for both
        assert molecule.conformers is None
        assert molecule2.conformers is None
        # The molecule was initialized from SMILES, so mol.partial_charges arrays should be None for both
        assert molecule.partial_charges is None
        assert molecule2.partial_charges is None

        assert molecule2.to_smiles(toolkit_registry=toolkit_wrapper) == expected_output_smiles

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
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
            oeatom.SetPartialCharge(float('nan'))
            break
        eth_from_oe = Molecule.from_openeye(oemol)
        assert math.isnan(eth_from_oe.partial_charges[0] / unit.elementary_charge)
        for pc in eth_from_oe.partial_charges[1:]:
            assert not math.isnan(pc / unit.elementary_charge)
        # Then, set all the OEMol's partial charges to nan, and ensure that
        # from_openeye produces an OFFMol with partial_charges = None
        for oeatom in oemol.GetAtoms():
            oeatom.SetPartialCharge(float('nan'))
        eth_from_oe = Molecule.from_openeye(oemol)
        assert eth_from_oe.partial_charges is None

        # Send the OFFMol with partial_charges = None back to OEMol, and
        # ensure that all its charges are nan
        oemol2 = eth_from_oe.to_openeye()
        for oeatom in oemol2.GetAtoms():
            assert math.isnan(oeatom.GetPartialCharge())


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
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

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_openeye_from_smiles_hydrogens_are_explicit(self):
        """
        Test to ensure that OpenEyeToolkitWrapper.from_smiles has the proper behavior with
        respect to its hydrogens_are_explicit kwarg
        """
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles_impl = "C#C"
        with pytest.raises(ValueError,
                           match="but OpenEye Toolkit interpreted SMILES 'C#C' as having implicit hydrogen") as excinfo:
            offmol = Molecule.from_smiles(smiles_impl,
                                          toolkit_registry=toolkit_wrapper,
                                          hydrogens_are_explicit=True)
        offmol = Molecule.from_smiles(smiles_impl,
                                      toolkit_registry=toolkit_wrapper,
                                      hydrogens_are_explicit=False)
        assert offmol.n_atoms == 4

        smiles_expl = "HC#CH"
        offmol = Molecule.from_smiles(smiles_expl,
                                      toolkit_registry=toolkit_wrapper,
                                      hydrogens_are_explicit=True)
        assert offmol.n_atoms == 4
        # It's debatable whether this next function should pass. Strictly speaking, the hydrogens in this SMILES
        # _are_ explicit, so allowing "hydrogens_are_explicit=False" through here is allowing a contradiction.
        # We might rethink the name of this kwarg.

        offmol = Molecule.from_smiles(smiles_expl,
                                      toolkit_registry=toolkit_wrapper,
                                      hydrogens_are_explicit=False)
        assert offmol.n_atoms == 4

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    @pytest.mark.parametrize('molecule', get_mini_drug_bank(OpenEyeToolkitWrapper))
    def test_to_inchi(self, molecule):
        """Test conversion to standard and non-standard InChI"""

        toolkit = OpenEyeToolkitWrapper()
        inchi = molecule.to_inchi(toolkit_registry=toolkit)
        non_standard = molecule.to_inchi(True, toolkit_registry=toolkit)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    @pytest.mark.parametrize('molecule', get_mini_drug_bank(OpenEyeToolkitWrapper))
    def test_to_inchikey(self, molecule):
        """Test the conversion to standard and non-standard InChIKey"""

        toolkit = OpenEyeToolkitWrapper()
        inchikey = molecule.to_inchikey(toolkit_registry=toolkit)
        non_standard_key = molecule.to_inchikey(True, toolkit_registry=toolkit)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_from_bad_inchi(self):
        """Test building a molecule from a bad InChI string"""

        toolkit = OpenEyeToolkitWrapper()
        inchi = 'InChI=1S/ksbfksfksfksbfks'
        with pytest.raises(RuntimeError):
            mol = Molecule.from_inchi(inchi, toolkit_registry=toolkit)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    @pytest.mark.parametrize('molecule', get_mini_drug_bank(OpenEyeToolkitWrapper))
    def test_non_standard_inchi_round_trip(self, molecule):
        """Test if a molecule can survive an InChi round trip test in some cases the standard InChI
        will not enough to ensure information is preserved so we test the non-standard inchi here."""

        from openforcefield.utils.toolkits import UndefinedStereochemistryError

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
                    assert molecule.is_isomorphic_with(mol2, bond_order_matching=False)
            else:
                assert molecule.is_isomorphic_with(mol2, bond_order_matching=False)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_get_sdf_coordinates(self):
        """Test OpenEyeToolkitWrapper for importing a single set of coordinates from a sdf file"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule.conformers) == 1
        assert molecule.conformers[0].shape == (15,3)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_load_multiconformer_sdf_as_separate_molecules(self):
        """
        Test OpenEyeToolkitWrapper for reading a "multiconformer" SDF, which the OFF
        Toolkit should treat as separate molecules
        """
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path('molecules/methane_multiconformer.sdf')
        molecules = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecules) == 2
        assert len(molecules[0].conformers) == 1
        assert len(molecules[1].conformers) == 1
        assert molecules[0].conformers[0].shape == (5, 3)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_load_multiconformer_sdf_as_separate_molecules_properties(self):
        """
        Test OpenEyeToolkitWrapper for reading a "multiconformer" SDF, which the OFF
        Toolkit should treat as separate molecules, and it should load their SD properties
        and partial charges separately
        """
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path('molecules/methane_multiconformer_properties.sdf')
        molecules = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecules) == 2
        assert len(molecules[0].conformers) == 1
        assert len(molecules[1].conformers) == 1
        assert molecules[0].conformers[0].shape == (5, 3)
        # The first molecule in the SDF has the following properties and charges:
        assert molecules[0].properties['test_property_key'] == 'test_property_value'
        np.testing.assert_allclose(molecules[0].partial_charges / unit.elementary_charge,
                                          [-0.108680, 0.027170, 0.027170, 0.027170, 0.027170])
        # The second molecule in the SDF has the following properties and charges:
        assert molecules[1].properties['test_property_key'] == 'test_property_value2'
        assert molecules[1].properties['another_test_property_key'] == 'another_test_property_value'
        np.testing.assert_allclose(molecules[1].partial_charges / unit.elementary_charge,
                                   [0.027170, 0.027170, 0.027170, 0.027170, -0.108680])

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_write_sdf_charges(self):
        """Test OpenEyeToolkitWrapper for writing partial charges to a sdf file"""
        from io import StringIO
        toolkit_wrapper = OpenEyeToolkitWrapper()
        ethanol = create_ethanol()
        sio = StringIO()
        ethanol.to_file(sio, 'SDF', toolkit_registry=toolkit_wrapper)
        sdf_text = sio.getvalue()
        # The output lines of interest here will look like
        # > <atom.dprop.PartialCharge>
        # -0.400000 -0.300000 -0.200000 -0.100000 0.000010 0.100000 0.200000 0.300000 0.400000
        # Parse the SDF text, grabbing the numeric line above
        sdf_split = sdf_text.split('\n')
        charge_line_found = False
        for line in sdf_split:
            if charge_line_found:
                charges = [float(i) for i in line.split()]
                break
            if '> <atom.dprop.PartialCharge>' in line:
                charge_line_found = True

        # Make sure that a charge line was ever found
        assert charge_line_found == True

        # Make sure that the charges found were correct
        assert_almost_equal(charges, [-0.4, -0.3, -0.2, -0.1, 0.00001, 0.1, 0.2, 0.3, 0.4])


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_write_sdf_no_charges(self):
        """Test OpenEyeToolkitWrapper for writing an SDF file without charges"""
        from io import StringIO
        toolkit_wrapper = OpenEyeToolkitWrapper()
        ethanol = create_ethanol()
        ethanol.partial_charges = None
        sio = StringIO()
        ethanol.to_file(sio, 'SDF', toolkit_registry=toolkit_wrapper)
        sdf_text = sio.getvalue()
        # In our current configuration, if the OFFMol doesn't have partial charges, we DO NOT want a partial charge
        # block to be written. For reference, it's possible to indicate that a partial charge is not known by writing
        # out "n/a" (or another placeholder) in the partial charge block atoms without charges.
        assert '<atom.dprop.PartialCharge>' not in sdf_text

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_sdf_properties_roundtrip(self):
        """Test OpenEyeToolkitWrapper for performing a round trip of a molecule with defined partial charges
        and entries in the properties dict to and from a sdf file"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        ethanol = create_ethanol()
        ethanol.properties['test_property'] = 'test_value'
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix='.sdf') as iofile:
            ethanol.to_file(iofile.name, file_format='SDF', toolkit_registry=toolkit_wrapper)
            ethanol2 = Molecule.from_file(iofile.name, file_format='SDF', toolkit_registry=toolkit_wrapper)
        np.testing.assert_allclose(ethanol.partial_charges / unit.elementary_charge,
                                   ethanol2.partial_charges / unit.elementary_charge)
        assert ethanol2.properties['test_property'] == 'test_value'

        # Now test with no properties or charges
        ethanol = create_ethanol()
        ethanol.partial_charges = None
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix='.sdf') as iofile:
            ethanol.to_file(iofile.name, file_format='SDF', toolkit_registry=toolkit_wrapper)
            ethanol2 = Molecule.from_file(iofile.name, file_format='SDF', toolkit_registry=toolkit_wrapper)
        assert ethanol2.partial_charges is None
        assert ethanol2.properties == {}



    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_write_multiconformer_mol_as_sdf(self):
        """
        Test OpenEyeToolkitWrapper for writing a multiconformer molecule to SDF. The OFF toolkit should only
        save the first conformer.
        """
        from io import StringIO

        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path('molecules/ethanol.sdf')
        ethanol = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        ethanol.partial_charges = np.array([-4., -3., -2., -1., 0., 1., 2., 3., 4.]) * unit.elementary_charge
        ethanol.properties['test_prop'] = 'test_value'
        new_conf = ethanol.conformers[0] + (np.ones(ethanol.conformers[0].shape) * unit.angstrom)
        ethanol.add_conformer(new_conf)
        sio = StringIO()
        ethanol.to_file(sio, 'sdf', toolkit_registry=toolkit_wrapper)
        data = sio.getvalue()
        # In SD format, each molecule ends with "$$$$"
        assert data.count('$$$$') == 1
        # A basic SDF for ethanol would be 27 lines, though the properties add three more
        assert len(data.split('\n')) == 30
        assert 'test_prop' in data
        assert '<atom.dprop.PartialCharge>' in data
        # Ensure the first conformer's first atom's X coordinate is in the file
        assert str(ethanol.conformers[0][0][0].value_in_unit(unit.angstrom))[:5] in data
        # Ensure the SECOND conformer's first atom's X coordinate is NOT in the file
        assert str(ethanol.conformers[1][0][0].in_units_of(unit.angstrom))[:5] not in data


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_get_mol2_coordinates(self):
        """Test OpenEyeToolkitWrapper for importing a single set of molecule coordinates"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.mol2')
        molecule1 = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule1.conformers) == 1
        assert molecule1.conformers[0].shape == (15, 3)
        assert_almost_equal(molecule1.conformers[0][5][1] / unit.angstrom, 22.98, decimal=2)

        # Test loading from file-like object
        with open(filename, 'r') as infile:
            molecule2 = Molecule(infile, file_format='MOL2', toolkit_registry=toolkit_wrapper)
        assert molecule1.is_isomorphic_with(molecule2)
        assert len(molecule2.conformers) == 1
        assert molecule2.conformers[0].shape == (15, 3)
        assert_almost_equal(molecule2.conformers[0][5][1] / unit.angstrom, 22.98, decimal=2)

        # Test loading from gzipped mol2
        import gzip
        with gzip.GzipFile(filename + '.gz', 'r') as infile:
            molecule3 = Molecule(infile, file_format='MOL2', toolkit_registry=toolkit_wrapper)
        assert molecule1.is_isomorphic_with(molecule3)
        assert len(molecule3.conformers) == 1
        assert molecule3.conformers[0].shape == (15, 3)
        assert_almost_equal(molecule3.conformers[0][5][1] / unit.angstrom, 22.98, decimal=2)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_get_mol2_charges(self):
        """Test OpenEyeToolkitWrapper for importing a mol2 file specifying partial charges"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path('molecules/toluene_charged.mol2')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule.conformers) == 1
        assert molecule.conformers[0].shape == (15,3)
        target_charges = unit.Quantity(np.array([-0.1342,-0.1271,-0.1271,-0.1310,
                                                 -0.1310,-0.0765,-0.0541, 0.1314,
                                                  0.1286, 0.1286, 0.1303, 0.1303,
                                                  0.0440, 0.0440, 0.0440]),
                                                unit.elementary_charge)
        for pc1, pc2 in zip(molecule._partial_charges, target_charges):
            pc1_ul = pc1 / unit.elementary_charge
            pc2_ul = pc2 / unit.elementary_charge
            assert_almost_equal(pc1_ul, pc2_ul, decimal=4)



    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_mol2_charges_roundtrip(self):
        """Test OpenEyeToolkitWrapper for performing a round trip of a molecule with partial charge to and from
        a mol2 file"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        ethanol = create_ethanol()
        # we increase the magnitude of the partial charges here, since mol2 is only
        # written to 4 digits of precision, and the default middle charge for our test ethanol is 1e-5
        ethanol.partial_charges *= 100
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix='.mol2') as iofile:
            ethanol.to_file(iofile.name, file_format='mol2', toolkit_registry=toolkit_wrapper)
            ethanol2 = Molecule.from_file(iofile.name, file_format='mol2', toolkit_registry=toolkit_wrapper)
        np.testing.assert_allclose(ethanol.partial_charges / unit.elementary_charge,
                                   ethanol2.partial_charges / unit.elementary_charge)

        # Now test with no properties or charges
        ethanol = create_ethanol()
        ethanol.partial_charges = None
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix='.mol2') as iofile:
            ethanol.to_file(iofile.name, file_format='mol2', toolkit_registry=toolkit_wrapper)
            ethanol2 = Molecule.from_file(iofile.name, file_format='mol2', toolkit_registry=toolkit_wrapper)
        assert ethanol2.partial_charges is None
        assert ethanol2.properties == {}


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_get_mol2_gaff_atom_types(self):
        """Test that a warning is raised OpenEyeToolkitWrapper when it detects GAFF atom types in a mol2 file."""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        mol2_file_path = get_data_file_path('molecules/AlkEthOH_test_filt1_ff.mol2')
        with pytest.warns(GAFFAtomTypeWarning, match='SYBYL'):
            Molecule.from_file(mol2_file_path, toolkit_registry=toolkit_wrapper)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_generate_conformers(self):
        """Test OpenEyeToolkitWrapper generate_conformers()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers()
        assert molecule.n_conformers != 0
        assert not(molecule.conformers[0] == (0.*unit.angstrom)).all()

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_generate_multiple_conformers(self):
        """Test OpenEyeToolkitWrapper generate_conformers() for generating multiple conformers"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = 'CCCCCCC'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(rms_cutoff=1*unit.angstrom,
                                     n_conformers=100,
                                     toolkit_registry=toolkit_wrapper)
        assert molecule.n_conformers > 1
        assert not(molecule.conformers[0] == (0.*unit.angstrom)).all()

        # Ensure rms_cutoff kwarg is working
        molecule2 = toolkit_wrapper.from_smiles(smiles)
        molecule2.generate_conformers(rms_cutoff=0.1*unit.angstrom,
                                      n_conformers=100,
                                      toolkit_registry=toolkit_wrapper)
        assert molecule2.n_conformers > molecule.n_conformers

        # Ensure n_conformers kwarg is working
        molecule2 = toolkit_wrapper.from_smiles(smiles)
        molecule2.generate_conformers(rms_cutoff=0.1*unit.angstrom,
                                      n_conformers=10,
                                      toolkit_registry=toolkit_wrapper)
        assert molecule2.n_conformers == 10

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_compute_partial_charges_am1bcc(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges_am1bcc()"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()
        molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry)  # , charge_model=charge_model)
        charge_sum = 0 * unit.elementary_charge
        abs_charge_sum = 0 * unit.elementary_charge
        for pc in molecule._partial_charges:
            charge_sum += pc
            abs_charge_sum += abs(pc)
        assert abs(charge_sum) < 0.005 * unit.elementary_charge
        assert abs_charge_sum > 0.25 * unit.elementary_charge

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_compute_partial_charges_am1bcc_net_charge(self):
        """Test OpenEyeToolkitWrapper assign_partial_charges() on a molecule with a net +1 charge"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_acetate()
        molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry)
        charge_sum = 0 * unit.elementary_charge
        for pc in molecule._partial_charges:
            charge_sum += pc
        assert -0.999 * unit.elementary_charge > charge_sum > -1.001 * unit.elementary_charge

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_compute_partial_charges_am1bcc_wrong_n_confs(self):
        """
        Test OpenEyeToolkitWrapper compute_partial_charges_am1bcc() when requesting to use an incorrect number of
        conformers. This test is a bit shorter than that for AmberToolsToolkitWrapper because OETK uses the
        ELF10 multiconformer method of AM1BCC, which doesn't have a maximum number of conformers.
        """
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=2,
                                     rms_cutoff=0.1*unit.angstrom,
                                     toolkit_registry=toolkit_registry)

        # Try again, with strict_n_confs as true, but not including use_confs, so the
        # recommended number of confs will be generated
        molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry,
                                                strict_n_conformers=True)



    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    @pytest.mark.parametrize("partial_charge_method", ['am1bcc', 'am1-mulliken', 'gasteiger'])
    def test_assign_partial_charges_neutral(self, partial_charge_method):
        """Test OpenEyeToolkitWrapper assign_partial_charges()"""
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method)
        charge_sum = 0. * unit.elementary_charge
        for pc in molecule.partial_charges:
            charge_sum += pc
        assert -1.e-5 < charge_sum.value_in_unit(unit.elementary_charge) < 1.e-5

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    @pytest.mark.parametrize("partial_charge_method", ['am1bcc', 'am1-mulliken'])
    def test_assign_partial_charges_conformer_dependence(self, partial_charge_method):
        """Test OpenEyeToolkitWrapper assign_partial_charges()'s use_conformers kwarg
        to ensure charges are really conformer dependent. Skip Gasteiger because it isn't
        conformer dependent."""
        from openforcefield.tests.test_forcefield import create_ethanol
        import copy
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=1)
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method,
                                        use_conformers=molecule.conformers)
        pcs1 = copy.deepcopy(molecule.partial_charges)
        molecule._conformers[0][0][0] += 0.2 * unit.angstrom
        molecule._conformers[0][1][1] -= 0.2 * unit.angstrom
        molecule._conformers[0][2][1] += 0.2 * unit.angstrom
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method,
                                        use_conformers=molecule.conformers)
        for pc1, pc2 in zip(pcs1, molecule.partial_charges):
            assert abs(pc1 - pc2) > 1.e-5 * unit.elementary_charge

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    @pytest.mark.parametrize("partial_charge_method", ['am1bcc', 'am1-mulliken', 'gasteiger'])
    def test_assign_partial_charges_net_charge(self, partial_charge_method):
        """
        Test OpenEyeToolkitWrapper assign_partial_charges() on a molecule with net charge.
        """
        from openforcefield.tests.test_forcefield import create_acetate
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_acetate()
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method)
        charge_sum = 0. * unit.elementary_charge
        for pc in molecule.partial_charges:
            charge_sum += pc
        assert -1.e-5 < charge_sum.value_in_unit(unit.elementary_charge) + 1. < 1.e-5

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_assign_partial_charges_bad_charge_method(self):
        """Test OpenEyeToolkitWrapper assign_partial_charges() for a nonexistent charge method"""
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()

        # Molecule.assign_partial_charges calls the ToolkitRegistry with raise_exception_types = [],
        # which means it will only ever return ValueError
        with pytest.raises(ValueError, match="is not available from OpenEyeToolkitWrapper") as excinfo:
            molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                            partial_charge_method="NotARealChargeMethod")

        # ToolkitWrappers raise a specific exception class, so we test that here
        with pytest.raises(ChargeMethodUnavailableError, match="is not available from OpenEyeToolkitWrapper") as excinfo:
            OETKW = OpenEyeToolkitWrapper()
            OETKW.assign_partial_charges(molecule=molecule,
                                         partial_charge_method="NotARealChargeMethod")

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    @pytest.mark.parametrize("partial_charge_method,expected_n_confs", [('am1bcc', 1),
                                                                        ('am1-mulliken', 1),
                                                                        ('gasteiger', 0)])
    def test_assign_partial_charges_wrong_n_confs(self, partial_charge_method, expected_n_confs):
        """
        Test OpenEyeToolkitWrapper assign_partial_charges() when requesting to use an incorrect number of
        conformers
        """
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[OpenEyeToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=2, rms_cutoff=0.01*unit.angstrom)

        # Try passing in the incorrect number of confs, but without specifying strict_n_conformers,
        # which should produce a warning
        with pytest.warns(IncorrectNumConformersWarning,
                          match=f"has 2 conformers, but charge method '{partial_charge_method}' "
                                f"expects exactly {expected_n_confs}."):
            molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                            partial_charge_method=partial_charge_method,
                                            use_conformers=molecule.conformers,
                                            strict_n_conformers=False)

        # Try again, with strict_n_confs as true, but not including use_confs, so the
        # recommended number of confs will be generated
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method,
                                        strict_n_conformers=True)

        # Test calling the ToolkitWrapper _indirectly_, though a ToolkitRegistry,
        # which should aggregate any exceptions and bundle all of the messages
        # in a failed task together in a single ValueError.
        with pytest.raises(ValueError,
                           match=f"has 2 conformers, but charge method '{partial_charge_method}' "
                                 f"expects exactly {expected_n_confs}."):
            molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                            partial_charge_method=partial_charge_method,
                                            use_conformers=molecule.conformers,
                                            strict_n_conformers=True)

        # Test calling the ToolkitWrapper _directly_, passing in the incorrect number of
        # confs, and specify strict_n_conformers, which should produce an IncorrectNumConformersError
        with pytest.raises(IncorrectNumConformersError,
                           match=f"has 2 conformers, but charge method '{partial_charge_method}' "
                                 f"expects exactly {expected_n_confs}."):
            OETKW = OpenEyeToolkitWrapper()
            OETKW.assign_partial_charges(molecule=molecule,
                                         partial_charge_method=partial_charge_method,
                                         use_conformers=molecule.conformers,
                                         strict_n_conformers=True)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_compute_partial_charges_failure(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges() on a molecule it cannot assign charges to"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[Li+1]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)

        # For now, I'm just testing AM1-BCC (will test more when the SMIRNOFF spec for other charges is finalized)
        with pytest.raises(Exception) as excinfo:
            molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_wrapper)
            assert "Unable to assign charges" in str(excinfo)
            assert "OE Error: " in str(excinfo)


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_compute_partial_charges_trans_cooh_am1bcc(self):
        """Test OpenEyeToolkitWrapper for computing partial charges for problematic molecules, as exemplified by
        Issue 346 (https://github.com/openforcefield/openforcefield/issues/346)"""

        lysine = Molecule.from_smiles("C(CC[NH3+])C[C@@H](C(=O)O)N")
        toolkit_wrapper = OpenEyeToolkitWrapper()
        lysine.generate_conformers(toolkit_registry=toolkit_wrapper)
        lysine.compute_partial_charges_am1bcc(toolkit_registry=toolkit_wrapper)


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_assign_fractional_bond_orders(self):
        """Test OpenEyeToolkitWrapper assign_fractional_bond_orders()"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        for bond_order_model in ['am1-wiberg', 'pm3-wiberg']:
            molecule.assign_fractional_bond_orders(toolkit_registry=toolkit_wrapper,
                                                    bond_order_model=bond_order_model)
            # TODO: Add test for equivalent Wiberg orders for equivalent bonds



    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_assign_fractional_bond_orders_neutral_charge_mol(self):
        """Test OpenEyeToolkitWrapper assign_fractional_bond_orders() for neutral and charged molecule"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        # Reading neutral molecule from file
        filename = get_data_file_path('molecules/CID20742535_neutral.sdf')
        molecule1 = Molecule.from_file(filename)
        # Reading negative molecule from file
        filename = get_data_file_path('molecules/CID20742535_anion.sdf')
        molecule2 = Molecule.from_file(filename)

        # Checking that only one additional bond is present in the neutral molecule
        assert (len(molecule1.bonds) == len(molecule2.bonds)+1)

        for bond_order_model in ['am1-wiberg']:
            molecule1.assign_fractional_bond_orders(toolkit_registry=toolkit_wrapper,
                                                    bond_order_model=bond_order_model,
                                                    use_conformers=molecule1.conformers)

            for i in molecule1.bonds:
                if i.is_aromatic:
                    # Checking aromatic bonds
                    assert (1.05 < i.fractional_bond_order < 1.65)
                elif (i.atom1.atomic_number == 1 or i.atom2.atomic_number == 1):
                    # Checking bond order of C-H or O-H bonds are around 1
                    assert (0.85 < i.fractional_bond_order < 1.05)
                elif (i.atom1.atomic_number == 8 or i.atom2.atomic_number == 8):
                    # Checking C-O single bond
                    wbo_C_O_neutral = i.fractional_bond_order
                    assert (1.0 < wbo_C_O_neutral < 1.5)
                else:
                    # Should be C-C single bond
                    assert (i.atom1_index == 4 and i.atom2_index == 6) or (i.atom1_index == 6 and i.atom2_index == 4)
                    wbo_C_C_neutral = i.fractional_bond_order
                    assert (1.0 < wbo_C_C_neutral < 1.3)

            molecule2.assign_fractional_bond_orders(toolkit_registry=toolkit_wrapper,
                                                    bond_order_model=bond_order_model,
                                                    use_conformers=molecule2.conformers)
            for i in molecule2.bonds:
                if i.is_aromatic:
                    # Checking aromatic bonds
                    assert (1.05 < i.fractional_bond_order < 1.65)
                elif (i.atom1.atomic_number == 1 or i.atom2.atomic_number == 1):
                    # Checking bond order of C-H or O-H bonds are around 1
                    assert (0.85 < i.fractional_bond_order < 1.05)
                elif (i.atom1.atomic_number == 8 or i.atom2.atomic_number == 8):
                    # Checking C-O single bond
                    wbo_C_O_anion = i.fractional_bond_order
                    assert (1.3 < wbo_C_O_anion < 1.8)
                else:
                    # Should be C-C single bond
                    assert(i.atom1_index == 4 and i.atom2_index == 6) or (i.atom1_index == 6 and i.atom2_index == 4)
                    wbo_C_C_anion = i.fractional_bond_order
                    assert (1.0 < wbo_C_C_anion < 1.3)

            # Wiberg bond order of C-C single bond is higher in the anion
            assert (wbo_C_C_anion > wbo_C_C_neutral)
            # Wiberg bond order of C-O bond is higher in the anion
            assert (wbo_C_O_anion > wbo_C_O_neutral)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_assign_fractional_bond_orders_charged(self):
        """Test OpenEyeToolkitWrapper assign_fractional_bond_orders() on a molecule with net charge +1"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        for bond_order_model in ['am1-wiberg', 'pm3-wiberg']:
            molecule.assign_fractional_bond_orders(toolkit_registry=toolkit_wrapper,
                                                   bond_order_model=bond_order_model)
            # TODO: Add test for equivalent Wiberg orders for equivalent bonds

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_assign_fractional_bond_orders_invalid_method(self):
        """
        Test that OpenEyeToolkitWrapper assign_fractional_bond_orders() raises the
        correct error if an invalid charge model is provided
        """
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        expected_error = "Bond order model 'not a real bond order model' is not supported by " \
                         "OpenEyeToolkitWrapper. Supported models are ([[]'am1-wiberg', 'pm3-wiberg'[]])"
        with pytest.raises(ValueError, match=expected_error) as excinfo:
            molecule.assign_fractional_bond_orders(toolkit_registry=toolkit_wrapper,
                                                    bond_order_model='not a real bond order model')


    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_assign_fractional_bond_orders_double_bond(self):
        """Test OpenEyeToolkitWrapper assign_fractional_bond_orders() on a molecule with a double bond"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = r'C\C(F)=C(/F)C[C@@](C)(Cl)Br'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        for bond_order_model in ['am1-wiberg', 'pm3-wiberg']:
            molecule.assign_fractional_bond_orders(toolkit_registry=toolkit_wrapper,
                                                   bond_order_model=bond_order_model)
            # TODO: Add test for equivalent Wiberg orders for equivalent bonds

        double_bond_has_wbo_near_2 = False
        for bond in molecule.bonds:
            if bond.bond_order == 2:
                if 1.75 < bond.fractional_bond_order < 2.25:
                    double_bond_has_wbo_near_2 = True
        assert double_bond_has_wbo_near_2

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_substructure_search_on_large_molecule(self):
        """Test OpenEyeToolkitWrapper substructure search when a large number hits are found"""

        tk = OpenEyeToolkitWrapper()
        smiles = "C"*600
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
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups='[#6:1]-[#8:2]')
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 6

        # now ignore the O-C bond, backwards
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups='[#8:1]-[#6:2]')
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 6

        # now ignore the C-C bond
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups='[#6:1]-[#6:2]')
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 8

        # ignore a list of searches, forward
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups=['[#6:1]-[#8:2]', '[#6:1]-[#6:2]'])
        assert bonds == []

        # ignore a list of searches, backwards
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups=['[#6:1]-[#6:2]', '[#8:1]-[#6:2]'])
        assert bonds == []

        # test  molecules that should have no rotatable bonds
        cyclohexane = create_cyclohexane()
        bonds = cyclohexane.find_rotatable_bonds()
        assert bonds == []

        methane = Molecule.from_smiles('C')
        bonds = methane.find_rotatable_bonds()
        assert bonds == []

        ethene = Molecule.from_smiles('C=C')
        bonds = ethene.find_rotatable_bonds()
        assert bonds == []

        terminal_forwards = '[*]~[*:1]-[X2H1,X3H2,X4H3:2]-[#1]'
        terminal_backwards = '[#1]-[X2H1,X3H2,X4H3:1]-[*:2]~[*]'
        # test removing terminal rotors
        toluene = Molecule.from_file(get_data_file_path('molecules/toluene.sdf'))
        bonds = toluene.find_rotatable_bonds()
        assert len(bonds) == 1
        assert toluene.atoms[bonds[0].atom1_index].atomic_number == 6
        assert toluene.atoms[bonds[0].atom2_index].atomic_number == 6

        # find terminal bonds forward
        bonds = toluene.find_rotatable_bonds(ignore_functional_groups=terminal_forwards)
        assert bonds == []

        # find terminal bonds backwards
        bonds = toluene.find_rotatable_bonds(ignore_functional_groups=terminal_backwards)
        assert bonds == []
        
        
        # TODO: Check partial charge invariants (total charge, charge equivalence)

        # TODO: Add test for aromaticity
        # TODO: Add test and molecule functionality for isotopes



class TestRDKitToolkitWrapper:
    """Test the RDKitToolkitWrapper"""
    
    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_smiles(self):
        """Test RDKitToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        # This differs from OE's expected output due to different canonicalization schemes
        smiles = '[H][C]([H])([H])[C]([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles,
                                        toolkit_registry=toolkit_wrapper)
        # When making a molecule from SMILES, partial charges should be initialized to None
        assert molecule.partial_charges is None
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        #print(smiles, smiles2)
        assert smiles == smiles2

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    @pytest.mark.parametrize("smiles,exception_regex", [
        (r"C\C(F)=C(/F)CC(C)(Cl)Br", "Undefined chiral centers"),
        (r"C\C(F)=C(/F)C[C@@](C)(Cl)Br", None),
        (r"CC(F)=C(F)C[C@@](C)(Cl)Br", "Bonds with undefined stereochemistry")
    ])
    def test_smiles_missing_stereochemistry(self, smiles, exception_regex):
        """Test RDKitToolkitWrapper to_smiles() and from_smiles() when given ambiguous stereochemistry"""
        toolkit_wrapper = RDKitToolkitWrapper()

        if exception_regex is not None:
            with pytest.raises(UndefinedStereochemistryError, match=exception_regex):
                Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)
            Molecule.from_smiles(smiles,
                                 toolkit_registry=toolkit_wrapper,
                                 allow_undefined_stereo=True)
        else:
            Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)

    # TODO: test_smiles_round_trip

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_smiles_add_H(self):
        """Test RDKitToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        input_smiles = 'CC'
        # This differs from OE's expected output due to different canonicalization schemes
        expected_output_smiles = '[H][C]([H])([H])[C]([H])([H])[H]'
        molecule = Molecule.from_smiles(input_smiles,
                                        toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles2 == expected_output_smiles

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_rdkit_from_smiles_hydrogens_are_explicit(self):
        """
        Test to ensure that RDKitToolkitWrapper.from_smiles has the proper behavior with
        respect to its hydrogens_are_explicit kwarg
        """
        toolkit_wrapper = RDKitToolkitWrapper()
        smiles_impl = "C#C"
        with pytest.raises(ValueError,
                           match="but RDKit toolkit interpreted SMILES 'C#C' as having implicit hydrogen") as excinfo:
            offmol = Molecule.from_smiles(smiles_impl,
                                          toolkit_registry=toolkit_wrapper,
                                          hydrogens_are_explicit=True)
        offmol = Molecule.from_smiles(smiles_impl,
                                      toolkit_registry=toolkit_wrapper,
                                      hydrogens_are_explicit=False)
        assert offmol.n_atoms == 4

        smiles_expl = "[H][C]#[C][H]"
        offmol = Molecule.from_smiles(smiles_expl,
                                      toolkit_registry=toolkit_wrapper,
                                      hydrogens_are_explicit=True)
        assert offmol.n_atoms == 4
        # It's debatable whether this next function should pass. Strictly speaking, the hydrogens in this SMILES
        # _are_ explicit, so allowing "hydrogens_are_explicit=False" through here is allowing a contradiction.
        # We might rethink the name of this kwarg.

        offmol = Molecule.from_smiles(smiles_expl,
                                      toolkit_registry=toolkit_wrapper,
                                      hydrogens_are_explicit=False)
        assert offmol.n_atoms == 4

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    @pytest.mark.parametrize('molecule', get_mini_drug_bank(RDKitToolkitWrapper))
    def test_to_inchi(self, molecule):
        """Test conversion to standard and non-standard InChI"""

        toolkit = RDKitToolkitWrapper()
        inchi = molecule.to_inchi(toolkit_registry=toolkit)
        non_standard = molecule.to_inchi(fixed_hydrogens=True,toolkit_registry=toolkit)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    @pytest.mark.parametrize('molecule', get_mini_drug_bank(RDKitToolkitWrapper))
    def test_to_inchikey(self, molecule):
        """Test the conversion to standard and non-standard InChIKey"""

        toolkit = RDKitToolkitWrapper()
        inchikey = molecule.to_inchikey(toolkit_registry=toolkit)
        non_standard_key = molecule.to_inchikey(fixed_hydrogens=True, toolkit_registry=toolkit)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_from_bad_inchi(self):
        """Test building a molecule from a bad InChI string"""

        toolkit = RDKitToolkitWrapper()
        inchi = 'InChI=1S/ksbfksfksfksbfks'
        with pytest.raises(RuntimeError):
            mol = Molecule.from_inchi(inchi, toolkit_registry=toolkit)

    inchi_data = [{'molecule': create_ethanol(), 'standard_inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
                   'fixed_hydrogen_inchi': 'InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3'},
                  {'molecule': create_reversed_ethanol(), 'standard_inchi': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
                   'fixed_hydrogen_inchi': 'InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3'},
                  {'molecule': create_acetaldehyde(), 'standard_inchi': 'InChI=1S/C2H4O/c1-2-3/h2H,1H3',
                   'fixed_hydrogen_inchi': 'InChI=1/C2H4O/c1-2-3/h2H,1H3'},
                  {'molecule': create_cyclohexane(), 'standard_inchi': 'InChI=1S/C6H12/c1-2-4-6-5-3-1/h1-6H2',
                   'fixed_hydrogen_inchi': 'InChI=1/C6H12/c1-2-4-6-5-3-1/h1-6H2'}
                  ]

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    @pytest.mark.parametrize('data', inchi_data)
    def test_from_inchi(self, data):
        """Test building a molecule from standard and non-standard InChI strings."""

        toolkit = RDKitToolkitWrapper()

        ref_mol = data['molecule']
        # make a molecule from inchi
        inchi_mol = Molecule.from_inchi(data['standard_inchi'], toolkit_registry=toolkit)
        assert inchi_mol.to_inchi(toolkit_registry=toolkit) == data['standard_inchi']

        def compare_mols(ref_mol, inchi_mol):
            assert ref_mol.n_atoms == inchi_mol.n_atoms
            assert ref_mol.n_bonds == inchi_mol.n_bonds
            assert ref_mol.n_angles == inchi_mol.n_angles
            assert ref_mol.n_propers == inchi_mol.n_propers
            assert ref_mol.is_isomorphic_with(inchi_mol) is True

        compare_mols(ref_mol, inchi_mol)

        # now make the molecule from the non-standard inchi and compare
        nonstandard_inchi_mol = Molecule.from_inchi(data['fixed_hydrogen_inchi'])
        assert nonstandard_inchi_mol.to_inchi(fixed_hydrogens=True, toolkit_registry=toolkit) == data['fixed_hydrogen_inchi']

        compare_mols(ref_mol, nonstandard_inchi_mol)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    @pytest.mark.parametrize('molecule', get_mini_drug_bank(RDKitToolkitWrapper))
    def test_non_standard_inchi_round_trip(self, molecule):
        """Test if a molecule can survive an InChi round trip test in some cases the standard InChI
        will not be enough to ensure information is preserved so we test the non-standard inchi here."""

        from openforcefield.utils.toolkits import UndefinedStereochemistryError

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
            if molecule.name in rdkit_inchi_isomorphic_fails:
                # Some molecules graphs change during the round trip testing
                # we test quite strict isomorphism here
                with pytest.raises(AssertionError):
                    assert molecule.is_isomorphic_with(mol2, bond_order_matching=False)
            else:
                assert molecule.is_isomorphic_with(mol2, bond_order_matching=False)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_smiles_charged(self):
        """Test RDKitWrapper functions for reading/writing charged SMILES"""
        toolkit_wrapper = RDKitToolkitWrapper()
        # This differs from OE's expected output due to different canonicalization schemes
        smiles = '[H][C]([H])([H])[N+]([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles,
                                        toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_to_from_rdkit_core_props_filled(self):
        """Test RDKitToolkitWrapper to_rdkit() and from_rdkit() when given populated core property fields"""
        toolkit_wrapper = RDKitToolkitWrapper()

        # Replacing with a simple molecule with stereochemistry
        input_smiles = r'C\C(F)=C(/F)C[C@@](C)(Cl)Br'
        expected_output_smiles = r'[H][C]([H])([H])/[C]([F])=[C](\[F])[C]([H])([H])[C@@]([Cl])([Br])[C]([H])([H])[H]'
        molecule = Molecule.from_smiles(input_smiles, toolkit_registry=toolkit_wrapper)
        assert molecule.to_smiles(toolkit_registry=toolkit_wrapper) == expected_output_smiles

        # Populate core molecule property fields
        molecule.name = 'Alice'
        partial_charges = unit.Quantity(np.array([-.9, -.8, -.7, -.6,
                                                  -.5, -.4, -.3, -.2,
                                                  -.1,  0.,  .1,  .2,
                                                   .3,  .4,  .5,  .6,
                                                   .7,  .8]), unit.elementary_charge)
        molecule.partial_charges = partial_charges
        coords = unit.Quantity(np.array([['0.0', '1.0', '2.0'],    ['3.0', '4.0', '5.0'],    ['6.0', '7.0', '8.0'],
                                         ['9.0', '10.0', '11.0'] , ['12.0', '13.0', '14.0'], ['15.0', '16.0', '17.0'],
                                         ['18.0', '19.0', '20.0'], ['21.0', '22.0', '23.0'], ['24.0', '25.0', '26.0'],
                                         ['27.0', '28.0', '29.0'], ['30.0', '31.0', '32.0'], ['33.0', '34.0', '35.0'],
                                         ['36.0', '37.0', '38.0'], ['39.0', '40.0', '41.0'], ['42.0', '43.0', '44.0'],
                                         ['45.0', '46.0', '47.0'], ['48.0', '49.0', '50.0'], ['51.0', '52.0', '53.0']]),
                                    unit.angstrom)
        molecule.add_conformer(coords)
        # Populate core atom property fields
        molecule.atoms[2].name = 'Bob'
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

        # Do a first conversion to/from oemol
        rdmol = molecule.to_rdkit()
        molecule2 = Molecule.from_rdkit(rdmol)

        # Test that properties survived first conversion
        #assert molecule.to_dict() == molecule2.to_dict()
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
            pc1_ul = pc1 / unit.elementary_charge
            pc2_ul = pc2 / unit.elementary_charge
            assert_almost_equal(pc1_ul, pc2_ul, decimal=6)
        assert molecule2.to_smiles(toolkit_registry=toolkit_wrapper) == expected_output_smiles
        # TODO: This should be its own test

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_to_from_rdkit_core_props_unset(self):
        """Test RDKitToolkitWrapper to_rdkit() and from_rdkit() when given empty core property fields"""
        toolkit_wrapper = RDKitToolkitWrapper()

        # Replacing with a simple molecule with stereochemistry
        input_smiles = r'C\C(F)=C(/F)C[C@](C)(Cl)Br'
        expected_output_smiles = r'[H][C]([H])([H])/[C]([F])=[C](\[F])[C]([H])([H])[C@]([Cl])([Br])[C]([H])([H])[H]'
        molecule = Molecule.from_smiles(input_smiles, toolkit_registry=toolkit_wrapper)
        assert molecule.to_smiles(toolkit_registry=toolkit_wrapper) == expected_output_smiles

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

        assert molecule2.to_smiles(toolkit_registry=toolkit_wrapper) == expected_output_smiles
        
    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_get_sdf_coordinates(self):
        """Test RDKitToolkitWrapper for importing a single set of coordinates from a sdf file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule.conformers) == 1
        assert molecule.conformers[0].shape == (15, 3)
        assert_almost_equal(molecule.conformers[0][5][1] / unit.angstrom, 2.0104, decimal=4)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_read_sdf_charges(self):
        """Test RDKitToolkitWrapper for importing a charges from a sdf file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/ethanol_partial_charges.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert molecule.partial_charges is not None
        assert molecule.partial_charges[0] == -0.4 * unit.elementary_charge
        assert molecule.partial_charges[-1] == 0.4 * unit.elementary_charge

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_write_sdf_charges(self):
        """Test RDKitToolkitWrapper for writing partial charges to a sdf file"""
        from io import StringIO
        toolkit_wrapper = RDKitToolkitWrapper()
        ethanol = create_ethanol()
        sio = StringIO()
        ethanol.to_file(sio, 'SDF', toolkit_registry=toolkit_wrapper)
        sdf_text = sio.getvalue()
        # The output lines of interest here will look like
        # >  <atom.dprop.PartialCharge>  (1)
        # -0.40000000000000002 -0.29999999999999999 -0.20000000000000001 -0.10000000000000001 0.01 0.10000000000000001 0.20000000000000001 0.29999999999999999 0.40000000000000002

        # Parse the SDF text, grabbing the numeric line above
        sdf_split = sdf_text.split('\n')
        charge_line_found = False
        for line in sdf_split:
            if charge_line_found:
                charges = [float(i) for i in line.split()]
                break
            if '>  <atom.dprop.PartialCharge>' in line:
                charge_line_found = True

        # Make sure that a charge line was ever found
        assert charge_line_found

        # Make sure that the charges found were correct
        assert_almost_equal(charges, [-0.4, -0.3, -0.2, -0.1, 0.00001, 0.1, 0.2, 0.3, 0.4])


    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_sdf_properties_roundtrip(self):
        """Test RDKitToolkitWrapper for performing a round trip of a molecule with defined partial charges
        and entries in the properties dict to and from a sdf file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        ethanol = create_ethanol()
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix='.sdf') as iofile:
            ethanol.to_file(iofile.name, file_format='SDF', toolkit_registry=toolkit_wrapper)
            ethanol2 = Molecule.from_file(iofile.name, file_format='SDF', toolkit_registry=toolkit_wrapper)
        assert (ethanol.partial_charges == ethanol2.partial_charges).all()

        # Now test with no properties or charges
        ethanol = create_ethanol()
        ethanol.partial_charges = None
        # Write ethanol to a temporary file, and then immediately read it.
        with NamedTemporaryFile(suffix='.sdf') as iofile:
            ethanol.to_file(iofile.name, file_format='SDF', toolkit_registry=toolkit_wrapper)
            ethanol2 = Molecule.from_file(iofile.name, file_format='SDF', toolkit_registry=toolkit_wrapper)
        assert ethanol2.partial_charges is None
        assert ethanol2.properties == {}


    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_write_sdf_no_charges(self):
        """Test RDKitToolkitWrapper for writing an SDF file with no charges"""
        from io import StringIO
        toolkit_wrapper = RDKitToolkitWrapper()
        ethanol = create_ethanol()
        ethanol.partial_charges = None
        sio = StringIO()
        ethanol.to_file(sio, 'SDF', toolkit_registry=toolkit_wrapper)
        sdf_text = sio.getvalue()
        # In our current configuration, if the OFFMol doesn't have partial charges, we DO NOT want a partial charge
        # block to be written. For reference, it's possible to indicate that a partial charge is not known by writing
        # out "n/a" (or another placeholder) in the partial charge block atoms without charges.
        assert '>  <atom.dprop.PartialCharge>' not in sdf_text


    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_load_multiconformer_sdf_as_separate_molecules(self):
        """
        Test RDKitToolkitWrapper for reading a "multiconformer" SDF, which the OFF
        Toolkit should treat as separate molecules
        """
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/methane_multiconformer.sdf')
        molecules = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecules) == 2
        assert len(molecules[0].conformers) == 1
        assert len(molecules[1].conformers) == 1
        assert molecules[0].conformers[0].shape == (5, 3)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_load_multiconformer_sdf_as_separate_molecules_properties(self):
        """
        Test RDKitToolkitWrapper for reading a "multiconformer" SDF, which the OFF
        Toolkit should treat as separate molecules
        """
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/methane_multiconformer_properties.sdf')
        molecules = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecules) == 2
        assert len(molecules[0].conformers) == 1
        assert len(molecules[1].conformers) == 1
        assert molecules[0].conformers[0].shape == (5, 3)
        # The first molecule in the SDF has the following properties and charges:
        assert molecules[0].properties['test_property_key'] == 'test_property_value'
        np.testing.assert_allclose(molecules[0].partial_charges / unit.elementary_charge,
                                          [-0.108680, 0.027170, 0.027170, 0.027170, 0.027170])
        # The second molecule in the SDF has the following properties and charges:
        assert molecules[1].properties['test_property_key'] == 'test_property_value2'
        assert molecules[1].properties['another_test_property_key'] == 'another_test_property_value'
        np.testing.assert_allclose(molecules[1].partial_charges / unit.elementary_charge,
                                   [0.027170, 0.027170, 0.027170, 0.027170, -0.108680])

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_write_multiconformer_mol_as_sdf(self):
        """
        Test RDKitToolkitWrapper for writing a multiconformer molecule to SDF. The OFF toolkit should only
        save the first conformer
        """
        from io import StringIO

        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/ethanol.sdf')
        ethanol = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        ethanol.partial_charges = np.array([-4., -3., -2., -1., 0., 1., 2., 3., 4.]) * unit.elementary_charge
        ethanol.properties['test_prop'] = 'test_value'
        new_conf = ethanol.conformers[0] + (np.ones(ethanol.conformers[0].shape) * unit.angstrom)
        ethanol.add_conformer(new_conf)
        sio = StringIO()
        ethanol.to_file(sio, 'sdf', toolkit_registry=toolkit_wrapper)
        data = sio.getvalue()
        # In SD format, each molecule ends with "$$$$"
        assert data.count('$$$$') == 1
        # A basic SDF for ethanol would be 27 lines, though the properties add three more
        assert len(data.split('\n')) == 30
        assert 'test_prop' in data
        assert '<atom.dprop.PartialCharge>' in data
        # Ensure the first conformer's first atom's X coordinate is in the file
        assert str(ethanol.conformers[0][0][0].value_in_unit(unit.angstrom))[:5] in data
        # Ensure the SECOND conformer's first atom's X coordinate is NOT in the file
        assert str(ethanol.conformers[1][0][0].in_units_of(unit.angstrom))[:5] not in data

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_write_milticonformer_pdb(self):
        """
        Make sure RDKit can write multi conformer PDB files.
        """
        from io import StringIO

        toolkit = RDKitToolkitWrapper()
        # load up a multiconformer pdb file and condense down the conformers
        molecules = Molecule.from_file(get_data_file_path('molecules/butane_multi.sdf'), toolkit_registry=toolkit)
        butane = molecules.pop(0)
        for mol in molecules:
            butane.add_conformer(mol.conformers[0])
        assert butane.n_conformers == 7
        sio = StringIO()
        butane.to_file(sio, 'pdb', toolkit_registry=toolkit)
        # we need to make sure each conformer is wrote to the file
        pdb = sio.getvalue()
        for i in range(1, 8):
            assert f'MODEL        {i}' in pdb

    # Unskip this when we implement PDB-reading support for RDKitToolkitWrapper
    @pytest.mark.skip
    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_get_pdb_coordinates(self):
        """Test RDKitToolkitWrapper for importing a single set of coordinates from a pdb file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.pdb')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule.conformers) == 1
        assert molecule.conformers[0].shape == (15,3)

    # Unskip this when we implement PDB-reading support for RDKitToolkitWrapper
    @pytest.mark.skip
    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_load_aromatic_pdb(self):
        """Test OpenEyeToolkitWrapper for importing molecule conformers"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.pdb')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule.conformers) == 1
        assert molecule.conformers[0].shape == (15,3)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_generate_conformers(self):
        """Test RDKitToolkitWrapper generate_conformers()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers()
        # TODO: Make this test more robust

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_generate_multiple_conformers(self):
        """Test RDKitToolkitWrapper generate_conformers() for generating multiple conformers"""
        toolkit_wrapper = RDKitToolkitWrapper()
        smiles = 'CCCCCCC'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(rms_cutoff=1 * unit.angstrom,
                                     n_conformers=100,
                                     toolkit_registry=toolkit_wrapper)
        assert molecule.n_conformers > 1
        assert not (molecule.conformers[0] == (0. * unit.angstrom)).all()

        # Ensure rms_cutoff kwarg is working
        molecule2 = toolkit_wrapper.from_smiles(smiles)
        molecule2.generate_conformers(rms_cutoff=0.1 * unit.angstrom,
                                      n_conformers=100,
                                      toolkit_registry=toolkit_wrapper)
        assert molecule2.n_conformers > molecule.n_conformers

        # Ensure n_conformers kwarg is working
        molecule2 = toolkit_wrapper.from_smiles(smiles)
        molecule2.generate_conformers(rms_cutoff=0.1 * unit.angstrom,
                                      n_conformers=10,
                                      toolkit_registry=toolkit_wrapper)
        assert molecule2.n_conformers == 10

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
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
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups='[#6:1]-[#8:2]')
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 6

        # now ignore the O-C bond, backwards
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups='[#8:1]-[#6:2]')
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 6

        # now ignore the C-C bond
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups='[#6:1]-[#6:2]')
        assert len(bonds) == 1
        assert ethanol.atoms[bonds[0].atom1_index].atomic_number == 6
        assert ethanol.atoms[bonds[0].atom2_index].atomic_number == 8

        # ignore a list of searches, forward
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups=['[#6:1]-[#8:2]', '[#6:1]-[#6:2]'])
        assert bonds == []

        # ignore a list of searches, backwards
        bonds = ethanol.find_rotatable_bonds(ignore_functional_groups=['[#6:1]-[#6:2]', '[#8:1]-[#6:2]'])
        assert bonds == []

        # test  molecules that should have no rotatable bonds
        cyclohexane = create_cyclohexane()
        bonds = cyclohexane.find_rotatable_bonds()
        assert bonds == []

        methane = Molecule.from_smiles('C')
        bonds = methane.find_rotatable_bonds()
        assert bonds == []

        ethene = Molecule.from_smiles('C=C')
        bonds = ethene.find_rotatable_bonds()
        assert bonds == []

        terminal_forwards = '[*]~[*:1]-[X2H1,X3H2,X4H3:2]-[#1]'
        terminal_backwards = '[#1]-[X2H1,X3H2,X4H3:1]-[*:2]~[*]'
        # test removing terminal rotors
        toluene = Molecule.from_file(get_data_file_path('molecules/toluene.sdf'))
        bonds = toluene.find_rotatable_bonds()
        assert len(bonds) == 1
        assert toluene.atoms[bonds[0].atom1_index].atomic_number == 6
        assert toluene.atoms[bonds[0].atom2_index].atomic_number == 6

        # find terminal bonds forward
        bonds = toluene.find_rotatable_bonds(ignore_functional_groups=terminal_forwards)
        assert bonds == []

        # find terminal bonds backwards
        bonds = toluene.find_rotatable_bonds(ignore_functional_groups=terminal_backwards)
        assert bonds == []

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_to_rdkit_losing_aromaticity_(self):
        # test the example given in issue #513
        # <https://github.com/openforcefield/openforcefield/issues/513>
        smiles = "[H]c1c(c(c(c(c1OC2=C(C(=C(N3C2=C(C(=C3[H])C#N)[H])[H])F)[H])OC([H])([H])C([H])([H])N4C(=C(C(=O)N(C4=O)[H])[H])[H])[H])F)[H]"

        mol = Molecule.from_smiles(smiles)
        rdmol = mol.to_rdkit()

        # now make sure the aromaticity matches for each atom
        for (offatom, rdatom) in zip(mol.atoms, rdmol.GetAtoms()):
            assert offatom.is_aromatic is rdatom.GetIsAromatic()

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_substructure_search_on_large_molecule(self):
        """Test RDKitToolkitWrapper substructure search when a large number hits are found"""

        tk = RDKitToolkitWrapper()
        smiles = "C"*3000
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




class TestAmberToolsToolkitWrapper:
    """Test the AmberToolsToolkitWrapper"""

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                    reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_compute_partial_charges_am1bcc(self):
        """Test AmberToolsToolkitWrapper compute_partial_charges_am1bcc()"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        molecule = create_ethanol()
        molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry)
        charge_sum = 0 * unit.elementary_charge
        abs_charge_sum = 0 * unit.elementary_charge
        for pc in molecule._partial_charges:
            charge_sum += pc
            abs_charge_sum += abs(pc)
        assert abs(charge_sum) < 0.001 * unit.elementary_charge
        assert abs_charge_sum > 0.25 * unit.elementary_charge

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_compute_partial_charges_am1bcc_net_charge(self):
        """Test AmberToolsToolkitWrapper assign_partial_charges() on a molecule with a net -1 charge"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        molecule = create_acetate()
        molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry)
        charge_sum = 0 * unit.elementary_charge
        for pc in molecule._partial_charges:
            charge_sum += pc
        assert -0.99 * unit.elementary_charge > charge_sum > -1.01 * unit.elementary_charge

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_compute_partial_charges_am1bcc_wrong_n_confs(self):
        """
        Test AmberToolsToolkitWrapper compute_partial_charges_am1bcc() when requesting to use an incorrect number of
        conformers
        """
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=2, rms_cutoff=0.1*unit.angstrom)

        # Try passing in the incorrect number of confs, but without specifying strict_n_conformers,
        # which should produce a warning
        with pytest.warns(IncorrectNumConformersWarning,
                          match="has 2 conformers, but charge method 'am1bcc' expects exactly 1."):
            molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry,
                                                    use_conformers=molecule.conformers,
                                                    strict_n_conformers=False)

        # Try again, with strict_n_confs as true, but not including use_confs, so the
        # recommended number of confs will be generated
        molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry,
                                                strict_n_conformers=True)

        # Test calling the ToolkitWrapper _indirectly_, though the Molecule API,
        # which should raise the first error encountered
        with pytest.raises(ValueError,
                           match=f"has 2 conformers, but charge method 'am1bcc' "
                                 f"expects exactly 1."):
            molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry,
                                                   use_conformers=molecule.conformers,
                                                   strict_n_conformers=True)

        # Test calling the ToolkitWrapper _indirectly_, though a ToolkitRegistry,
        # specifying raise_exception_types=[]
        # which should aggregate any exceptions and bundle all of the messages
        # in a failed task together in a single ValueError.
        with pytest.raises(ValueError,
                           match=f"has 2 conformers, but charge method 'am1bcc' "
                                 f"expects exactly 1."):
            toolkit_registry.call('compute_partial_charges_am1bcc',
                                  molecule=molecule,
                                  use_conformers=molecule.conformers,
                                  strict_n_conformers=True,
                                  raise_exception_types=[])

        # Test calling the ToolkitWrapper _directly_, passing in the incorrect number of
        # confs, and specify strict_n_conformers, which should produce an IncorrectNumConformersError
        with pytest.raises(IncorrectNumConformersError,
                           match=f"has 2 conformers, but charge method 'am1bcc' "
                                 f"expects exactly 1."):
            ATTKW = AmberToolsToolkitWrapper()
            ATTKW.compute_partial_charges_am1bcc(molecule=molecule,
                                                 use_conformers=molecule.conformers,
                                                 strict_n_conformers=True)


    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    @pytest.mark.parametrize("partial_charge_method", ['am1bcc', 'am1-mulliken', 'gasteiger'])
    def test_assign_partial_charges_neutral(self, partial_charge_method):
        """Test AmberToolsToolkitWrapper assign_partial_charges()"""
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        molecule = create_ethanol()
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method)
        charge_sum = 0. * unit.elementary_charge
        for pc in molecule.partial_charges:
            charge_sum += pc
        assert -1.e-5 < charge_sum.value_in_unit(unit.elementary_charge) < 1.e-5

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    @pytest.mark.parametrize("partial_charge_method", ['am1bcc', 'am1-mulliken'])
    def test_assign_partial_charges_conformer_dependence(self, partial_charge_method):
        """Test AmberToolsToolkitWrapper assign_partial_charges()'s use_conformers kwarg
        to ensure charges are really conformer dependent. Skip Gasteiger because it isn't
        conformer dependent."""
        from openforcefield.tests.test_forcefield import create_ethanol
        import copy
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=1)
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method,
                                        use_conformers=molecule.conformers)
        pcs1 = copy.deepcopy(molecule.partial_charges)
        # This test case needs a pretty extreme coordinate change since ambertools only
        # stores partial charges to 1e-3
        molecule._conformers[0][0][0] += 3. * unit.angstrom
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method,
                                        use_conformers=molecule.conformers)
        for pc1, pc2 in zip(pcs1, molecule.partial_charges):
            assert abs(pc1 - pc2) > 1.e-3 * unit.elementary_charge

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    @pytest.mark.parametrize("partial_charge_method", ['am1bcc', 'am1-mulliken', 'gasteiger'])
    def test_assign_partial_charges_net_charge(self, partial_charge_method):
        """
        Test AmberToolsToolkitWrapper assign_partial_charges().
        """
        from openforcefield.tests.test_forcefield import create_acetate
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        molecule = create_acetate()
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method)
        charge_sum = 0. * unit.elementary_charge
        for pc in molecule.partial_charges:
            charge_sum += pc
        assert -1.01 < charge_sum.value_in_unit(unit.elementary_charge) < -0.99


    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_assign_partial_charges_bad_charge_method(self):
        """Test AmberToolsToolkitWrapper assign_partial_charges() for a nonexistent charge method"""
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        molecule = create_ethanol()

        # For now, ToolkitRegistries lose track of what exception type
        # was thrown inside them, so we just check for a ValueError here
        with pytest.raises(ValueError, match="is not available from AmberToolsToolkitWrapper") as excinfo:
            molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                            partial_charge_method="NotARealChargeMethod")

        # ToolkitWrappers raise a specific exception class, so we test that here
        with pytest.raises(ChargeMethodUnavailableError, match="is not available from AmberToolsToolkitWrapper") as excinfo:
            ATTKW = AmberToolsToolkitWrapper()
            ATTKW.assign_partial_charges(molecule=molecule,
                                         partial_charge_method="NotARealChargeMethod")


    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    @pytest.mark.parametrize("partial_charge_method,expected_n_confs", [('am1bcc', 1),
                                                                        ('am1-mulliken', 1),
                                                                        ('gasteiger', 0)])
    def test_assign_partial_charges_wrong_n_confs(self, partial_charge_method, expected_n_confs):
        """
        Test AmberToolsToolkitWrapper assign_partial_charges() when requesting to use an incorrect number of
        conformers
        """
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=2, rms_cutoff=0.01*unit.angstrom)

        # Try passing in the incorrect number of confs, but without specifying strict_n_conformers,
        # which should produce a warning
        with pytest.warns(IncorrectNumConformersWarning,
                          match=f"has 2 conformers, but charge method '{partial_charge_method}' "
                                f"expects exactly {expected_n_confs}."):
            molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                            partial_charge_method=partial_charge_method,
                                            use_conformers=molecule.conformers,
                                            strict_n_conformers=False)

        # Try again, with strict_n_confs as true, but not including use_confs, so the
        # recommended number of confs will be generated
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method,
                                        strict_n_conformers=True)

        # Test calling the ToolkitWrapper _indirectly_, though the Molecule API
        # which should aggregate any exceptions and bundle all of the messages
        # in a failed task together in a single ValueError.
        with pytest.raises(ValueError,
                           match=f"has 2 conformers, but charge method '{partial_charge_method}' "
                                 f"expects exactly {expected_n_confs}."):
            molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                            partial_charge_method=partial_charge_method,
                                            use_conformers=molecule.conformers,
                                            strict_n_conformers=True)

        # Test calling the ToolkitWrapper _directly_, passing in the incorrect number of
        # confs, and specify strict_n_conformers, which should produce an IncorrectNumConformersError
        with pytest.raises(IncorrectNumConformersError,
                           match=f"has 2 conformers, but charge method '{partial_charge_method}' "
                                 f"expects exactly {expected_n_confs}."):
            ATTKW = AmberToolsToolkitWrapper()
            ATTKW.assign_partial_charges(molecule=molecule,
                                         partial_charge_method=partial_charge_method,
                                         use_conformers=molecule.conformers,
                                         strict_n_conformers=True)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_assign_fractional_bond_orders(self):
        """Test OpenEyeToolkitWrapper assign_fractional_bond_orders()"""

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = toolkit_registry.call('from_smiles', smiles)
        for bond_order_model in ['am1-wiberg']:
            molecule.assign_fractional_bond_orders(toolkit_registry=toolkit_registry, bond_order_model=bond_order_model)
            # TODO: Add test for equivalent Wiberg orders for equivalent bonds

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_assign_fractional_bond_orders_neutral_charge_mol(self):
        """Test OpenEyeToolkitWrapper assign_fractional_bond_orders() for neutral and charged molecule.
        Also tests using existing conformers"""

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        # Reading neutral molecule from file
        filename = get_data_file_path('molecules/CID20742535_neutral.sdf')
        molecule1 = Molecule.from_file(filename)
        # Reading negative molecule from file
        filename = get_data_file_path('molecules/CID20742535_anion.sdf')
        molecule2 = Molecule.from_file(filename)

        # Checking that only one additional bond is present in the neutral molecule
        assert (len(molecule1.bonds) == len(molecule2.bonds) + 1)

        for bond_order_model in ['am1-wiberg']:
            molecule1.assign_fractional_bond_orders(toolkit_registry=toolkit_registry,
                                                    bond_order_model=bond_order_model,
                                                    use_conformers=molecule1.conformers)

            for i in molecule1.bonds:
                if i.is_aromatic:
                    # Checking aromatic bonds
                    assert (1.05 < i.fractional_bond_order < 1.65)
                elif (i.atom1.atomic_number == 1 or i.atom2.atomic_number == 1):
                    # Checking bond order of C-H or O-H bonds are around 1
                    assert (0.85 < i.fractional_bond_order < 1.05)
                elif (i.atom1.atomic_number == 8 or i.atom2.atomic_number == 8):
                    # Checking C-O single bond
                    wbo_C_O_neutral = i.fractional_bond_order
                    assert (1.0 < wbo_C_O_neutral < 1.5)
                else:
                    # Should be C-C single bond
                    assert (i.atom1_index == 4 and i.atom2_index == 6) or (i.atom1_index == 6 and i.atom2_index == 4)
                    wbo_C_C_neutral = i.fractional_bond_order
                    assert (1.0 < wbo_C_C_neutral < 1.3)

            molecule2.assign_fractional_bond_orders(toolkit_registry=toolkit_registry,
                                                    bond_order_model=bond_order_model,
                                                    use_conformers=molecule2.conformers)
            for i in molecule2.bonds:
                if i.is_aromatic:
                    # Checking aromatic bonds
                    assert (1.05 < i.fractional_bond_order < 1.65)

                elif (i.atom1.atomic_number == 1 or i.atom2.atomic_number == 1):
                    # Checking bond order of C-H or O-H bonds are around 1
                    assert (0.85 < i.fractional_bond_order < 1.05)
                elif (i.atom1.atomic_number == 8 or i.atom2.atomic_number == 8):
                    # Checking C-O single bond
                    wbo_C_O_anion = i.fractional_bond_order
                    assert (1.3 < wbo_C_O_anion < 1.8)
                else:
                    # Should be C-C single bond
                    assert (i.atom1_index == 4 and i.atom2_index == 6) or (i.atom1_index == 6 and i.atom2_index == 4)
                    wbo_C_C_anion = i.fractional_bond_order
                    assert (1.0 < wbo_C_C_anion < 1.3)

            # Wiberg bond order of C-C single bond is higher in the anion
            assert (wbo_C_C_anion > wbo_C_C_neutral)
            # Wiberg bond order of C-O bond is higher in the anion
            assert (wbo_C_O_anion > wbo_C_O_neutral)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_assign_fractional_bond_orders_charged(self):
        """Test OpenEyeToolkitWrapper assign_fractional_bond_orders() on a molecule with net charge +1"""

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = toolkit_registry.call('from_smiles', smiles)
        for bond_order_model in ['am1-wiberg']:
            molecule.assign_fractional_bond_orders(toolkit_registry=toolkit_registry,
                                                    bond_order_model=bond_order_model)
            # TODO: Add test for equivalent Wiberg orders for equivalent bonds

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_assign_fractional_bond_orders_invalid_method(self):
        """
        Test that AmberToolsToolkitWrapper.assign_fractional_bond_orders() raises the
        correct error if an invalid charge model is provided
        """

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = toolkit_registry.call('from_smiles', smiles)

        expected_error = "Bond order model 'not a real charge model' is not supported by " \
                         "AmberToolsToolkitWrapper. Supported models are ([[]'am1-wiberg'[]])"
        with pytest.raises(ValueError, match=expected_error) as excinfo:
            molecule.assign_fractional_bond_orders(toolkit_registry=AmberToolsToolkitWrapper(),
                                                   bond_order_model='not a real charge model')

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_assign_fractional_bond_orders_double_bond(self):
        """Test OpenEyeToolkitWrapper assign_fractional_bond_orders() on a molecule with a double bond"""

        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        smiles = r'C\C(F)=C(/F)C[C@@](C)(Cl)Br'
        molecule = toolkit_registry.call('from_smiles', smiles)
        for bond_order_model in ['am1-wiberg']:
            molecule.assign_fractional_bond_orders(toolkit_registry=toolkit_registry,
                                                   bond_order_model=bond_order_model)
            # TODO: Add test for equivalent Wiberg orders for equivalent bonds

        double_bond_has_wbo_near_2 = False
        for bond in molecule.bonds:
            if bond.bond_order == 2:
                if 1.75 < bond.fractional_bond_order < 2.25:
                    double_bond_has_wbo_near_2 = True
        assert double_bond_has_wbo_near_2


class TestBuiltInToolkitWrapper:
    """Test the BuiltInToolkitWrapper"""
    @pytest.mark.parametrize("partial_charge_method", ['zeros', 'formal_charge'])
    def test_assign_partial_charges_neutral(self, partial_charge_method):
        """Test BuiltInToolkitWrapper assign_partial_charges()"""
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[BuiltInToolkitWrapper])
        molecule = create_ethanol()
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method)
        charge_sum = 0. * unit.elementary_charge
        for pc in molecule.partial_charges:
            charge_sum += pc
        assert -1.e-6 < charge_sum.value_in_unit(unit.elementary_charge) < 1.e-6

    @pytest.mark.parametrize("partial_charge_method", ['formal_charge'])
    def test_assign_partial_charges_net_charge(self, partial_charge_method):
        """
        Test BuiltInToolkitWrapper assign_partial_charges(). Only formal_charge is tested, since zeros will not
        sum up to the proper number
        """
        from openforcefield.tests.test_forcefield import create_acetate
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[BuiltInToolkitWrapper])
        molecule = create_acetate()
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method=partial_charge_method)
        charge_sum = 0. * unit.elementary_charge
        for pc in molecule.partial_charges:
            charge_sum += pc
        assert -1.e-6 < charge_sum.value_in_unit(unit.elementary_charge) + 1. < 1.e-6


    def test_assign_partial_charges_bad_charge_method(self):
        """Test BuiltInToolkitWrapper assign_partial_charges() for a nonexistent charge method"""
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[BuiltInToolkitWrapper])
        molecule = create_ethanol()

        # For now, the Molecule API passes raise_exception_types=[] to ToolkitRegistry.call,
        # which loses track of what exception type
        # was thrown inside them, so we just check for a ValueError here
        with pytest.raises(ValueError, match="is not supported by the Built-in toolkit") as excinfo:
            molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                            partial_charge_method="NotARealChargeMethod")

        # ToolkitWrappers raise a specific exception class, so we test that here
        with pytest.raises(ChargeMethodUnavailableError, match="is not supported by the Built-in toolkit") as excinfo:
            BITKW = BuiltInToolkitWrapper()
            BITKW.assign_partial_charges(molecule=molecule,
                                         partial_charge_method="NotARealChargeMethod")

    def test_assign_partial_charges_wrong_n_confs(self):
        """
        Test BuiltInToolkitWrapper assign_partial_charges() when requesting to use an incorrect number of
        conformers
        """
        from openforcefield.tests.test_forcefield import create_ethanol
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[BuiltInToolkitWrapper])
        molecule = create_ethanol()
        molecule.generate_conformers(n_conformers=1)
        with pytest.warns(IncorrectNumConformersWarning,
                          match="has 1 conformers, but charge method 'zeros' expects exactly 0."):
            molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                            partial_charge_method="zeros",
                                            use_conformers=molecule.conformers,
                                            strict_n_conformers=False)

        # Specify strict_n_conformers=True, but not use_conformers, so a recommended number of
        # conformers will be generated internally
        molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                        partial_charge_method="zeros",
                                        strict_n_conformers=True)

        # For now, the Molecule API passes raise_exception_types=[] to ToolkitRegistry.call,
        # which loses track of what exception type
        # was thrown inside them, so we just check for a ValueError here
        with pytest.raises(ValueError,
                           match=f"has 1 conformers, but charge method 'zeros' "
                                 f"expects exactly 0."):
            molecule.assign_partial_charges(toolkit_registry=toolkit_registry,
                                            partial_charge_method="zeros",
                                            use_conformers=molecule.conformers,
                                            strict_n_conformers=True)

        # Test calling the ToolkitWrapper _directly_, passing in the incorrect number of
        # confs, and specify strict_n_conformers, which should produce an IncorrectNumConformersError
        with pytest.raises(IncorrectNumConformersError,
                           match=f"has 1 conformers, but charge method 'zeros' "
                                 f"expects exactly 0."):
            BITKW = BuiltInToolkitWrapper()
            BITKW.assign_partial_charges(molecule=molecule,
                                         partial_charge_method="zeros",
                                         use_conformers=molecule.conformers,
                                         strict_n_conformers=True)

class TestToolkitWrapper:
    """Test the ToolkitWrapper class"""
    def test_check_n_conformers(self):
        """Ensure that _check_n_conformers is working properly"""
        tkw = ToolkitWrapper()
        mol = create_ethanol()

        ## Test molecule with no conformers
        # Check with no min or max should pass
        tkw._check_n_conformers(mol, 'nocharge')
        # Check with min=1 should warn
        with pytest.warns(IncorrectNumConformersWarning,
                           match="has 0 conformers, but charge method 'nocharge' expects at least 1"):
            tkw._check_n_conformers(mol, 'nocharge', min_confs=1)
        # Check with min=1 and strict_n_conformers should raise an error
        with pytest.raises(IncorrectNumConformersError,
                           match="has 0 conformers, but charge method 'nocharge' expects at least 1"):
            tkw._check_n_conformers(mol, 'nocharge', min_confs=1, strict_n_conformers=True)
        # Check with min=1, max=1 and strict_n_conformers should raise an error
        with pytest.raises(IncorrectNumConformersError,
                           match="has 0 conformers, but charge method 'nocharge' expects exactly 1"):
            tkw._check_n_conformers(mol, 'nocharge', min_confs=1, max_confs=1, strict_n_conformers=True)
        # Check with min=1, max=2 and strict_n_conformers should raise an error
        with pytest.raises(IncorrectNumConformersError,
                           match="has 0 conformers, but charge method 'nocharge' expects between 1 and 2"):
            tkw._check_n_conformers(mol, 'nocharge', min_confs=1, max_confs=2, strict_n_conformers=True)
        # Check with max=1 should pass
        tkw._check_n_conformers(mol, 'nocharge', max_confs=1, strict_n_conformers=True)

        ## Test molecule with conformers
        # Add some conformers
        mol.generate_conformers(n_conformers=1)
        for _ in range(9):
            mol.add_conformer(mol.conformers[0])

        # Check with no min or max should pass
        tkw._check_n_conformers(mol, 'nocharge')

        ## min_confs checks
        # Check with min=1 should be fine
        tkw._check_n_conformers(mol, 'nocharge', min_confs=1)
        # Check with min=10 should be fine
        tkw._check_n_conformers(mol, 'nocharge', min_confs=10)
        # Check with min=11 should warn
        with pytest.warns(IncorrectNumConformersWarning,
                           match="has 10 conformers, but charge method 'nocharge' expects at least 11"):
            tkw._check_n_conformers(mol, 'nocharge', min_confs=11)
        # Check with min=11 and strict_n_conformers should raise an error
        with pytest.raises(IncorrectNumConformersError,
                           match="has 10 conformers, but charge method 'nocharge' expects at least 11"):
            tkw._check_n_conformers(mol, 'nocharge', min_confs=11, strict_n_conformers=True)

        ## max_confs checks
        # Check with max=1 and strict_n_conformers should raise an error
        with pytest.raises(IncorrectNumConformersError,
                           match="has 10 conformers, but charge method 'nocharge' expects at most 1"):
            tkw._check_n_conformers(mol, 'nocharge', max_confs=1, strict_n_conformers=True)
        # Check with max=10 and strict_n_conformers should be OK
        tkw._check_n_conformers(mol, 'nocharge', max_confs=10, strict_n_conformers=True)
        # Check with max=11 and strict_n_conformers should be OK
        tkw._check_n_conformers(mol, 'nocharge', max_confs=11, strict_n_conformers=True)

        ## min_confs and max_confs checks
        # Check with max=10 and min=10 and strict_n_conformers should be OK
        tkw._check_n_conformers(mol, 'nocharge', min_confs=10, max_confs=10, strict_n_conformers=True)
        # Check with max=10 and min=9 and strict_n_conformers should be OK
        tkw._check_n_conformers(mol, 'nocharge', min_confs=9, max_confs=10, strict_n_conformers=True)
        # Check with max=11 and min=10 and strict_n_conformers should be OK
        tkw._check_n_conformers(mol, 'nocharge', min_confs=10, max_confs=11, strict_n_conformers=True)
        # Check with max=11 and min=9 and strict_n_conformers should be OK
        tkw._check_n_conformers(mol, 'nocharge', min_confs=9, max_confs=11, strict_n_conformers=True)
        # Check with min=9 and max=9 and strict_n_conformers should raise an error
        with pytest.raises(IncorrectNumConformersError,
                           match="has 10 conformers, but charge method 'nocharge' expects exactly 9"):
            tkw._check_n_conformers(mol, 'nocharge', min_confs=9, max_confs=9, strict_n_conformers=True)
        # Check with min=1 and max=9 and strict_n_conformers should raise an error
        with pytest.raises(IncorrectNumConformersError,
                           match="has 10 conformers, but charge method 'nocharge' expects between 1 and 9"):
            tkw._check_n_conformers(mol, 'nocharge', min_confs=1, max_confs=9, strict_n_conformers=True)
        # Check with min=11 and max=12 and strict_n_conformers should raise an error
        with pytest.raises(IncorrectNumConformersError,
                           match="has 10 conformers, but charge method 'nocharge' expects between 11 and 12"):
            tkw._check_n_conformers(mol, 'nocharge', min_confs=11, max_confs=12, strict_n_conformers=True)


class TestToolkitRegistry:
    """Test the ToolkitRegistry class"""

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_add_bad_toolkit(self):
        registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])
        with pytest.raises(InvalidToolkitError):
            registry.add_toolkit('rdkit as a string')

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    @pytest.mark.skipif(OpenEyeToolkitWrapper.is_available(), reason='Skipping while OpenEye is available')
    def test_register_unavailable_toolkit(self):
        registry = ToolkitRegistry(toolkit_precedence=[RDKitToolkitWrapper])
        with pytest.raises(ToolkitUnavailableException):
            registry.register_toolkit(toolkit_wrapper=OpenEyeToolkitWrapper, exception_if_unavailable=True)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_register_openeye(self):
        """Test creation of toolkit registry with OpenEye toolkit"""
        # Test registration of OpenEyeToolkitWrapper
        toolkit_precedence = [OpenEyeToolkitWrapper]
        registry = ToolkitRegistry(toolkit_precedence=toolkit_precedence, register_imported_toolkit_wrappers=False)

        assert {type(c) for c in registry.registered_toolkits} == {OpenEyeToolkitWrapper}

        # Test ToolkitRegistry.resolve()
        assert registry.resolve('to_smiles') == registry.registered_toolkits[0].to_smiles

        # Test ToolkitRegistry.call()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = registry.call('from_smiles', smiles)
        smiles2 = registry.call('to_smiles', molecule)
        assert smiles == smiles2

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_register_rdkit(self):
        """Test creation of toolkit registry with RDKit toolkit"""
        # Test registration of RDKitToolkitWrapper
        toolkit_precedence = [RDKitToolkitWrapper]
        registry = ToolkitRegistry(toolkit_precedence=toolkit_precedence,
                                   register_imported_toolkit_wrappers=False)

        assert set([ type(c) for c in registry.registered_toolkits]) == {RDKitToolkitWrapper}

        # Test ToolkitRegistry.resolve()
        assert registry.resolve('to_smiles') == registry.registered_toolkits[0].to_smiles

        # Test ToolkitRegistry.call()
        smiles = '[H][C]([H])([H])[C]([H])([H])[H]'
        molecule = registry.call('from_smiles', smiles)
        smiles2 = registry.call('to_smiles', molecule)
        assert smiles == smiles2

    @pytest.mark.skipif(
        not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_register_ambertools(self):
        """Test creation of toolkit registry with AmberToolsToolkitWrapper and RDKitToolkitWrapper
        """
        # Test registration of AmberToolsToolkitWrapper
        toolkit_precedence = [AmberToolsToolkitWrapper, RDKitToolkitWrapper]
        registry = ToolkitRegistry(toolkit_precedence=toolkit_precedence,
                                   register_imported_toolkit_wrappers=False)
        #registry.register_toolkit(AmberToolsToolkitWrapper)
        assert set([ type(c) for c in registry.registered_toolkits]) == set([AmberToolsToolkitWrapper,RDKitToolkitWrapper])

        # Test ToolkitRegistry.resolve()
        registry.resolve('assign_partial_charges')
        assert registry.resolve('assign_partial_charges') == registry.registered_toolkits[0].assign_partial_charges

        # Test ToolkitRegistry.call()
        registry.register_toolkit(RDKitToolkitWrapper)
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = registry.call('from_smiles', smiles)
        smiles2 = registry.call('to_smiles', molecule)
        assert smiles == smiles2

    @pytest.mark.skipif(
        not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
        reason='RDKitToolkit or AmberToolsToolkit is not available')
    def test_deregister_toolkit(self):
        """Test removing an instantiated toolkit from the registry"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])

        assert any([isinstance(tk, AmberToolsToolkitWrapper) for tk in toolkit_registry._toolkits])
        assert any([isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits])

        toolkit_registry.deregister_toolkit(toolkit_registry._toolkits[-1])
        assert any([isinstance(tk, AmberToolsToolkitWrapper) for tk in toolkit_registry._toolkits])
        assert not any([isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits])

        toolkit_registry.deregister_toolkit(toolkit_registry._toolkits[-1])
        assert not any([isinstance(tk, AmberToolsToolkitWrapper) for tk in toolkit_registry._toolkits])
        assert not any([isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits])

    @pytest.mark.skipif(
        not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_deregister_toolkit_by_class(self):
        """Test removing a toolkit from the registry by matching class types"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])

        assert any([isinstance(tk, AmberToolsToolkitWrapper) for tk in toolkit_registry._toolkits])
        assert any([isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits])

        toolkit_registry.deregister_toolkit(RDKitToolkitWrapper)
        assert any([isinstance(tk, AmberToolsToolkitWrapper) for tk in toolkit_registry._toolkits])
        assert not any([isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits])

        toolkit_registry.deregister_toolkit(AmberToolsToolkitWrapper)
        assert not any([isinstance(tk, AmberToolsToolkitWrapper) for tk in toolkit_registry._toolkits])
        assert not any([isinstance(tk, RDKitToolkitWrapper) for tk in toolkit_registry._toolkits])

    @pytest.mark.skipif(
        not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_deregister_toolkit_bad_inputs(self):
        """Test bad inputs to deregister_toolkit"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])

        with pytest.raises(InvalidToolkitError):
            toolkit_registry.deregister_toolkit('rdkit as a string')

    def test_register_builtintoolkit(self):
        """Test creation of toolkit registry with Built-in toolkit"""
        # Test registration of BuiltInToolkitWrapper
        toolkit_precedence = [BuiltInToolkitWrapper]
        registry = ToolkitRegistry(toolkit_precedence=toolkit_precedence,
                                   register_imported_toolkit_wrappers=False)
        #registry.register_toolkit(BuiltInToolkitWrapper)
        assert set([ type(c) for c in registry.registered_toolkits]) == set([BuiltInToolkitWrapper])

        # Test ToolkitRegistry.resolve()
        assert registry.resolve('assign_partial_charges') == registry.registered_toolkits[0].assign_partial_charges

    @pytest.mark.skipif(not AmberToolsToolkitWrapper.is_available(), reason='AmberTools Toolkit not available')
    def test_call_raise_first_error(self):
        """Test to ensure proper behavior of raise_first_error kwarg to ToolkitRegistry.call"""
        toolkit_precedence = [BuiltInToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper]
        registry = ToolkitRegistry(toolkit_precedence=toolkit_precedence,
                                   register_imported_toolkit_wrappers=False)
        mol = registry.call('from_smiles', 'C')
        # Specify that the ToolkitRegistry should raise the first ChargeMethodUnavailableError it encounters
        with pytest.raises(ChargeMethodUnavailableError, match='"notarealchargemethod"" is not supported by the Built-in toolkit.'):
            registry.call('assign_partial_charges',
                          molecule=mol,
                          partial_charge_method="NotARealChargeMethod",
                          raise_exception_types=[ChargeMethodUnavailableError])
        # Specify that the ToolkitRegistry should collect all the errors it encounters and
        # ensure it raises a single ValueError when no ToolkitWrappers succeed
        with pytest.raises(ValueError, match="partial_charge_method \'notarealchargemethod\' is not available from AmberToolsToolkitWrapper"):
            registry.call('assign_partial_charges',
                          molecule=mol,
                          partial_charge_method="NotARealChargeMethod",
                          raise_exception_types=[])
