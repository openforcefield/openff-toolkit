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

import pytest
from openforcefield.utils.toolkits import (OpenEyeToolkitWrapper, RDKitToolkitWrapper,
                                           AmberToolsToolkitWrapper, ToolkitRegistry,
                                           GAFFAtomTypeWarning, UndefinedStereochemistryError)
from openforcefield.utils import get_data_file_path
from openforcefield.topology.molecule import Molecule


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
                    molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)
            else:
                molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_wrapper)

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
        assert (molecule._conformers[0] == molecule2._conformers[0]).all()
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
        rdmol = molecule.to_openeye()
        molecule2 = Molecule.from_openeye(rdmol)

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
        assert (molecule._conformers == None)
        assert (molecule2._conformers == None)
        for pc1, pc2 in zip(molecule._partial_charges, molecule2._partial_charges):
            pc1_ul = pc1 / unit.elementary_charge
            pc2_ul = pc2 / unit.elementary_charge
            assert_almost_equal(pc1_ul, pc2_ul, decimal=6)
        assert molecule2.to_smiles(toolkit_registry=toolkit_wrapper) == expected_output_smiles

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_get_sdf_coordinates(self):
        """Test OpenEyeToolkitWrapper for importing a single set of coordinates from a sdf file"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    @pytest.mark.skip
    def test_get_multiconformer_sdf_coordinates(self):
        """Test OpenEyeToolkitWrapper for importing multiple sets of coordinates from a sdf file"""
        raise NotImplementedError
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_get_mol2_coordinates(self):
        """Test OpenEyeToolkitWrapper for importing a single set of molecule coordinates"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.mol2')
        molecule1 = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule1._conformers) == 1
        assert molecule1._conformers[0].shape == (15, 3)
        assert_almost_equal(molecule1.conformers[0][5][1] / unit.angstrom, 22.98, decimal=2)

        # Test loading from file-like object
        with open(filename, 'r') as infile:
            molecule2 = Molecule(infile, file_format='MOL2', toolkit_registry=toolkit_wrapper)
        assert molecule1.is_isomorphic(molecule2)
        assert len(molecule2._conformers) == 1
        assert molecule2._conformers[0].shape == (15, 3)
        assert_almost_equal(molecule2.conformers[0][5][1] / unit.angstrom, 22.98, decimal=2)

        # Test loading from gzipped mol2
        import gzip
        with gzip.GzipFile(filename + '.gz', 'r') as infile:
            molecule3 = Molecule(infile, file_format='MOL2', toolkit_registry=toolkit_wrapper)
        assert molecule1.is_isomorphic(molecule3)
        assert len(molecule3._conformers) == 1
        assert molecule3._conformers[0].shape == (15, 3)
        assert_almost_equal(molecule3.conformers[0][5][1] / unit.angstrom, 22.98, decimal=2)

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_get_mol2_charges(self):
        """Test OpenEyeToolkitWrapper for importing a mol2 file specifying partial charges"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_file_path('molecules/toluene_charged.mol2')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)
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

        # TODO: Make this test more robust

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_compute_partial_charges(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        # Ensure that an exception is raised if no conformers are provided
        with pytest.raises(Exception) as excinfo:
            molecule.compute_partial_charges(toolkit_registry=toolkit_wrapper)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        # Ensure that an exception is raised if an invalid charge model is passed in
        with pytest.raises(Exception) as excinfo:
            charge_model = 'notARealChargeModel'
            molecule.compute_partial_charges(toolkit_registry=toolkit_wrapper, charge_model=charge_model)

        # TODO: Test all supported charge models
        # Note: "amber" and "amberff94" only work for a subset of residue types, so we'll need to find testing data for
        # those
        # charge_model = [,'amber','amberff94']
        # TODO: 'mmff' and 'mmff94' often assign charges of 0, presumably if the molecule is unrecognized.
        # charge_model = ['mmff', 'mmff94']
        for charge_model in ['noop', 'am1bcc', 'am1bccnosymspt', 'am1bccelf10']:
            with pytest.raises(NotImplementedError) as excinfo:
                molecule.compute_partial_charges(toolkit_registry=toolkit_wrapper)  # , charge_model=charge_model)
                charge_sum = 0 * unit.elementary_charge
                for pc in molecule._partial_charges:
                    charge_sum += pc
                assert charge_sum < 0.001 * unit.elementary_charge

        # For now, just test AM1-BCC while the SMIRNOFF spec for other charge models gets worked out
        molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_wrapper)  # , charge_model=charge_model)
        charge_sum = 0 * unit.elementary_charge
        for pc in molecule._partial_charges:
            charge_sum += pc
        assert charge_sum < 0.001 * unit.elementary_charge

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_compute_partial_charges_net_charge(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges() on a molecule with a net +1 charge"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)


        with pytest.raises(NotImplementedError) as excinfo:
            charge_model = 'notARealChargeModel'
            molecule.compute_partial_charges(toolkit_registry=toolkit_wrapper)#, charge_model=charge_model)

        # TODO: Test all supported charge models
        # TODO: "amber" and "amberff94" only work for a subset of residue types, so we'll need to find testing data for
        # those
        # charge_model = [,'amber','amberff94']
        # The 'noop' charge model doesn't add up to the formal charge, so we shouldn't test it
        # charge_model = ['noop']
        for charge_model in ['mmff', 'mmff94', 'am1bcc', 'am1bccnosymspt', 'am1bccelf10']:
            with pytest.raises(NotImplementedError) as excinfo:
                molecule.compute_partial_charges(toolkit_registry=toolkit_wrapper) #, charge_model=charge_model)
                charge_sum = 0 * unit.elementary_charge
                for pc in molecule._partial_charges:
                    charge_sum += pc
                assert 0.999 * unit.elementary_charge < charge_sum < 1.001 * unit.elementary_charge
        # For now, I'm just testing AM1-BCC (will test more when the SMIRNOFF spec for other charges is finalized)
        molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_wrapper)
        charge_sum = 0 * unit.elementary_charge
        for pc in molecule._partial_charges:
            charge_sum += pc
        assert 0.999 * unit.elementary_charge < charge_sum < 1.001 * unit.elementary_charge

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_compute_wiberg_bond_orders(self):
        """Test OpenEyeToolkitWrapper compute_wiberg_bond_orders()"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        for charge_model in ['am1','pm3']:
            molecule.compute_wiberg_bond_orders(toolkit_registry=toolkit_wrapper, charge_model=charge_model)
            print([bond.fractional_bond_order for bond in molecule.bonds])
            # TODO: Add test for equivalent Wiberg orders for equivalent bonds

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_compute_wiberg_bond_orders_charged(self):
        """Test OpenEyeToolkitWrapper compute_wiberg_bond_orders() on a molecule with net charge +1"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        for charge_model in ['am1', 'pm3']:
            molecule.compute_wiberg_bond_orders(toolkit_registry=toolkit_wrapper, charge_model=charge_model)
            # TODO: Add test for equivalent Wiberg orders for equivalent bonds

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(),
    reason='OpenEye Toolkit not available')
    def test_compute_wiberg_bond_orders_double_bond(self):
        """Test OpenEyeToolkitWrapper compute_wiberg_bond_orders() on a molecule with a double bond"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = r'C\C(F)=C(/F)C[C@@](C)(Cl)Br'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        for charge_model in ['am1', 'pm3']:
            molecule.compute_wiberg_bond_orders(toolkit_registry=toolkit_wrapper, charge_model=charge_model)
            # TODO: Add test for equivalent Wiberg orders for equivalent bonds

        double_bond_has_wbo_near_2 = False
        for bond in molecule.bonds:
            if bond.bond_order == 2:
                if 1.75 < bond.fractional_bond_order < 2.25:
                    double_bond_has_wbo_near_2 = True
        assert double_bond_has_wbo_near_2





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
        assert (molecule._conformers[0] == molecule2._conformers[0]).all()
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

        # Do a first conversion to/from oemol
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
        assert (molecule._conformers == None)
        assert (molecule2._conformers == None)
        for pc1, pc2 in zip(molecule._partial_charges, molecule2._partial_charges):
            pc1_ul = pc1 / unit.elementary_charge
            pc2_ul = pc2 / unit.elementary_charge
            assert_almost_equal(pc1_ul, pc2_ul, decimal=6)
        assert molecule2.to_smiles(toolkit_registry=toolkit_wrapper) == expected_output_smiles
        
    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_get_sdf_coordinates(self):
        """Test RDKitToolkitWrapper for importing a single set of coordinates from a sdf file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15, 3)
        assert_almost_equal(molecule.conformers[0][5][1] / unit.angstrom, 2.0104, decimal=4)

    # Find a multiconformer SDF file
    @pytest.mark.skip
    #@pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_get_multiconformer_sdf_coordinates(self):
        """Test RDKitToolkitWrapper for importing a single set of coordinates from a sdf file"""
        raise NotImplementedError
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)

    # Unskip this when we implement PDB-reading support for RDKitToolkitWrapper
    @pytest.mark.skip
    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_get_pdb_coordinates(self):
        """Test RDKitToolkitWrapper for importing a single set of coordinates from a pdb file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.pdb')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)

    # Unskip this when we implement PDB-reading support for RDKitToolkitWrapper
    @pytest.mark.skip
    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_load_aromatic_pdb(self):
        """Test OpenEyeToolkitWrapper for importing molecule conformers"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_file_path('molecules/toluene.pdb')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)

    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available(), reason='RDKit Toolkit not available')
    def test_generate_conformers(self):
        """Test RDKitToolkitWrapper generate_conformers()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers()
        # TODO: Make this test more robust
        




        
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
    def test_compute_partial_charges(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges()"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])

        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_registry)
        molecule.generate_conformers(toolkit_registry=toolkit_registry)

        # TODO: Implementation of these tests is pending a decision on the API for our charge model
        with pytest.raises(NotImplementedError) as excinfo:
            charge_model = 'notARealChargeModel'
            molecule.compute_partial_charges(toolkit_registry=toolkit_registry)#, charge_model=charge_model)

        # ['cm1', 'cm2']
        for charge_model in ['gas', 'mul', 'bcc']:
            with pytest.raises(NotImplementedError) as excinfo:
                molecule.compute_partial_charges(toolkit_registry=toolkit_registry)#, charge_model=charge_model)
                charge_sum = 0 * unit.elementary_charge
                for pc in molecule._partial_charges:
                    charge_sum += pc
                assert charge_sum < 0.01 * unit.elementary_charge

        # For now, just test AM1-BCC while the SMIRNOFF spec for other charge models gets worked out
        molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry)  # , charge_model=charge_model)
        charge_sum = 0 * unit.elementary_charge
        for pc in molecule._partial_charges:
            charge_sum += pc
        assert charge_sum < 0.002 * unit.elementary_charge



    @pytest.mark.skipif(not RDKitToolkitWrapper.is_available() or not AmberToolsToolkitWrapper.is_available(),
                        reason='RDKitToolkit and AmberToolsToolkit not available')
    def test_compute_partial_charges_net_charge(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges() on a molecule with a net +1 charge"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_registry)
        molecule.generate_conformers(toolkit_registry=toolkit_registry)


        with pytest.raises(NotImplementedError) as excinfo:
            charge_model = 'notARealChargeModel'
            molecule.compute_partial_charges(toolkit_registry=toolkit_registry)#, charge_model=charge_model)

        # TODO: Figure out why ['cm1', 'cm2'] fail
        for charge_model in  ['gas', 'mul', 'bcc']:
            with pytest.raises(NotImplementedError) as excinfo:
                molecule.compute_partial_charges(toolkit_registry=toolkit_registry)#, charge_model=charge_model)
                charge_sum = 0 * unit.elementary_charge
                for pc in molecule._partial_charges:
                    charge_sum += pc
                assert 0.99 * unit.elementary_charge < charge_sum < 1.01 * unit.elementary_charge

        # For now, I'm just testing AM1-BCC (will test more when the SMIRNOFF spec for other charges is finalized)
        molecule.compute_partial_charges_am1bcc(toolkit_registry=toolkit_registry)
        charge_sum = 0 * unit.elementary_charge
        for pc in molecule._partial_charges:
            charge_sum += pc
        assert 0.999 * unit.elementary_charge < charge_sum < 1.001 * unit.elementary_charge


class TestToolkitRegistry:
    """Test the ToolkitRegistry"""

    @pytest.mark.skipif(not OpenEyeToolkitWrapper.is_available(), reason='OpenEye Toolkit not available')
    def test_register_openeye(self):
        """Test creation of toolkit registry with OpenEye toolkit"""
        # Test registration of OpenEyeToolkitWrapper
        toolkit_precedence = [OpenEyeToolkitWrapper]
        registry = ToolkitRegistry(toolkit_precedence=toolkit_precedence, register_imported_toolkit_wrappers=False)
        #registry.register_toolkit(OpenEyeToolkitWrapper)
        assert set([type(c) for c in registry.registered_toolkits]) == set([OpenEyeToolkitWrapper])

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
        #registry.register_toolkit(RDKitToolkitWrapper)
        assert set([ type(c) for c in registry.registered_toolkits]) == set([RDKitToolkitWrapper])

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
        registry.resolve('compute_partial_charges')
        assert registry.resolve('compute_partial_charges') == registry.registered_toolkits[0].compute_partial_charges

        # Test ToolkitRegistry.call()
        registry.register_toolkit(RDKitToolkitWrapper)
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = registry.call('from_smiles', smiles)
        #partial_charges = registry.call('compute_partial_charges', molecule)

