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

from functools import partial
from unittest import TestCase
from openforcefield import utils
from simtk import unit
import numpy as np

import pytest
from openforcefield.utils.toolkits import ToolkitWrapper, OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, ToolkitRegistry

from openforcefield.utils import get_data_filename

from openforcefield.topology.molecule import Molecule

#=============================================================================================
# TESTS
#=============================================================================================

class TestToolkitWrapper(TestCase):
    """Test the ToolkitWrapper abstract base class"""
    pass


class TestOpenEyeToolkitWrapper(TestCase):
    """Test the OpenEyeToolkitWrapper"""

    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_smiles(self):
        """Test OpenEyeToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        # This differs from RDKit's SMILES due to different canonicalization schemes
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles,
                                        toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2

    @OpenEyeToolkitWrapper.requires_toolkit()
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

    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_smiles_charged(self):
        """Test OpenEyeToolkitWrapper functions for reading/writing charged SMILES"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        # This differs from RDKit's expected output due to different canonicalization schemes
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles,
                                        toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2


    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_to_from_openeye(self):
        """Test OpenEyeToolkitWrapper to_openeye() and from_openeye()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        # This differs from RDKit's expected output due to different canonicalization schemes
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles)
        oemol = molecule.to_openeye()
        molecule2 = Molecule.from_openeye(oemol)
        new_conf1 = unit.Quantity(np.array([[1,2,3],[4,5,6],[7,8,9],
                                           [10,11,12],[13,14,15],[16,17,18],
                                           [19,20,21],[22,23,24],[25,26,27]],
                                          dtype=np.float),
                                 unit.angstrom)
        new_conf2 = unit.Quantity(np.array([[101,102,103],[104,105,106],[107,108,109],
                                            [110,111,112],[113,114,115],[116,117,118],
                                            [119,120,121],[122,123,124],[125,126,127]],
                                          dtype=np.float),
                                 unit.angstrom)
        molecule2.add_conformer(new_conf1)
        molecule2.add_conformer(new_conf2)
        
        
        smiles2 = molecule2.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2
        oemol2 = Molecule.to_openeye(molecule2)
        assert oemol2.NumConfs() == 2
        molecule3 = Molecule.from_openeye(oemol2)
        assert len(molecule3._conformers) == 2
        assert (molecule2._conformers[0] == molecule3._conformers[0]).all()
        assert (molecule2._conformers[1] == molecule3._conformers[1]).all()
        
    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_get_sdf_coordinates(self):
        """Test OpenEyeToolkitWrapper for importing a single set of coordinates from a sdf file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_filename('molecules/toluene.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)

    
    @pytest.mark.skip #Implement this
    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_get_multiconformer_sdf_coordinates(self):
        """Test OpenEyeToolkitWrapper for importing a single set of coordinates from a sdf file"""
        raise NotImplementedError
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_filename('molecules/toluene.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)
        

    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_get_mol2_coordinates(self):
        """Test OpenEyeToolkitWrapper for importing a single set of molecule coordinates"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        filename = get_data_filename('molecules/toluene.mol2')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)


    #@pytest.mark.skipif( not OpenEyeToolkitWrapper.toolkit_is_available() )
    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_generate_conformers(self):
        """Test OpenEyeToolkitWrapper generate_conformers()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers()
        # TODO: Make this test more robust

    # @pytest.mark.skipif( not OpenEyeToolkitWrapper.toolkit_is_available() )
    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_compute_partial_charges(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        # Ensure that an exception is raised if no conformers are provided
        with self.assertRaises(Exception) as context:
            molecule.compute_partial_charges(toolkit_registry=toolkit_wrapper)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        # Ensure that an exception is raised if an invalid charge model is passed in
        with self.assertRaises(Exception) as context:
            charge_model = 'notARealChargeModel'
            molecule.compute_partial_charges(toolkit_registry=toolkit_wrapper, charge_model=charge_model)

        # TODO: Test all supported charge models
        # Note: "amber" and "amberff94" only work for a subset of residue types, so we'll need to find testing data for
        # those
        # charge_model = [,'amber','amberff94']
        # TODO: 'mmff' and 'mmff94' often assign charges of 0, presumably if the molecule is unrecognized.
        # charge_model = ['mmff', 'mmff94']
        for charge_model in ['noop', 'am1bcc', 'am1bccnosymspt', 'am1bccelf10']:
            molecule.compute_partial_charges(toolkit_registry=toolkit_wrapper, charge_model=charge_model)
            charge_sum = 0 * unit.elementary_charge
            for pc in molecule._partial_charges:
                charge_sum += pc
            assert charge_sum < 0.001 * unit.elementary_charge

    # @pytest.mark.skipif( not OpenEyeToolkitWrapper.toolkit_is_available() )
    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_compute_partial_charges_net_charge(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges() on a molecule with a net +1 charge"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)

        with self.assertRaises(Exception) as context:
            charge_model = 'notARealChargeModel'
            molecule.compute_partial_charges(toolkit_registry=toolkit_wrapper, charge_model=charge_model)

        # TODO: Test all supported charge models
        # TODO: "amber" and "amberff94" only work for a subset of residue types, so we'll need to find testing data for
        # those
        # charge_model = [,'amber','amberff94']
        # The 'noop' charge model doesn't add up to the formal charge, so we shouldn't test it
        # charge_model = ['noop']
        for charge_model in ['mmff', 'mmff94', 'am1bcc', 'am1bccnosymspt', 'am1bccelf10']:
            molecule.compute_partial_charges(toolkit_registry=toolkit_wrapper, charge_model=charge_model)
            charge_sum = 0 * unit.elementary_charge
            for pc in molecule._partial_charges:
                charge_sum += pc
            assert 0.999 * unit.elementary_charge < charge_sum < 1.001 * unit.elementary_charge

    @OpenEyeToolkitWrapper.requires_toolkit()
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


    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_compute_wiberg_bond_orders_charged(self):
        """Test OpenEyeToolkitWrapper compute_wiberg_bond_orders() on a molecule with net charge +1"""

        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        molecule.generate_conformers(toolkit_registry=toolkit_wrapper)
        for charge_model in ['am1','pm3']:
            molecule.compute_wiberg_bond_orders(toolkit_registry=toolkit_wrapper, charge_model=charge_model)
            # TODO: Add test for equivalent Wiberg orders for equivalent bonds





        # TODO: Check partial charge invariants (total charge, charge equivalence)
        # TODO: Add test for partial charge calculation with formally charged atoms/molecules
        
        # TODO: Add test for higher bonds orders
        # TODO: Add test for aromaticity
        # TODO: Add test and molecule functionality for isotopes


        
class TestRDKitToolkitWrapper(TestCase):
    """Test the RDKitToolkitWrapper"""
    
    #@pytest.mark.skipif( not RDKitToolkitWrapper.toolkit_is_available() )
    @RDKitToolkitWrapper.requires_toolkit()
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

    @RDKitToolkitWrapper.requires_toolkit()
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


    @RDKitToolkitWrapper.requires_toolkit()
    def test_smiles_charged(self):
        """Test RDKitWrapper functions for reading/writing charged SMILES"""
        toolkit_wrapper = RDKitToolkitWrapper()
        # This differs from OE's expected output due to different canonicalization schemes
        smiles = '[H][C]([H])([H])[N+]([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles,
                                        toolkit_registry=toolkit_wrapper)
        smiles2 = molecule.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2

    #@pytest.mark.skipif( not RDKitToolkitWrapper.toolkit_is_available() )
    @RDKitToolkitWrapper.requires_toolkit()
    def test_to_from_rdkit(self):
        """Test RDKitToolkitWrapper to_rdkit() and from_rdkit()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        # This differs from OE's expected output due to different canonicalization schemes
        smiles = '[H][C]([H])([H])[C]([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles,
                                        toolkit_registry=toolkit_wrapper)
        rdmol = molecule.to_rdkit()
        molecule2 = Molecule.from_rdkit(rdmol)
        smiles2 = molecule2.to_smiles(toolkit_registry=toolkit_wrapper)
        new_conf1 = unit.Quantity(np.array([[1,2,3],[4,5,6],[7,8,9],
                                           [10,11,12],[13,14,15],[16,17,18],
                                           [19,20,21],[22,23,24],[25,26,27]],
                                          dtype=np.float),
                                 unit.angstrom)
        new_conf2 = unit.Quantity(np.array([[101,102,103],[104,105,106],[107,108,109],
                                            [110,111,112],[113,114,115],[116,117,118],
                                            [119,120,121],[122,123,124],[125,126,127]],
                                          dtype=np.float),
                                 unit.angstrom)
        molecule2.add_conformer(new_conf1)
        molecule2.add_conformer(new_conf2)
        
        
        smiles2 = molecule2.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2
        oemol2 = Molecule.to_openeye(molecule2)
        assert oemol2.NumConfs() == 2
        molecule3 = Molecule.from_openeye(oemol2)
        assert len(molecule3._conformers) == 2
        assert (molecule2._conformers[0] == molecule3._conformers[0]).all()
        assert (molecule2._conformers[1] == molecule3._conformers[1]).all()
        
        
    @RDKitToolkitWrapper.requires_toolkit()
    def test_get_sdf_coordinates(self):
        """Test RDKitToolkitWrapper for importing a single set of coordinates from a sdf file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_filename('molecules/toluene.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)
    
    @pytest.mark.skip #Implement this
    @RDKitToolkitWrapper.requires_toolkit()
    def test_get_multiconformer_sdf_coordinates(self):
        """Test RDKitToolkitWrapper for importing a single set of coordinates from a sdf file"""
        raise NotImplementedError
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_filename('molecules/toluene.sdf')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)
    
        
    # Unskip this when we implement PDB-reading support for RDKitToolkitWrapper
    @pytest.mark.skip
    @RDKitToolkitWrapper.requires_toolkit()
    def test_get_pdb_coordinates(self):
        """Test RDKitToolkitWrapper for importing a single set of coordinates from a pdb file"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_filename('molecules/toluene.pdb')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)


    # Unskip this when we implement PDB-reading support for RDKitToolkitWrapper
    @pytest.mark.skip
    @RDKitToolkitWrapper.requires_toolkit()
    def test_load_aromatic_pdb(self):
        """Test OpenEyeToolkitWrapper for importing molecule conformers"""
        toolkit_wrapper = RDKitToolkitWrapper()
        filename = get_data_filename('molecules/toluene.pdb')
        molecule = Molecule.from_file(filename, toolkit_registry=toolkit_wrapper)
        assert len(molecule._conformers) == 1
        assert molecule._conformers[0].shape == (15,3)



    #@pytest.mark.skipif( not OpenEyeToolkitWrapper.toolkit_is_available() )
    @RDKitToolkitWrapper.requires_toolkit()
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



        
class TestAmberToolsWrapper(TestCase):
    """Test the AmberToolsWraper"""



    # @pytest.mark.skipif( not OpenEyeToolkitWrapper.toolkit_is_available() )
    @RDKitToolkitWrapper.requires_toolkit()
    @AmberToolsToolkitWrapper.requires_toolkit()
    def test_compute_partial_charges(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges()"""
        toolkit_registry = ToolkitRegistry(toolkit_precedence=[AmberToolsToolkitWrapper, RDKitToolkitWrapper])

        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_registry)
        molecule.generate_conformers(toolkit_registry=toolkit_registry)

        with self.assertRaises(Exception) as context:
            charge_model = 'notARealChargeModel'
            molecule.compute_partial_charges(toolkit_registry=toolkit_registry, charge_model=charge_model)

        # ['cm1', 'cm2']
        for charge_model in ['gas', 'mul', 'bcc']:
            molecule.compute_partial_charges(toolkit_registry=toolkit_registry, charge_model=charge_model)
            charge_sum = 0 * unit.elementary_charge
            for pc in molecule._partial_charges:
                charge_sum += pc
            assert charge_sum < 0.01 * unit.elementary_charge


    # @pytest.mark.skipif( not OpenEyeToolkitWrapper.toolkit_is_available() )
    @RDKitToolkitWrapper.requires_toolkit()
    @AmberToolsToolkitWrapper.requires_toolkit()
    def test_compute_partial_charges_net_charge(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges() on a molecule with a net +1 charge"""
        toolkit_registry = ToolkitRegistry([AmberToolsToolkitWrapper, RDKitToolkitWrapper])
        smiles = '[H]C([H])([H])[N+]([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles, toolkit_registry=toolkit_registry)
        molecule.generate_conformers(toolkit_registry=toolkit_registry)

        with self.assertRaises(Exception) as context:
            charge_model = 'notARealChargeModel'
            molecule.compute_partial_charges(toolkit_registry=toolkit_registry, charge_model=charge_model)

        # TODO: Figure out why ['cm1', 'cm2'] fail
        for charge_model in  ['gas', 'mul', 'bcc']:
            molecule.compute_partial_charges(toolkit_registry=toolkit_registry, charge_model=charge_model)
            charge_sum = 0 * unit.elementary_charge
            for pc in molecule._partial_charges:
                charge_sum += pc
            assert 0.99 * unit.elementary_charge < charge_sum < 1.01 * unit.elementary_charge



class TestToolkitRegistry(TestCase):
    """Test the ToolkitRegistry"""

    @OpenEyeToolkitWrapper.requires_toolkit()
    #@pytest.mark.skipif( not OpenEyeToolkitWrapper.toolkit_is_available() )
    def test_register_openeye(self):
        """Test creation of toolkit registry with OpenEye toolkit"""
        # Test registration of OpenEyeToolkitWrapper
        toolkit_precedence = [OpenEyeToolkitWrapper]
        registry = ToolkitRegistry(toolkit_precedence=toolkit_precedence, register_imported_toolkit_wrappers=False)
        #registry.register_toolkit(OpenEyeToolkitWrapper)
        assert set([ type(c) for c in registry.registered_toolkits]) == set([OpenEyeToolkitWrapper])

        # Test ToolkitRegistry.resolve()
        assert registry.resolve('to_smiles') == registry.registered_toolkits[0].to_smiles

        # Test ToolkitRegistry.call()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = registry.call('from_smiles', smiles)
        smiles2 = registry.call('to_smiles', molecule)
        assert smiles == smiles2

    @RDKitToolkitWrapper.requires_toolkit()
    #@pytest.mark.skipif( not RDKitToolkitWrapper.toolkit_is_available() )
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

    @RDKitToolkitWrapper.requires_toolkit()
    @AmberToolsToolkitWrapper.requires_toolkit()
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

