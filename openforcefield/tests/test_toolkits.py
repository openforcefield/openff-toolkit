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



import pytest
from openforcefield.utils.toolkits import ToolkitWrapper, OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsToolkitWrapper, ToolkitRegistry

## From Jeff: Where should I get molecules for ToolkitWrapper tests?
from openforcefield.topology.molecule import Molecule

#=============================================================================================
# TESTS
#=============================================================================================

class TestToolkitWrapper(TestCase):
    """Test the ToolkitWrapper abstract base class"""

    def test_compute_partial_charges(self):
        """Test that exceptions are properly raised by supported methods"""
        molecule = Molecule()
        with self.assertRaises(NotImplementedError) as context:
            toolkit_wrapper = ToolkitWrapper()
            toolkit_wrapper.compute_partial_charges(molecule)

    def test_to_smiles(self):
        """Test that exceptions are properly raised by supported methods"""
        molecule = Molecule()
        with self.assertRaises(NotImplementedError) as context:
            toolkit_wrapper = ToolkitWrapper()
            toolkit_wrapper.to_smiles(molecule)

    def test_from_smiles(self):
        """Test that exceptions are properly raised by supported methods"""
        with self.assertRaises(NotImplementedError) as context:
            toolkit_wrapper = ToolkitWrapper()
            toolkit_wrapper.from_smiles('CC')

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
    def test_openeye(self):
        """Test OpenEyeToolkitWrapper to_openeye() and from_openeye()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        # This differs from RDKit's expected output due to different canonicalization schemes
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles)
        oemol = molecule.to_openeye()
        molecule2 = Molecule.from_openeye(oemol)
        smiles2 = molecule2.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2

    #@pytest.mark.skipif( not OpenEyeToolkitWrapper.toolkit_is_available() )
    @OpenEyeToolkitWrapper.requires_toolkit()
    def test_compute_partial_charges(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = toolkit_wrapper.from_smiles(smiles)
        # TODO: Test all supported charge models
        partial_charges = toolkit_wrapper.compute_partial_charges(molecule)
        # TODO: Check partial charge invariants (total charge, charge equivalence)

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
    def test_rdkit(self):
        """Test RDKitToolkitWrapper to_rdkit() and from_rdkit()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        # This differs from OE's expected output due to different canonicalization schemes
        smiles = '[H][C]([H])([H])[C]([H])([H])[H]'
        molecule = Molecule.from_smiles(smiles,
                                        toolkit_registry=toolkit_wrapper)
        rdmol = molecule.to_rdkit()
        molecule2 = Molecule.from_rdkit(rdmol)
        smiles2 = molecule2.to_smiles(toolkit_registry=toolkit_wrapper)
        assert smiles == smiles2

class TestAmberToolsWrapper(TestCase):
    """Test the AmberToolsWraper"""

    @RDKitToolkitWrapper.requires_toolkit()
    @AmberToolsToolkitWrapper.requires_toolkit()
    def test_compute_partial_charges(self):
        """Test AmberTools compute_partial_charges()"""
        rdkit_toolkit_wrapper = RDKitToolkitWrapper()
        amber_toolkit_wrapper = AmberToolsToolkitWrapper()
        smiles = '[H]C([H])([H])C([H])([H])[H]'
        molecule = rdkit_toolkit_wrapper.from_smiles(smiles)
        # TODO: Test all supported charge models
        partial_charges = amber_toolkit_wrapper.compute_partial_charges(molecule)
        # TODO: Check partial charge invariants (total charge, charge equivalence)

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
