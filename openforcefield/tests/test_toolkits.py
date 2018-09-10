#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for cheminformatics toolkit interfaces

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

from functools import partial
from unittest import TestCase
from openforcefield import utils

from openforcefield.utils.toolkits import ToolkitWrapper, OpenEyeToolkitWrapper, RDKitToolkitWrapper, AmberToolsWrapper, ToolkitRegistry

#=============================================================================================
# TESTS
#=============================================================================================

class TestToolkitWrapper(TestCase):
    """Test the ToolkitWrapper abstract base class"""

    def test_compute_partial_charges(self):
        """Test that exceptions are properly raised by supported methods"""
        with self.assertRaises(NotImplementedException) as context:
            wrapper = ToolkitWrapper()
            wrapper.compute_partial_charges()

    def test_to_smiles(self):
        """Test that exceptions are properly raised by supported methods"""
        with self.assertRaises(NotImplementedException) as context:
            wrapper = ToolkitWrapper()
            wrapper.to_smiles()

    def test_from_smiles(self):
        """Test that exceptions are properly raised by supported methods"""
        with self.assertRaises(NotImplementedException) as context:
            wrapper = ToolkitWrapper()
            wrapper.from_smiles()

class TestOpenEyeToolkitWrapper(TestCase):
    """Test the OpenEyeToolkitWrapper"""

    @requires_openeye('oechem') # TODO: Is this needed?
    def test_smiles(self):
        """Test OpenEyeToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = 'CC'
        molecule = toolkit_wrapper.from_smiles(smiles)
        smiles2 = molecule.to_smiles()
        assert smiles == smiles2

    @requires_openeye('oechem') # TODO: Is this needed?
    def test_openeye(self):
        """Test OpenEyeToolkitWrapper to_openeye() and from_openeye()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = 'CC'
        molecule = Molecule.from_smiles(smiles)
        oemol = toolkit_wrapper.to_openeye(molecule)
        molecule2 = toolkit_wrapper.from_openeye(oemol)
        smiles2 = molecule2.to_smiles()
        assert smiles == smiles2

    @requires_openeye('oechem') # TODO: Is this needed?
    def test_compute_partial_charges(self):
        """Test OpenEyeToolkitWrapper compute_partial_charges()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        smiles = 'CC'
        molecule = toolkit_wrapper.from_smiles(smiles)
        # TODO: Test all supported charge models
        partial_charges = toolkit_wrapper.compute_partial_charges(molecule)
        # TODO: Check partial charge invariants (total charge, charge equivalence)

class RDKitToolkitWrapper(TestCase):
    """Test the RDKitToolkitWrapper"""

    @requires_rdkit() # TODO: Is this needed?
    def test_smiles(self):
        """Test RDKitToolkitWrapper to_smiles() and from_smiles()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        smiles = 'CC'
        molecule = toolkit_wrapper.from_smiles(smiles)
        smiles2 = molecule.to_smiles()
        assert smiles == smiles2

    @requires_rdkit() # TODO: Is this needed?
    def test_rdkit(self):
        """Test RDKitToolkitWrapper to_rdkit() and from_rdkit()"""
        toolkit_wrapper = RDKitToolkitWrapper()
        smiles = 'CC'
        molecule = Molecule.from_smiles(smiles)
        rdmol = toolkit_wrapper.to_openeye(molecule)
        molecule2 = toolkit_wrapper.from_openeye(rdmol)
        smiles2 = molecule2.to_smiles()
        assert smiles == smiles2

class TestAmberToolsWrapper(TestCase):
    """Test the AmberToolsWraper"""

    @requires_ambertools # TODO: Is this needed?
    @requires_rdkit # TODO: Is this needed?
    def test_compute_partial_charges(self):
        """Test AmberTools compute_partial_charges()"""
        rdkit_toolkit_wrapper = RDKitToolkitWrapper()
        amber_toolkit_wrapper = AmberToolsToolkitWrapper()
        smiles = 'CC'
        molecule = rdkit_toolkit_wrapper.from_smiles(smiles)
        # TODO: Test all supported charge models
        partial_charges = amber_toolkit_wrapper.compute_partial_charges(molecule)
        # TODO: Check partial charge invariants (total charge, charge equivalence)

class TestToolkitRegistry(TestCase):
    """Test the ToolkitRegistry"""

    @requires_openeye('oechem') # TODO: Is this needed?
    def test_register_openeye(self):
        """Test creation of toolkit registry with OpenEye toolkit"""
        # Test registration of OpenEyeToolkitWrapper
        registry = ToolkitRegistry(register_imported_toolkit_wrappers=False)
        registry.register_toolkit(OpenEyeToolkitWrapper)
        assert set([ type(c) for c in registry.registered_toolkits]) == set([OpenEyeToolkitWrapper])

        # Test ToolkitRegistry.resolve()
        assert registry.resolve('to_smiles') == registry.registered_toolkits[0].to_smiles

        # Test ToolkitRegistry.call()
        smiles = 'CC'
        molecule = registry.call('from_smiles', smiles)
        smiles2 = registry.call('to_smiles', molecule)
        assert smiles == smiles2

    @requires_rdkit # TODO: Is this needed?
    def test_register_rdkit(self):
        """Test creation of toolkit registry with RDKit toolkit"""
        # Test registration of RDKitToolkitWrapper
        registry = ToolkitRegistry(register_imported_toolkit_wrappers=False)
        registry.register_toolkit(RDKitToolkitWrapper)
        assert set([ type(c) for c in registry.registered_toolkits]) == set([RDKitToolkitWrapper])

        # Test ToolkitRegistry.resolve()
        assert registry.resolve('to_smiles') == registry.registered_toolkits[0].to_smiles

        # Test ToolkitRegistry.call()
        smiles = 'CC'
        molecule = registry.call('from_smiles', smiles)
        smiles2 = registry.call('to_smiles', molecule)
        assert smiles == smiles2

    @requires_ambertools # TODO: Is this needed?
    @requires_rdkit # TODO: Is this needed?
    def test_register_ambertools(self):
        """Test creation of toolkit registry with RDKit toolkit"""
        # Test registration of AmberToolsToolkitWrapper
        registry = ToolkitRegistry(register_imported_toolkit_wrappers=False)
        registry.register_toolkit(AmberToolsToolkitWrapper)
        assert set([ type(c) for c in registry.registered_toolkits]) == set([AmberToolsToolkitWrapper])

        # Test ToolkitRegistry.resolve()
        assert registry.resolve('compute_partial_charges') == registry.registered_toolkits[0].compute_partial_charges

        # Test ToolkitRegistry.call()
        registry.register_toolkit(RDKitToolkitWrapper)
        smiles = 'CC'
        molecule = registry.call('from_smiles', smiles)
        partial_charges = registry.call('compute_partial_charges', molecule)
