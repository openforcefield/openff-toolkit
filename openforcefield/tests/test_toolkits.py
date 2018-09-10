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
    # TODO: Test that exceptions are properly raised by supported methods
    pass

class TestOpenEyeToolkitWrapper(TestCase):
    """Test the OpenEyeToolkitWrapper"""

    @requires_openeye('oechem')
    def test_to_openeye(self):
        """Test OpenEyeToolkitWrapper.to_openeye()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        molecule = Molecule.from_smiles('CC')
        oemol = toolkit_wrapper.to_openeye(molecule)
        # TODO: Assert equality

    @requires_openeye('oechem')
    def test_from_openeye(self):
        """Test OpenEyeToolkitWrapper.from_openeye()"""
        toolkit_wrapper = OpenEyeToolkitWrapper()
        molecule = Molecule.from_smiles('CC')
        oemol = toolkit_wrapper.to_openeye(molecule)
        molecule2 = toolkit_wrapper.from_openeye(oemol)
        assert molecule == molecule2

    @requires_openeye('oechem')
    def test_to_smiles(self):
        """Test OpenEyeToolkitWrapper.to_smiles()"""
        smiles = 'CC'
        molecule = Molecule.from_smiles(smiles)
        smiles2 = toolkit_wrapper.to_smiles(molecule)
        assert smiles == smiles2

class RDKitToolkitWrapper(TestCase):
    """Test the RDKitToolkitWrapper"""
    pass

class TestAmberToolsWrapper(TestCase):
    """Test the AmberToolsWraper"""
    pass

class TestToolkitRegistry(TestCase):
    """Test the ToolkitRegistry"""
    pass
