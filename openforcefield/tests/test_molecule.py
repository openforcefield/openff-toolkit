#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for Topology

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import pickle
from functools import partial
from unittest import TestCase

import pytest

from openforcefield import utils, topology
from openforcefield.topology import Molecule
from openforcefield.utils import get_data_filename
from openforcefield.utils import RDKIT_UNAVAILABLE, OPENEYE_UNAVAILABLE, SUPPORTED_TOOLKITS

#=============================================================================================
# TESTS
#=============================================================================================


_OPENEYE_UNAVAILABLE_MESSAGE = 'requires the OpenEye toolkit'
_RDKIT_UNAVAILABLE_MESSAGE = 'requires RDKit'

# TODO: Add tests comparing RDKit and OpenEye aromaticity perception

def assert_molecule_is_equal(molecule1, molecule2, msg):
    """Compare whether two Molecule objects are equal

    Parameters
    ----------
    molecule1, molecule2 : openforcefield.topology.Molecule
        Molecules to be compared
    msg : str
        Message to include if molecules fail to match.

    """
    # TODO:
    pass

def test_cheminformatics_toolkit_is_installed():
    """Ensure that at least one supported cheminformatics toolkit is installed."""
    if RDKIT_UNAVAILABLE and OPENEYE_UNAVAILABLE:
        msg = 'No supported cheminformatics toolkits are installed. Please install a supported toolkit:\n'
        msg += str(SUPPORTED_TOOLKITS)
        raise Exception(msg)

@pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE) # TODO: Remove this when setUp is openeye-independent
class TestMolecule(TestCase):
    from openforcefield.topology import Molecule

    @pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE) # TODO: Remove this when setUp is openeye-independent
    def setUp(self):
        # TODO: Serialize the offmols instead so that we can run this test without OpenEye
        #self.molecules = pickle.load('zinc-subset-offmols.pkl')
        filename = get_data_filename('molecules/zinc-subset-tripos.mol2.gz')
        self.molecules = Molecule.from_file(filename)

    def test_serialize(self):
        serialized = pickle.dumps(self.molecules)
        molecules_copy = pickle.loads(serialized)

    @pytest.mark.skipif(RDKIT_UNAVAILABLE, reason=_RDKIT_UNAVAILABLE_MESSAGE)
    def test_rdkit_roundtrip(self):
        for molecule in self.molecules:
            rdmol = molecule.to_rdkit()
            molecule2 = Molecule.from_rdmol(rdmol)
            assert_molecule_is_equal(molecule, molecule2, "Molecule.to_rdmol()/from_rdmol() round trip failed")
            molecule3 = Molecule(rdmol)
            assert_molecule_is_equal(molecule, molecule3, "Molecule(rdmol) constructor failed")

    @pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason='requires the OpenEye toolkit')
    def test_oemol_roundtrip(self):
        """Test creation of Molecule object from OpenEye OEMol
        """
        for molecule in self.molecules:
            oemol = molecule.to_openeye()
            molecule2 = Molecule.from_openeye(oemol)
            assert_molecule_is_equal(molecule, molecule2, "Molecule.to_openeye()/from_openeye() round trip failed")
            molecule3 = Molecule(oemol)
            assert_molecule_is_equal(molecule, molecule3, "Molecule(oemol) constructor failed")

    @pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE)
    def test_assign_partial_charges(self):
        """Test assignment of partial charges
        """
        # TODO: Only the OpenEye toolkit currently supports charge models
        # TODO: Does RDKit support other charge models?
        for molecule in self.molecules:
            for charge_model in topology.ALLOWED_CHARGE_MODELS:
                molecule.assign_partial_charges(method=charge_model)

    @pytest.mark.skipif(OPENEYE_UNAVAILABLE, reason=_OPENEYE_UNAVAILABLE_MESSAGE)
    def test_assign_fractional_bond_orders(self):
        """Test assignment of fractional bond orders
        """
        for molecule in self.molecules:
            for charge_model in topology.ALLOWED_FRACTIONAL_BONDORDER_MODELS:
                molecule.assign_fractional_bond_orders(method=charge_model)

    def test_bonds(self):
        """Test iteration over bonds
        """
        for molecule in self.molecules:
            bonds = molecule.bonds()
            # TODO: Check known cases

    def test_angles(self):
        """Test iteration over angles
        """
        for molecule in self.molecules:
            angles = molecule.angles()
            # TODO: Check known cases

    def test_torsions(self):
        """Test iteration over torsions
        """
        for molecule in self.molecules:
            torsion = molecule.torsions()
            # TODO: Check known cases

    def test_propers(self):
        """Test iteration over proper torsions
        """
        for molecule in self.molecules:
            torsion = molecule.propers()
            # TODO: Check known cases

    def test_impropers(self):
        """Test iteration over improper torsions
        """
        for molecule in self.molecules:
            torsion = molecule.impropers()
            # TODO: Check known cases
