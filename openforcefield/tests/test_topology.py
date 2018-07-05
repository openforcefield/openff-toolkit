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

from openforcefield import utils, topology

from pytest.mark import skipif
from openforcefield.utils import RDKIT_UNAVAILABLE, OPENEYE_UNAVAILABLE, SUPPORTED_TOOLKITS
from openforcefield.tests.utils.utils import get_amber_system, get_packmol_pdbfile, get_monomer_mol2file

#=============================================================================================
# TESTS
#=============================================================================================

def test_cheminformatics_toolkit_is_installed():
    """Ensure that at least one supported cheminformatics toolkit is installed."""
    if RDKIT_UNAVAILABLE and OPENEYE_UNAVAILABLE:
        msg = 'No supported cheminformatics toolkits are installed. Please install a supported toolkit:\n'
        msg += str(SUPPORTED_TOOLKITS)
        raise Exception(msg)

class TestTopology(TestCase):
    from openforcefield.topology import Topology

    def setUp(self):
        pass

    def test_from_openmm(self):
        """Test creation of an openforcefield Topology object from an OpenMM Topology and component molecules"""
        pdbfile = app.PDBFile(get_packmol_pdbfile('cyclohexane_ethanol_0.4_0.6.pdb'))
        molecules = [ Molecule.from_file(get_monomer_mol2file(name)) for name in ('ethanol', 'cyclohexane') ]
        topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
