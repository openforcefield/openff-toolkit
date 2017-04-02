from functools import partial
from unittest import TestCase

from openforcefield import utils

class TestUtils(TestCase):
    def test_read_molecules(self):
        molecules = utils.read_molecules('zinc-subset-tripos.mol2.gz', verbose=False)
    def test_positions(self):
        """Test ability to extract and set positions."""
        molecules = utils.read_molecules('zinc-subset-tripos.mol2.gz', verbose=False)
        positions = utils.extractPositionsFromOEMol(molecules[0])
        utils.setPositionsInOEMol(molecules[0], positions)
