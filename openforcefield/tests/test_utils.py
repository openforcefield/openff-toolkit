from functools import partial
from unittest import TestCase

from openforcefield import utils

class TestUtils(TestCase):
    def test_read_molecules(self):
        molecules = utils.read_molecules('zinc-subset-tripos.mol2.gz', verbose=False)
    def test_parse_odds_file(self):
        odds = utils.parse_odds_file('odds_files/atom_index_odds.smarts')
        odds = utils.parse_odds_file('odds_files/bond_OR_bases.smarts')
        self.assertIsNone(odds[1], msg = "Parsing odds file with no odds should give None as the second entry")
    def test_positions(self):
        """Test ability to extract and set positions."""
        molecules = utils.read_molecules('zinc-subset-tripos.mol2.gz', verbose=False)
        positions = utils.extractPositionsFromOEMol(molecules[0])
        utils.setPositionsInOEMol(molecules[0], positions)
