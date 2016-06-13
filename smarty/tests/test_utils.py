from functools import partial
from smarty.utils import get_data_filename
from unittest import TestCase

import smarty

class TestUtils(TestCase):
    def test_read_molecules(self):
        molecules = smarty.utils.read_molecules(get_data_filename('molecules/zinc-subset-tripos.mol2.gz'), verbose=False)
