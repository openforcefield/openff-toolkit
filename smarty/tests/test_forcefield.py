from functools import partial
from smarty import ForceField
import smarty
from smarty.utils import get_data_filename
from unittest import TestCase

class TestForceField(TestCase):
    def test_read_ffxml(self):
        forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.ffxml'))
