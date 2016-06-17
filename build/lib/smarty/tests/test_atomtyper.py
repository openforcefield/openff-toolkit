from functools import partial
from smarty import AtomTyper
import smarty
from smarty.utils import get_data_filename
from unittest import TestCase

class TestAtomTyper(TestCase):
    def test_read_typelist(self):
        atomtypes = AtomTyper.read_typelist(get_data_filename('atomtypes/basetypes.smarts'))
        decorators = AtomTyper.read_typelist(get_data_filename('atomtypes/decorators.smarts'))
        replacements = AtomTyper.read_typelist(get_data_filename('atomtypes/replacements.smarts'))

    def test_atomtyper(self):
        typetag = 'atomtype'
        atomtypes = AtomTyper.read_typelist(get_data_filename('atomtypes/basetypes.smarts'))
        replacements = AtomTyper.read_typelist(get_data_filename('atomtypes/replacements.smarts'))
        molecules = smarty.utils.read_molecules(get_data_filename('molecules/zinc-subset-tripos.mol2.gz'), verbose=False)

        atomtyper = AtomTyper(atomtypes, typetag, replacements=replacements)
        for molecule in molecules:
            atomtyper.assignTypes(molecule)
