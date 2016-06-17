from functools import partial
from smarty import AtomTyper, AtomTypeSampler
import smarty
from smarty.utils import get_data_filename
from unittest import TestCase

class TestAtomTypeSampler(TestCase):
    def test_atomtyper(self):
        basetypes_filename = get_data_filename('atomtypes/basetypes.smarts')
        decorators_filename = get_data_filename('atomtypes/decorators.smarts')
        replacements_filename = get_data_filename('atomtypes/replacements.smarts')
        molecules = smarty.utils.read_molecules(get_data_filename('molecules/zinc-subset-tripos.mol2.gz'), verbose=False)
        reference_typed_molecules = smarty.utils.read_molecules(get_data_filename('molecules/zinc-subset-parm@frosst.mol2.gz'), verbose=False)

        # Construct atom type sampler.
        atomtype_sampler = smarty.AtomTypeSampler(molecules, basetypes_filename, decorators_filename, replacements_filename=replacements_filename, reference_typed_molecules=reference_typed_molecules, verbose=False)

        # Start sampling atom types.
        atomtype_sampler.run(2)
