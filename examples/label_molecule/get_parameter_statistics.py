#!/bin/env python

from openforcefield.utils import get_data_filename
#from openforcefield.typing.engines.smirnoff import get_molecule_parameterIDs
from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField

#oemols = read_molecules(get_data_filename('molecules/AlkEthOH_test_filt1_tripos.mol2'))
mols = Molecule.from_file(get_data_filename('molecules/AlkEthOH_test_filt1_tripos.mol2'))
ff = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.offxml'))
#ffxml = get_data_filename('forcefield/Frosst_AlkEthOH.offxml')



parameters_by_molecule, parameters_by_ID = get_molecule_parameterIDs( oemols, ffxml)

# Print some info
print "Parameters by molecule:"
for smi in parameters_by_molecule.keys():
    print smi, parameters_by_molecule[smi]


print "Molecules with parameter IDs:"
for pid in parameters_by_ID.keys():
    print pid, parameters_by_ID[pid]
