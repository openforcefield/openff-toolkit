#!/bin/env python

from smarty.utils import *
from smarty.forcefield import *
from smarty.forcefield_utils import get_molecule_parameterIDs

oemols = read_molecules(get_data_filename('molecules/AlkEtOH-tripos.mol2.gz'))
ffxml = get_data_filename('forcefield/Frosst_AlkEtOH.ffxml')

parameters_by_molecule, parameters_by_ID = get_molecule_parameterIDs( oemols, ffxml)

# Print some info
print "Parameters by molecule:"
for smi in parameters_by_molecule.keys():
    print smi, parameters_by_molecule[smi]


print "Molecules with parameter IDs:"
for pid in parameters_by_ID.keys():
    print pid, parameters_by_ID[pid]
