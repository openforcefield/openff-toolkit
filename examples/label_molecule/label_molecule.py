#!/bin/env python


from openforcefield.utils import *
from openforcefield.typing.engines.smirnoff import get_molecule_parameterIDs, ForceField

# Create an oemol
mol = OEMol()
OEParseSmiles(mol, 'CCC')
OEAddExplicitHydrogens(mol)


ff = ForceField( get_data_filename('forcefield/Frosst_AlkEthOH.ffxml') )

labels= ff.labelMolecules( [mol], verbose = True )
print labels
for mol_entry in range(len(labels)):
    for force in labels[mol_entry].keys():
        print("\n%s:" % force)
        for (atom_indices, pid, smirks) in labels[mol_entry][force]:
            atomstr=''
            for idx in atom_indices:
                atomstr += '%6s' % idx
            print("%s : %s \t smirks %s" % (atomstr, pid, smirks) )
