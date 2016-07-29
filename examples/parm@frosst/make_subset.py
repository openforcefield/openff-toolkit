#!/bin/env python

"""Take the ZINC subset here and make a smaller subset of it for testing purposes."""

from openeye.oechem import *

nmols = 500 #Number of molecules to retain out of full ~7500
# Currently the above are taken as the first 500. We could also take randomly.



# Read set with tripos types, write subset
ifs = oemolistream( 'molecules/zinc-subset-tripos.mol2.gz')
ofs = oemolostream( 'molecules/zinc-subset-%s-tripos.mol2.gz' % nmols )
mol = OEMol()
ct=0
while OEReadMolecule(ifs, mol) and ct < nmols:
    OEWriteConstMolecule(ofs, mol)
    ct += 1


# Read set with parm@frosst types, write subset
# Use flavors here to ensure writing doesn't mangle atom types
ifs = oemolistream( 'molecules/zinc-subset-parm@frosst.mol2.gz')
flavor = OEIFlavor_Generic_Default | OEIFlavor_MOL2_Default | OEIFlavor_MOL2_Forcefield
ifs.SetFlavor(OEFormat_MOL2, flavor)
ofs = oemolostream( 'molecules/zinc-subset-%s-parm@frosst.mol2.gz' % nmols )
ofs.SetFlavor(OEFormat_MOL2, flavor)
mol = OEMol()
ct=0
while OEReadMolecule(ifs, mol) and ct < nmols:
    OEWriteConstMolecule(ofs, mol)
    ct+=1

