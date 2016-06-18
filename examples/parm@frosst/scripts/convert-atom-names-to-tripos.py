#!/usr/bin/env python
"""
Convert file of molecules from forcefield atom types to Tripos atom types.

Example:

> python ../convert-atom-names-to-tripos.py zinc-subset-parm@frosst.mol2.gz zinc-subset-tripos.mol2.gz
"""
################################################################
#  Copyright (C) 2006-2015 OpenEye Scientific Software, Inc.
################################################################
from __future__ import division
from __future__ import print_function
import os,sys
import openeye.oechem as oechem

def main(argv=sys.argv):
    if len(argv) != 3:
        oechem.OEThrow.Usage("%s <infile (forcefield types)> <outfile (Tripos types)>" % argv[0])

    ifs = oechem.oemolistream()
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor(oechem.OEFormat_MOL2, flavor)
    if not ifs.open(argv[1]):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % argv[1])

    ofs = oechem.oemolostream()
    if not ofs.open(argv[2]):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % argv[2])

    for mol in ifs.GetOEMols():
        oechem.OETriposAtomNames(mol)
        oechem.OEWriteConstMolecule(ofs, mol)

    ifs.close()
    ofs.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))#!/usr/bin/env python
