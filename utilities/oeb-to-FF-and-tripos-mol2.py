#!/usr/bin/env python
"""
Convert oeb file of molecules into Forcefield atom typed and Tripos atom typed mol2 files.

Example:

> python ../oeb-to-FF-and-tripos-mol2.py test_filt1
"""
from __future__ import division
from __future__ import print_function
import os,sys
import openeye.oechem as oechem

def main(argv=sys.argv):
    if len(argv) != 2:
        oechem.OEThrow.Usage("%s <infile (oeb file prefix)>" % argv[0])

    ifs = oechem.oemolistream()
    if not ifs.open(argv[1]+'.oeb'):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % argv[1]+'.oeb' )

    ofsff = oechem.oemolostream()
    ofsff.SetFlavor( oechem.OEFormat_MOL2, oechem.OEOFlavor_MOL2_Forcefield )
    if not ofsff.open(argv[1]+'_ff.mol2'):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % argv[1]+'_ff.mol2')

    ofsTri = oechem.oemolostream()
    ofsTri.SetFlavor( oechem.OEFormat_MOL2, oechem.OEOFlavor_MOL2_Forcefield )
    if not ofsTri.open(argv[1]+'_tripos.mol2'):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % argv[1]+'_tripos.mol2')

    for mol in ifs.GetOEMols():
        oechem.OETriposAtomNames(mol)
        oechem.OEWriteConstMolecule(ofsff, mol)
        oechem.OETriposAtomTypeNames(mol)
        oechem.OEWriteConstMolecule(ofsTri, mol)

    ifs.close()
    ofsff.close()
    ofsTri.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))#!/usr/bin/env python
