#!/usr/bin/env python
"""
Convert oeb file of molecules into force field atom typed and Tripos atom typed mol2 files.

Example:

> python ../oeb-to-FF-and-tripos-mol2.py test_filt1
"""
from __future__ import division
from __future__ import print_function
import sys
import openeye.oechem as oechem

def main(argv=sys.argv):
    if len(argv) != 2:
        oechem.OEThrow.Usage(f"{argv[0]} <infile (oeb file prefix)>")

    ifs = oechem.oemolistream()
    if not ifs.open(f"{argv[1]}.oeb"):
        oechem.OEThrow.Fatal(f"Unable to open {argv[1]}.oeb for reading")

    ofsff = oechem.oemolostream()
    ofsff.SetFlavor( oechem.OEFormat_MOL2, oechem.OEOFlavor_MOL2_Forcefield )
    if not ofsff.open(f"{argv[1]}_ff.mol2"):
        oechem.OEThrow.Fatal(f"Unable to open {argv[1]}_ff.mol2 for writing")

    ofsTri = oechem.oemolostream()
    ofsTri.SetFlavor( oechem.OEFormat_MOL2, oechem.OEOFlavor_MOL2_Forcefield )
    if not ofsTri.open(f"{argv[1]}_tripos.mol2"):
        oechem.OEThrow.Fatal(f"Unable to open {argv[1]}_tripos.mol2 for writing")

    for mol in ifs.GetOEMols():
        oechem.OETriposAtomNames(mol)
        oechem.OEWriteConstMolecule(ofsff, mol)
        oechem.OETriposAtomTypeNames(mol)
        oechem.OEWriteConstMolecule(ofsTri, mol)

    ifs.close()
    ofsff.close()
    ofsTri.close()

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
