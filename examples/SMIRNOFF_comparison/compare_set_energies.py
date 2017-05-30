#!/bin/env python

from openeye import oechem
from openforcefield.utils import get_data_filename
from openforcefield.typing.engines.smirnoff import ForceField, compare_molecule_energies

import os
import glob

# Cross-check energies of molecules from AlkEthOH set using SMIRNOFF xml file
# versus energies from AMBER .prmtop and .crd files (parm@frosst params)

#datapath = './AlkEthOH_inputfiles/AlkEthOH_chain_filt1'
datapath = './AlkEthOH_inputfiles/AlkEthOH_rings_filt1'

#Check if this is a directory, if not, extract it
# Check if we have this data file; if not we have to extract the archive
if not os.path.isdir( datapath ):
    print "Extracting archived molecule files."
    tarfile = 'AlkEthOH_inputfiles.tar.gz'
    os.system('tar -xf AlkEthOH_inputfiles.tar.gz')

#Obtain list of molecules
mol_filenames = glob.glob( datapath+'/*.mol2')
mol_filenames = [ fnm for fnm in mol_filenames if not 'c1302' in fnm ] #Skip water

for mol_filename in mol_filenames:
    molname = os.path.basename( mol_filename).replace('.mol2','')
    print("Comparing %s (%s/%s)..." % (molname, mol_filenames.index(mol_filename), len(mol_filenames) ) )

    # Load OEMol
    mol = oechem.OEGraphMol()
    ifs = oechem.oemolistream(mol_filename)
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
    oechem.OEReadMolecule(ifs, mol )
    oechem.OETriposAtomNames(mol)

    # Load forcefield
    forcefield = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH_parmAtFrosst.ffxml'))

    # Compare energies
    prmtop = os.path.join( datapath, molname+'.top')
    crd = os.path.join( datapath, molname+'.crd')
    results = compare_molecule_energies( prmtop, crd, forcefield, mol)
