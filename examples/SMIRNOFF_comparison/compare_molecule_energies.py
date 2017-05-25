#!/bin/env python

from openeye import oechem
from openforcefield.utils import get_data_filename
from openforcefield.typing.engines.smirnoff import ForceField, compare_molecule_energies

import os

# Cross-check energies of molecules from AlkEthOH set using SMIRNOFF xml file
# versus energies from AMBER .prmtop and .crd files (parm@frosst params)

datapath = './AlkEthOH_inputfiles/AlkEthOH_rings_filt1'
#molname = 'AlkEthOH_r0' #That fails, but it's complicated. Try cyclobutane
molname = 'AlkEthOH_r51'
mol_filename = os.path.join( datapath, molname+'.mol2')

# Check if we have this data file; if not we have to extract the archive
if not os.path.isfile( mol_filename):
    print "Extracting archived molecule files."
    tarfile = 'AlkEthOH_inputfiles.tar.gz'
    os.system('tar -xf AlkEthOH_inputfiles.tar.gz')

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
