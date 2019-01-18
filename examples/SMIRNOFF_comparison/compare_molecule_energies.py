#!/bin/env python
import os

# Cross-check energies of molecules from AlkEthOH set using SMIRNOFF xml file
# versus energies from AMBER .prmtop and .crd files (parm@frosst params)

datapath = './AlkEthOH_inputfiles/AlkEthOH_rings_filt1'
#molname = 'AlkEthOH_r0' #That fails, but it's complicated. Try cyclobutane
molname = 'AlkEthOH_r51'
mol_filepath = os.path.join(datapath, molname + '.mol2')
prmtop_filepath = os.path.join(datapath, molname + '.top')
inpcrd_filepath = os.path.join(datapath, molname + '.crd')

# Check if we have this data file; if not we have to extract the archive.
if not os.path.isfile(mol_filepath):
    print("Extracting archived molecule files.")
    tarfile = 'AlkEthOH_inputfiles.tar.gz'
    os.system('tar -xf AlkEthOH_inputfiles.tar.gz')

# Load molecule
from openforcefield.topology import Molecule
molecule = Molecule.from_file(mol_filepath)

# Load forcefield
from openforcefield.typing.engines.smirnoff import ForceField
forcefield = ForceField('forcefield/Frosst_AlkEthOH_parmAtFrosst.offxml')

# Compare energies
from openforcefield.tests.utils import compare_molecule_energies
results = compare_molecule_energies(prmtop_filepath, inpcrd_filepath, forcefield, molecule)
