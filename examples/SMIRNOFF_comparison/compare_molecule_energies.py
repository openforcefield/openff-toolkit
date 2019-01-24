#!/bin/env python
import os

# Cross-check energies of molecules from AlkEthOH set using SMIRNOFF xml file
# versus energies from AMBER .prmtop and .crd files (parm@frosst params)

datapath = './AlkEthOH_tripos/AlkEthOH_rings_filt1'
#molname = 'AlkEthOH_r0' #That fails, but it's complicated. Try cyclobutane
molname = 'AlkEthOH_r51'
mol_filepath = os.path.join(datapath, molname + '_tripos.mol2')
prmtop_filepath = os.path.join(datapath, molname + '.top')
inpcrd_filepath = os.path.join(datapath, molname + '.crd')

# Check if we have this data file; if not we have to extract the archive.
if not os.path.isdir(datapath):
    print("Extracting archived molecule files.")
    # Extract the AlkEthOH dataset shipped with the toolkit in data/molecules/ in the working directory.
    from openforcefield.tests.utils import get_data_filename
    tarfile_path = os.path.join(get_data_filename('molecules'), 'AlkEthOH_tripos.tar.gz')
    import tarfile
    with tarfile.open(tarfile_path, 'r:gz') as tar:
        tar.extractall()

# Load molecule
from openforcefield.topology import Molecule
molecule = Molecule.from_file(mol_filepath)

# Load forcefield
from openforcefield.typing.engines.smirnoff import ForceField
forcefield = ForceField('Frosst_AlkEthOH_parmAtFrosst.offxml')

# Compare energies
from openforcefield.tests.utils import compare_amber_smirnoff
results = compare_amber_smirnoff(prmtop_filepath, inpcrd_filepath, forcefield, molecule)

print(results)
