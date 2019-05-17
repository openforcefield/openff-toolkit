#!/usr/bin/env python

"""
Cross-check energies of molecules from AlkEthOH set using SMIRNOFF xml file
versus energies from AMBER .prmtop and .crd files (parm@frosst params).
"""

import os

# datapath = './AlkEthOH_tripos/AlkEthOH_chain_filt1'
# datapath = './AlkEthOH_tripos/AlkEthOH_rings_filt1'
datapath = './AlkEthOH_tripos/AlkEthOH_test_filt1'

molname = 'AlkEthOH_r22'

mol_filepath = os.path.join(datapath, molname + '_tripos.mol2')
prmtop_filepath = os.path.join(datapath, molname + '.top')
inpcrd_filepath = os.path.join(datapath, molname + '.crd')

# Check if we have this data file; if not we have to extract the archive.
if not os.path.isdir(datapath):
    print("Extracting archived molecule files.")
    # Extract the AlkEthOH dataset shipped with the toolkit in data/molecules/ in the working directory.
    from openforcefield.tests.utils import get_data_file_path
    tarfile_path = os.path.join(get_data_file_path('molecules'), 'AlkEthOH_tripos.tar.gz')
    import tarfile
    with tarfile.open(tarfile_path, 'r:gz') as tar:
        tar.extractall()

# Load molecule
from openforcefield.topology import Molecule
molecule = Molecule.from_file(mol_filepath)

# Load forcefield
from openforcefield.typing.engines.smirnoff import ForceField
forcefield = ForceField('test_forcefields/Frosst_AlkEthOH_parmAtFrosst.offxml')

# Compare energies
from openforcefield.tests.utils import compare_amber_smirnoff
# We ignore the charges as they are not included in the force field.
# TODO: Reactivate this check when we'll be able to load charges from the file.
energies = compare_amber_smirnoff(prmtop_filepath, inpcrd_filepath,
                                  forcefield, molecule,
                                  ignore_charges=True)

# Pretty-print the result.
from pprint import pprint
pprint(energies)
