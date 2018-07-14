#!/bin/env python
import os
import glob

# Cross-check energies of molecules from AlkEthOH set using SMIRNOFF xml file
# versus energies from AMBER .prmtop and .crd files (parm@frosst params)

#datapath = './AlkEthOH_inputfiles/AlkEthOH_chain_filt1'
datapath = './AlkEthOH_inputfiles/AlkEthOH_rings_filt1'

from openforcefield.utils import get_data_filename
offxml_filename = get_data_filename('forcefield/Frosst_AlkEthOH_parmAtFrosst.offxml')

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

    prmtop_filename = os.path.join(datapath, molname+'.top')
    inpcrd_filename = os.path.join(datapath, molname+'.crd')

    # Load molecule
    from openmmtools.topology import Molecule
    molecule = Molecule.from_file(mol_filename)

    # Load forcefield
    from openforcefield.typing.engines.smirnoff import ForceField
    forcefield = ForceField(offxml_filename)

    # Compare energies
    from openforcefield.tests.utils import compare_molecule_energies
    results = compare_molecule_energies(prmtop_filename, inpcrd_filename, forcefield, molecule)
