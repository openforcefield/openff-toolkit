#!/bin/env python

from simtk.openmm import app
import openforcefield.utils as utils
from openeye import oechem

molecule = oechem.OEMol()
molpdb = utils.get_data_filename('molecules/toluene.pdb')
with oechem.oemolistream(molpdb) as ifs:
    oechem.OEReadMolecule(ifs, molecule)

molecule_structure = utils.generateSMIRNOFFStructure(molecule)
print('Molecule:', molecule_structure)

pdbfile = utils.get_data_filename('proteins/T4-protein.pdb')
proteinpdb = app.PDBFile(pdbfile)
protein_structure = utils.generateProteinStructure(proteinpdb)
print('Protein:', protein_structure)

structure = utils.mergeStructure(protein_structure, molecule_structure)
print('Complex:', structure)
