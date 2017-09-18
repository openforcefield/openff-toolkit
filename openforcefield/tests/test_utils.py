from functools import partial
from unittest import TestCase
import parmed
from openforcefield import utils

class TestUtils(TestCase):
    def test_read_molecules(self):
        molecules = utils.read_molecules('zinc-subset-tripos.mol2.gz', verbose=False)
    def test_positions(self):
        """Test ability to extract and set positions."""
        molecules = utils.read_molecules('zinc-subset-tripos.mol2.gz', verbose=False)
        positions = utils.extractPositionsFromOEMol(molecules[0])
        utils.setPositionsInOEMol(molecules[0], positions)
    def test_smirnoff_structure(self):
       """Test function that generates the SMIRNOFF parmed.Structure."""
       molecule = utils.read_molecules('toluene.pdb', verbose=False)[0]
       molecule_structure = utils.generateSMIRNOFFStructure(molecule)
       self.assertIsInstance(molecule_structure, parmed.structure.Structure)
    def test_protein_structure(self):
        from simtk.openmm import app
        pdbfile = utils.get_data_filename('proteins/T4-protein.pdb')
        proteinpdb = app.PDBFile(pdbfile)
        protein_structure = utils.generateProteinStructure(proteinpdb)
        self.assertIsInstance(protein_structure, parmed.structure.Structure)
