#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Tests for utility methods that involve parmed Structure manipulation

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

from functools import partial
from unittest import TestCase
import parmed
from openforcefield import utils

#=============================================================================================
# TESTS
#=============================================================================================

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

def test_merge_system():
    """Test merging of a system created from AMBER and another created from SMIRNOFF."""

    # Create System from AMBER
    # TODO: Add create_system_from_amber convenience function
    from openforcefield.typing.engines.smirnoff import create_system_from_amber
    prmtop_filename, inpcrd_filename = get_amber_filepaths('cyclohexane_ethanol_0.4_0.6')
    topology0, system0, positions0 = create_system_from_amber(prmtop_filename, inpcrd_filename)

    # TODO:
    from openeye import oechem
    # Load simple OEMol
    ifs = oechem.oemolistream(get_testdata_filename('molecules/AlkEthOH_c100.mol2'))
    mol = oechem.OEMol()
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor( oechem.OEFormat_MOL2, flavor)
    oechem.OEReadMolecule(ifs, mol)
    oechem.OETriposAtomNames(mol)

    # Load forcefield file
    forcefield = ForceField(AlkEthOH_offxml_filename)
    topology1, system1, positions1 = create_system_from_molecule(forcefield, mol)

    merge_system( topology0, topology1, system0, system1, positions0, positions1, verbose=True )

def test_component_combination():
    """Test that a system still yields the same energy after rebuilding it out of its components
    """

    # We've had issues where subsequent instances of a molecule might have zero charges
    # Here we'll try to catch this (and also explicitly check the charges) by re-building
    # a system out of its components

    # Create an OpenMM System from a PDB file containing a cyclohexane-ethanol mixture
    forcefield = ForceField(AlkEthOH_offxml_filename)
    pdbfile = app.PDBFile(get_packmol_pdbfile('cyclohexane_ethanol_0.4_0.6.pdb'))
    molecules = [ Molecule.from_file(get_monomer_mol2file(name)) for name in ('ethanol', 'cyclohexane') ]
    topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
    system = forcefield.createSystem(topology)

    # Convert System to a ParmEd Structure
    structure = parmed.openmm.topsystem.load_topology(topology.to_openmm(), system, pdbfile.positions)

    # Split the Structure into components, then re-compose it out of its components
    tmp = structure.split()
    strs, nums = [], []
    for s, n in tmp:
        strs.append(s)
        nums.append(n)
    nums = [ len(n) for n in nums ]

    # Re-compose Structure from components
    new_structure = strs[0]*nums[0]
    for idx in range(1,len(nums)):
        new_structure += strs[idx]*nums[idx]
    # Swap in coordinates again
    new_structure.positions = structure.positions

    # Create System
    newsys = new_structure.createSystem(nonbondedMethod= app.NoCutoff, constraints=None, implicitSolvent=None)

    # Cross check energies
    groups0, groups1, energy0, energy1 = compare_system_energies(pdbfile.topology, pdbfile.topology, system, newsys, pdbfile.positions, verbose = False)

    # Also check that that the number of components is equal to the number I expect
    if not len(nums)==2:
        print("Error: Test system has incorrect number of components.")
        raise Exception('Incorrect number of components in cyclohexane/ethanol test system.')

    # Also check that none of residues have zero charge
    for resnr in range(len(structure.residues)):
        abscharges = [ abs(structure.residues[resnr].atoms[idx].charge) for idx in range(len(structure.residues[resnr].atoms))]
        if sum(abscharges)==0:
            raise Exception('Error: Residue %s in cyclohexane-ethanol test system has a charge of zero, which is incorrect.' % resnr)
