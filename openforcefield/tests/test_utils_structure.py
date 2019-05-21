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

import os

import parmed
import pytest

from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule, Topology
from openforcefield import utils
import openforcefield.utils.structure as structure


#=============================================================================================
# TESTS
#=============================================================================================

requires_openeye_mol2 = pytest.mark.skipif(not utils.OpenEyeToolkitWrapper.is_available(),
                                               reason='OpenEye is required to parse mol2 files')

class TestUtilsStructure:

    @requires_openeye_mol2
    def test_read_molecules(self):
        molecules = structure.read_molecules('zinc-subset-tripos.mol2.gz', verbose=False)

    @requires_openeye_mol2
    def test_positions(self):
        """Test ability to extract and set positions."""
        molecules = structure.read_molecules('zinc-subset-tripos.mol2.gz', verbose=False)
        positions = structure.extractPositionsFromOEMol(molecules[0])
        structure.setPositionsInOEMol(molecules[0], positions)

    @requires_openeye_mol2
    def test_smirnoff_structure(self):
        """Test function that generates the SMIRNOFF parmed.Structure."""
        molecule = structure.read_molecules('toluene.pdb', verbose=False)[0]
        molecule_structure = structure.generateSMIRNOFFStructure(molecule)
        assert isinstance(molecule_structure, parmed.structure.Structure)

    def test_protein_structure(self):
        from simtk.openmm import app
        pdbfile = utils.get_data_file_path('proteins/T4-protein.pdb')
        proteinpdb = app.PDBFile(pdbfile)
        protein_structure = structure.generateProteinStructure(proteinpdb)
        assert isinstance(protein_structure, parmed.structure.Structure)


@requires_openeye_mol2
def test_merge_system():
    """Test merging of a system created from AMBER and another created from SMIRNOFF."""
    from .utils import create_system_from_amber, get_amber_file_path, get_alkethoh_file_path

    # Create System from AMBER
    prmtop_filename, inpcrd_filename = get_amber_file_path('cyclohexane_ethanol_0.4_0.6')
    system0, topology0, positions0 = create_system_from_amber(prmtop_filename, inpcrd_filename)

    # TODO:
    from openeye import oechem
    # Load simple OEMol
    alkethoh_mol2_filepath = get_alkethoh_file_path('AlkEthOH_c100')[0]
    ifs = oechem.oemolistream(alkethoh_mol2_filepath)
    mol = oechem.OEMol()
    flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
    ifs.SetFlavor(oechem.OEFormat_MOL2, flavor)
    oechem.OEReadMolecule(ifs, mol)
    oechem.OETriposAtomNames(mol)

    # Load forcefield file.
    AlkEthOH_offxml_filename = utils.get_data_file_path('test_forcefields/Frosst_AlkEthOH.offxml')
    forcefield = ForceField(AlkEthOH_offxml_filename)

    # Create OpenMM System and Topology.
    off_mol = Molecule.from_openeye(mol, allow_undefined_stereo=True)
    off_top = Topology.from_molecules([off_mol])
    system1 = forcefield.create_openmm_system(off_top)
    topology1 = structure.generateTopologyFromOEMol(mol)
    positions1 = structure.extractPositionsFromOEMol(mol)

    structure.merge_system(topology0, topology1, system0, system1, positions0, positions1, verbose=True)


@pytest.mark.slow
@requires_openeye_mol2
def test_component_combination():
    """Test that a system still yields the same energy after rebuilding it out of its components
    """
    from simtk import openmm
    from .utils import compare_system_energies, get_packmol_pdb_file_path

    # We've had issues where subsequent instances of a molecule might have zero charges
    # Here we'll try to catch this (and also explicitly check the charges) by re-building
    # a system out of its components

    # Create an OpenMM System from mol2 files containing a cyclohexane-ethanol mixture.
    AlkEthOH_offxml_filename = utils.get_data_file_path('test_forcefields/Frosst_AlkEthOH.offxml')
    forcefield = ForceField(AlkEthOH_offxml_filename)
    pdbfile = openmm.app.PDBFile(get_packmol_pdb_file_path('cyclohexane_ethanol_0.4_0.6'))
    sdf_file_paths = [utils.get_data_file_path(os.path.join('systems', 'monomers', name+'.sdf'))
                      for name in ('ethanol', 'cyclohexane')]
    molecules = [Molecule.from_file(file_path) for file_path in sdf_file_paths]
    topology = Topology.from_openmm(pdbfile.topology, unique_molecules=molecules)
    system = forcefield.create_openmm_system(topology)

    # Convert System to a ParmEd Structure
    structure = parmed.openmm.topsystem.load_topology(topology.to_openmm(), system, pdbfile.positions)

    # Split the Structure into components, then re-compose it out of its components
    tmp = structure.split()
    strs, nums = [], []
    for s, n in tmp:
        strs.append(s)
        nums.append(n)
    nums = [len(n) for n in nums]

    # Re-compose Structure from components
    new_structure = strs[0]*nums[0]
    for idx in range(1,len(nums)):
        new_structure += strs[idx]*nums[idx]
    # Swap in coordinates again
    new_structure.positions = structure.positions

    # Create System
    newsys = new_structure.createSystem(nonbondedMethod=openmm.app.NoCutoff, constraints=None, implicitSolvent=None)

    # Cross check energies
    groups0, groups1, energy0, energy1 = compare_system_energies(pdbfile.topology, pdbfile.topology, system, newsys, pdbfile.positions, verbose = False)

    # Also check that that the number of components is equal to the number I expect
    if not len(nums) == 2:
        print("Error: Test system has incorrect number of components.")
        raise Exception('Incorrect number of components in cyclohexane/ethanol test system.')

    # Also check that none of residues have zero charge
    for resnr in range(len(structure.residues)):
        abscharges = [ abs(structure.residues[resnr].atoms[idx].charge) for idx in range(len(structure.residues[resnr].atoms))]
        if sum(abscharges)==0:
            raise Exception('Error: Residue %s in cyclohexane-ethanol test system has a charge of zero, which is incorrect.' % resnr)
