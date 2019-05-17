#!/bin/env python

"""
Illustrate how to combine a SMIRNOFF parameterized small molecule with an AMBER parameterized protein using ParmEd.

"""

#
# Load and parameterize the small molecule
#

# Load the small molecule
from openforcefield.utils import get_data_file_path
ligand_filename = get_data_file_path('molecules/toluene.mol2')
molecule = Molecule.from_file(ligand_filename)

# Load the smirnoff99Frosst force field
from openforcefield.typing.engines import smirnoff
forcefield = smirnoff.ForceField('test_forcefields/smirnoff99Frosst.offxml')

# Create a ParmEd structure for the molecule
molecule_structure = forcefield.create_parmed_structure(topology=molecule.to_topology(), positions=molecule.positions)
print('Molecule:', molecule_structure)

#
# Load and parameterize the protein
#

# Load the protein topology
protein_pdb_filename = get_data_file_path('proteins/T4-protein.pdb')
protein_pdbfile = app.PDBFile(protein_pdb_filename)

# Load the AMBER protein force field, along with a solvent force field
from simtk.openmm import app
protein_forcefield = 'amber99sbildn.xml'
solvent_forcefield = 'tip3p.xml'
forcefield = app.ForceField(protein_forcefield, solvent_forcefield)

# Parameterize the protein
protein_system = forcefield.createSystem(proteinpdb.topology)

# Create a ParmEd Structure for the protein
protein_structure = parmed.openmm.load_topology(proteinpdb.topology,
                                                protein_system,
                                                xyz=proteinpdb.positions)
print('Protein:', protein_structure)

# Combine the ParmEd Structure objects to produce a fully parameterized complex
# that can now be exported to AMBER, CHARMM, OpenMM, and other packages
# Note that the addition should always add the small molecule second so the box vectors if the first item (the protein) are to be preserved
complex_structure = protein_structure + molecule_structure
print('Complex:', complex_structure)

# TODO: How can we add solvent while ensuring the ligand doesn't overlap with solvent molecules?
# TODO: Can we have SMIRNOFF ForceField create an OpenMM ffxml file for the ligand, and then use the OpenMM pipeline?
# TODO: Or can OpenMM just use dummy parameters?
