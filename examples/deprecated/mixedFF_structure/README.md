## Combining a SMIRNOFF parameterized small molecule with an AMBER parameterized protein using ParmEd

This example illustrates how the [ParmEd](http://parmed.github.io/ParmEd/html/index.html) utility can be used to merge a small molecule parameterized by SMIRNOFF with a traditionally parameterized protein (or other biopolymer) to create a fully parameterized protein-ligand system.

The essential idea is to first create a ParmEd structure from the SMIRNOFF parameterized small molecule
```python
# Load the small molecule
from openforcefield.utils import get_data_file_path
ligand_filename = get_data_file_path('molecules/toluene.mol2')
molecule = Molecule.from_file(ligand_filename)

# Load the smirnoff99Frosst force field
from openforcefield.typing.engines import smirnoff
forcefield = smirnoff.ForceField('test_forcefields/smirnoff99Frosst.offxml')

# Create a ParmEd structure for the molecule
molecule_structure = forcefield.create_parmed_structure(topology=molecule.to_topology(), positions=molecule.positions)
```
Next, a ParmEd structure is created for the protein (parameterized with AMBER using the OpenMM `app.ForceField` implementation):
```python
# Load the protein topology
from openforcefield.utils import generate_protein_structure
protein_pdb_filename = get_data_file_path('proteins/T4-protein.pdb')
protein_pdbfile = app.PDBFile(protein_pdb_filename)

# Load the AMBER protein force field, along with a solvent force field
from simtk.openmm import app
protein_forcefield = 'amber99sbildn.xml'
solvent_forcefield = 'tip3p.xml'
forcefield = app.ForceField(protein_forcefield, solvent_forcefield)

# Parameterize the protein
# TODO: Can we create an OpenMM ffxml file for the ligand?
protein_system = forcefield.createSystem(proteinpdb.topology)

# Create a ParmEd Structure for the protein
protein_structure = parmed.openmm.load_topology(proteinpdb.topology,
                                                protein_system,
                                                xyz=proteinpdb.positions)
```
Finally, the two ParmEd `Structure` objects can be combined by concatenation:
```python
# Combine the ParmEd Structure objects to produce a fully parameterized complex
complex_structure = protein_structure + molecule_structure
```
