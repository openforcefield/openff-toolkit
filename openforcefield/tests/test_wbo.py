# will merge this file with test_toolkits.py later, currently doing it separately


from simtk.openmm.app import PDBFile
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils import get_data_file_path


## Create a simple molecule from SMILES and turn it into a topology.
#mol1 = Molecule.from_smiles("C1=NC=CC(=C1)C2=CC=C(C=C2)O")
#mol1.generate_conformers(n_conformers=1)
#mol1.compute_wiberg_bond_orders()
#
#mol2 = Molecule.from_smiles("C1=NC=CC(=C1)C2=CC=C(C=C2)[O-]")
#mol2.generate_conformers(n_conformers=1)
#mol2.compute_wiberg_bond_orders()

# Reading neutral molecule from file
filename = get_data_file_path('molecules/molecule_neutral.sdf')
molecule1 = Molecule.from_file(filename)  # not using toolkitwrapper flag
molecule1.compute_wiberg_bond_orders()
# Fetching wiberg bond order of of C-C bond
for i in molecule1.atoms[4].bonds:
    if i.atom1_index == 6 or i.atom2_index == 6:
        wbo_neutral = i.fractional_bond_order

# Reading negative molecule from file
filename = get_data_file_path('molecules/molecule_anion.sdf')
molecule2 = Molecule.from_file(filename)  # not using toolkitwrapper flag
molecule2.compute_wiberg_bond_orders()
# Fetching wiberg bond order of of C-C bond
for i in molecule2.atoms[4].bonds:
    if i.atom1_index == 6 or i.atom2_index == 6:
        wbo_anion = i.fractional_bond_order

# Checking that only one additional bond is present ihe neutral molecule
assert (len(molecule1.bonds)==len(molecule2.bonds)+1)
# Wiberg bond order is higher in the anion
assert (wbo_anion > wbo_neutral)



#from IPython import embed
#embed();
#exit()
