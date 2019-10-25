# will merge this file with test_toolkits.py later, currently doing it separately


from simtk.openmm.app import PDBFile
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField

def test_compute_partial_charges(self):
    """Test computation/retrieval of partial charges"""
    # TODO: Test only one molecule for speed?
    # TODO: Do we need to deepcopy each molecule, or is setUp called separately for each test method?
    from simtk import unit
    import numpy as np

    return


# Create a simple molecule from SMILES and turn it into a topology.
mol1 = Molecule.from_smiles("C1=NC=CC(=C1)C2=CC=C(C=C2)O")
mol1.generate_conformers(n_conformers=1)
mol1.compute_wiberg_bond_orders()
for i in mol1.bonds: 
   print(i.fractional_bond_order)
print("Total number of bonds:", len(mol1.bonds))

mol2 = Molecule.from_smiles("C1=NC=CC(=C1)C2=CC=C(C=C2)[O-]")
mol2.generate_conformers(n_conformers=1)
mol2.compute_wiberg_bond_orders()
for i in mol2.bonds: 
   print(i.fractional_bond_order)
print("Total number of bonds:", len(mol1.bonds))


#from IPython import embed
#embed();
#exit()
