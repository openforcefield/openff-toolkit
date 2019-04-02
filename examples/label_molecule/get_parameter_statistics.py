#!/bin/env python

from openforcefield.utils import get_data_filename
#from openforcefield.typing.engines.smirnoff import get_molecule_parameterIDs
from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils.structure import get_molecule_parameterIDs
#
# def get_molecule_parameterIDs(molecules, offxml_filename):
#     # Create storage
#     parameters_by_molecule = {}
#     parameters_by_ID = {}
#
#     # Generate isomeric SMILES
#     isosmiles = list()
#     for mol in molecules:
#         smi = mol.to_smiles()
#         if not smi in isosmiles:
#             isosmiles.append(smi)
#         # If the molecule is already here, raise exception
#         else:
#             raise ValueError(
#                 "Error: get_molecule_parameterIDs has been provided a list of oemols which contains the same molecule, having isomeric smiles %s, more than once." % smi)
#     # Label molecules
#     ff = ForceField(ffxml)
#     labels = ff.labelMolecules(oemols)
#
#     # Organize labels into output dictionary by looping over all molecules/smiles
#     for idx in range(len(isosmiles)):
#         # Pull smiles, initialize storage
#         smi = isosmiles[idx]
#         parameters_by_molecule[smi] = []
#
#         # Organize data for this molecule
#         data = labels[idx]
#         for force_type in data.keys():
#             for (atom_indices, pid, smirks) in data[force_type]:
#                 # Store pid to molecule
#                 parameters_by_molecule[smi].append(pid)
#
#                 # Store which molecule this pid occurred in
#                 if pid not in parameters_by_ID:
#                     parameters_by_ID[pid] = set()
#                     parameters_by_ID[pid].add(smi)
#                 else:
#                     parameters_by_ID[pid].add(smi)
#
#     return parameters_by_molecule, parameters_by_ID

#oemols = read_molecules(get_data_filename('molecules/AlkEthOH_test_filt1_tripos.mol2'))
mols = Molecule.from_file(get_data_filename('molecules/AlkEthOH_test_filt1_tripos.mol2'))
offxml_filename = ForceField(get_data_filename('forcefield/Frosst_AlkEthOH.offxml'))
#ffxml = get_data_filename('forcefield/Frosst_AlkEthOH.offxml')



parameters_by_molecule, parameters_by_ID = get_molecule_parameterIDs( mols, offxml_filename)

# Print some info
print("Parameters by molecule:")
for smi in parameters_by_molecule.keys():
    print(smi, parameters_by_molecule[smi])


print("Molecules with parameter IDs:")
for pid in parameters_by_ID.keys():
    print(pid, parameters_by_ID[pid])
