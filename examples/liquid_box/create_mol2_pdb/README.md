## Contact:
Lee-Ping Wang leeping@ucdavis.edu

## Description:
Provide a text file containing a single SMILES string and three-letter residue name, as:
`python create_mol2_pdb.py input.smi ETH`
Receive `ETH.pdb` and `ETH.mol2` files containing a single molecule with 1 conformation.
Receive `ETH-box.pdb` containing a box with specified number of molecules.

## Dependencies:
- OpenEye tools (for creating molecule from SMILES and conformer generation)
- openmoltools (OpenEye conformer generation wrapper)
- Gromacs 4.6.7 (for calling genbox to create solvent box)
- ForceBalance 1.5.x (for putting information back that was thrown away by genbox)
