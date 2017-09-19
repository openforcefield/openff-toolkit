# Setting up a host-guest simulation using SMIRNOFF parameters

This attempts to provide an example of how easy it is to set up rather nontrivial simulations with SMIRNOFF force fields.
In this case we will begin with a guest SMILES string and a 3D structure of a host, dock the host into the guest, solvate, and proceed to simulations.

## Authorship

Example provided by David Mobley (UCI)

## Source materials
- Based in part on a host-guest [docking Jupyter notebook](https://github.com/MobleyLab/SAMPL6/blob/master/host_guest/GenerateInputs.ipynb) I constructed for SAMPL6


## Manifest
- smirnoff_host_guest.ipynb
- hbonds.ffxml: Force field XML file for constraining hbonds (from https://github.com/MobleyLab/SMIRNOFF_paper_code/blob/master/scripts/hbonds.ffxml, courtesy Andrea Rizzi )
- OA.mol2: mol2 file of Octa Acid with AM1-BCC charges, as distributed with SAMPL5 (where it was labeled OAH.mol2). Atom types converted from GAFF to SYBYL via Antechember. 
