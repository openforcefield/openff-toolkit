## Setting up a host-guest simulation using SMIRNOFF parameters

This example illustrates how easy it is to set up rather nontrivial simulations with SMIRNOFF force fields.
In this case, we will begin with a guest SMILES string and a 3D structure of a host, dock the host into the guest, solvate, and proceed to simulations.

### Authorship

This example was provided by David Mobley (UCI)

### Source materials
- Based in part on a host-guest [docking Jupyter notebook](https://github.com/MobleyLab/SAMPL6/blob/master/host_guest/GenerateInputs.ipynb) I constructed for SAMPL6

### Example files
* [smirnoff_host_guest.ipynb](https://github.com/openforcefield/openforcefield/blob/master/examples/host_guest_simulation/smirnoff_host_guest.ipynb): Jupyter notebook illustrating setup of the host-guest system
* [hbonds.offxml](https://github.com/openforcefield/openforcefield/blob/master/examples/host_guest_simulation/hbonds.ffxml): Force field XML file for constraining hbonds (from https://github.com/MobleyLab/SMIRNOFF_paper_code/blob/master/scripts/hbonds.offxml, courtesy Andrea Rizzi)
* `OA.mol2`: mol2 file of Octa Acid with AM1-BCC charges, as distributed with SAMPL5 (where it was labeled `OAH.mol2`). Atom types converted from GAFF to SYBYL via Antechember.
