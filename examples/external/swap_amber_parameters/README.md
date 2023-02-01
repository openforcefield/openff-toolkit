## Taking an AMBER system and replacing the ligand parameters with OpenFF parameters

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/openforcefield/openff-toolkit/stable?filepath=examples%2Fswap_amber_parameters%2F)

These examples illustrate how the [ParmEd](http://parmed.github.io/ParmEd/html/index.html) utility can be used to take a fully parameterized (and solvated) protein-ligand system and replace the ligand parameters with those from OpenFF. In this case, the system uses an AMBER protein force field, ions, and TIP3P water, and we replace the GAFF ligand parameters with OpenFF 1.0 parameters, writing out the resulting fully parameterized system.

### BRD4:inhibitor complex

* [`swap_existing_ligand_parameters.ipynb`](swap_existing_ligand_parameters.ipynb) contains an example illustrating taking a fully parameterized BRD4 protein-ligand system, with an AMBER protein force field and GAFF ligand parameters, and replacing the ligand parameters with OpenFF parameters from SMIRNOFF format. The BRD4:inhibitor complex is taken from the [free energy benchmark systems living review](https://www.annualreviews.org/doi/abs/10.1146/annurev-biophys-070816-033654) [GitHub repo](https://github.com/MobleyLab/benchmarksets/tree/master/input_files/BRD4).

* [`swap_existing_ligand_parameters_with_openmmforcefields.ipynb`](swap_existing_ligand_parameters_with_openmmforcefields.ipynb) contains an example that uses the [`openmmforcefields`](http://github.com/openmm/openmmforcefields) package to take a fully parameterized BRD4 protein-ligand system, with an AMBER protein force field and GAFF ligand parameters, and replacing the ligand parameters with OpenFF parameters from SMIRNOFF format. The BRD4:inhibitor complex is taken from the [free energy benchmark systems living review](https://www.annualreviews.org/doi/abs/10.1146/annurev-biophys-070816-033654) [GitHub repo](https://github.com/MobleyLab/benchmarksets/tree/master/input_files/BRD4).
