## Taking an AMBER system and replacing the ligand parameters with OpenFF parameters

These examples illustrate how the [OpenMM Force Fields](http://github.com/openmm/openmmforcefields) utility can be used to parametrize a system by combining Amber and OpenFF force fields.

### BRD4:inhibitor complex


* [`swap_existing_ligand_parameters_with_openmmforcefields.ipynb`](swap_existing_ligand_parameters_with_openmmforcefields.ipynb) contains an example that uses the [`openmmforcefields`](http://github.com/openmm/openmmforcefields) package to take a fully parameterized BRD4 protein-ligand system, with an AMBER protein force field and GAFF ligand parameters, and replace the ligand parameters with OpenFF parameters from SMIRNOFF format. The BRD4:inhibitor complex is taken from the [free energy benchmark systems living review](https://www.annualreviews.org/doi/abs/10.1146/annurev-biophys-070816-033654) [GitHub repo](https://github.com/MobleyLab/benchmarksets/tree/master/input_files/BRD4).
