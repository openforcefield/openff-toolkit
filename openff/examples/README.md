# Examples using SMIRNOFF with the toolkit

The following examples are available in [the OpenFF toolkit repository](https://github.com/openforcefield/openff-toolkit/tree/master/examples). Each can be run interactively in the browser with [binder](https://mybinder.org/v2/gh/openforcefield/openforcefield/master?filepath=%2Fexamples%2F), without installing anything on your computer.

## Running the examples

We provide a tool named `openff-examples` to install the examples and their dependencies for you. This can be used to copy the examples to a local, user-writable directory, and manage their dependencies.

List available examples:

```
openff-examples list
```

Copy the `visualization` example to a subfolder of the working directory:

```
openff-examples install visualization
```

Copy the `visualization` example to a subfolder of the working directory and install its dependencies in a new virtual environment:

```
openff-examples install --create-environment visualization
```

Other options, including specifying the target directory, installing dependencies to an existing environment, or specifying the name or location of the new environment, are available with the `--help` switch. When finished, the tool will provide commands to run the examples. `openff-examples` will not be able to install dependencies if it unless `conda` or `mamba` is installed and available in your system `$PATH`.

## Index of provided examples

* [conformer_energies](https://github.com/openforcefield/openff-toolkit/tree/master/examples/conformer_energies) - compute conformer energies of one or more small molecules using a SMIRNOFF force field
* [SMIRNOFF_simulation](https://github.com/openforcefield/openff-toolkit/tree/master/examples/SMIRNOFF_simulation) - simulation of a molecule in the gas phase with the SMIRNOFF force field format
* [forcefield_modification](https://github.com/openforcefield/openff-toolkit/tree/master/examples/forcefield_modification) - modify force field parameters and evaluate how system energy changes
* [using_smirnoff_in_amber_or_gromacs](https://github.com/openforcefield/openff-toolkit/tree/master/examples/using_smirnoff_in_amber_or_gromacs) - convert a System generated with the Open Force Field Toolkit, which can be simulated natively with OpenMM, into AMBER prmtop/inpcrd and GROMACS top/gro input files through the ParmEd library.
* [swap_amber_parameters](https://github.com/openforcefield/openff-toolkit/tree/master/examples/swap_amber_parameters) - take a prepared AMBER protein-ligand system (prmtop and crd) along with a structure file of the ligand, and replace ligand parameters with OpenFF parameters.
* [inspect_assigned_parameters](https://github.com/openforcefield/openff-toolkit/tree/master/examples/inspect_assigned_parameters) - check which parameters are used in which molecules and generate parameter usage statistics.
* [using_smirnoff_with_amber_protein_forcefield](https://github.com/openforcefield/openff-toolkit/tree/master/examples/using_smirnoff_with_amber_protein_forcefield) - use SMIRNOFF parameters for small molecules in combination with more conventional force fields for proteins and other components of your system (using ParmEd to combine parameterized structures)
* [check_dataset_parameter_coverage](https://github.com/openforcefield/openff-toolkit/tree/master/examples/check_dataset_parameter_coverage) - shows how to use the Open Force Field Toolkit to ingest a dataset of molecules, and generate a report summarizing any chemistry that can not be parameterized.
* [visualization](https://github.com/openforcefield/openff-toolkit/tree/master/examples/visualization) - shows how rich representation of `Molecule` objects work in the context of Jupyter Notebooks.
