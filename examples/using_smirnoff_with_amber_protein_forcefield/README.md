## Combining a SMIRNOFF parameterized small molecule with an Amber parameterized protein using ParmEd

These examples illustrate how the [ParmEd](http://parmed.github.io/ParmEd/html/index.html) utility can be used to merge a small molecule parameterized by SMIRNOFF with a traditionally parameterized protein (or other biopolymer) to create a fully parameterized protein-ligand system.

### BRD4:inhibitor complex

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/openforcefield/openforcefield/master?filepath=examples%2Fusing_smirnoff_with_amber_protein_forcefield%2FBRD4_inhibitor_benchmark.ipynb)

[`BRD4_inhibitor_benchmark.ipynb`](BRD4_inhibitor_benchmark.ipynb) contains an example illustrating applying SMIRNOFF and Amber14 parameters to a BRD4:inhibitor complex taken from the [free energy benchmark systems living review](https://www.annualreviews.org/doi/abs/10.1146/annurev-biophys-070816-033654) [GitHub repo](https://github.com/MobleyLab/benchmarksets/tree/master/input_files/BRD4).

### Toluene in complex with T4 lysozyme L99A in TIP3P-FB water

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/openforcefield/openforcefield/master?filepath=examples%2Fusing_smirnoff_with_amber_protein_forcefield%2Ftoluene_in_T4_lysozyme.ipynb)

[`toluene_in_T4_lysozyme.ipynb`](toluene_in_T4_lysozyme.ipynb) contains an example illustrating applying SMIRNOFF and Amber99SB-ILDN parameters to toluene complexed with T4 lysozyme L99A in a TIP3P-FB water box.
