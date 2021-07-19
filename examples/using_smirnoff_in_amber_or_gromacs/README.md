## Export OpenFF-generated OpenMM System to AMBER and GROMACS files

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/openforcefield/openff-toolkit/latest?filepath=examples%2Fusing_smirnoff_in_amber_or_gromacs%2Fconvert_to_amber_gromacs.ipynb)

This example shows how you can convert an OpenMM `System` generated with the Open Forcefield Toolkit, which can be simulated natively with OpenMM, into AMBER prmtop/inpcrd and GROMACS top/gro input files through the ParmEd library.

#### Manifest:

- `convert_to_amber_gromacs.ipynb`: IPython notebook showing how to generate AMBER and GROMACS topology and coordinates files starting from a PDB.
- `1_cyclohexane_1_ethanol.pdb`: PDB file containing one molecule of cyclohexane and one of ethanol in vacuum.
