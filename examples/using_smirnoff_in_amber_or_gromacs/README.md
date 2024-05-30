## Export OpenFF-generated OpenMM System to AMBER and GROMACS files

This example shows how you can convert an OpenMM `System` generated with the Open Forcefield Toolkit, which can be simulated natively with OpenMM, into AMBER prmtop/inpcrd and GROMACS top/gro input files through the ParmEd library.

#### Manifest:

- `export_with_interchange.ipynyb`: IPython notebook showing how to generate AMBER and GROMACS topology and coordinates files starting from a PDB using Interchange.
- `1_cyclohexane_1_ethanol.pdb`: PDB file containing one molecule of cyclohexane and one of ethanol in vacuum.
