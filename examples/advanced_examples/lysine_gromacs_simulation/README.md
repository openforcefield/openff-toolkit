## Create and simulate a SMIRNOFF-parametrized system in GROMACS

This example shows how to use the Open Force Field Toolkit to produce an OpenMM system for a small organic molecule (lysine), convert it to a pair of GROMACS input files (`top`/`gro`), add waters and ions, and run a GROMACS simulation. 

#### Manifest:

- `Lysine2Gromacs.ipynb`: IPython notebook showing how to generate GROMACS topology and coordinates files starting from a PDB, and simulate them using GROMACS.
- `LYS_IDEAL.pdb`: PDB file containing a single conformer of lysine.
