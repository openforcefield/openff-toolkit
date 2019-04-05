# Examples using SMIRNOFF with the toolkit

The following examples are available in [the openforcefield toolkit repository](https://github.com/openforcefield/openforcefield/tree/master/examples):

### Index of provided examples

* `SMIRNOFF_simulation/` - simulation of a molecule in the gas phase with the SMIRNOFF forcefield format
* `forcefield_modification/` - modify forcefield parameters and evaluate how system energy changes
* `conversion_amber_gromacs` - convert a System generated with the Open Forcefield Toolkit, which can be simulated natively with OpenMM, into AMBER prmtop/inpcrd and GROMACS top/gro input files through the ParmEd library.`
* `label_molecule/` - check which parameters are used in which molecules, do parameter usage statistics, etc.
* `mixedFF_structure/` - use SMIRNOFF parameters for small molecules in combination with more conventional force fields for proteins and other components of your system (using ParmEd to combine parameterized structures)
