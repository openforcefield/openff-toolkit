# Examples for `openforcefield` tools

## Manifest

* `SMIRNOFF_simulation/` - simulation of a molecule in the gas phase with the SMIRNOFF forcefield format
* `SMIRNOFF_comparison/` - compare molecule energies from SMIRNOFF/OpenMM and AMBER
* `forcefield_modification/` - Jupyter notebook example of modifying a forcefield parameter and evaluating how it changes an energy
* `chemicalEnvironments/` - example and documentation of using chemical environment objects to manipulate environment being considered, generate example SMIRKS, etc. Also contains IPython notebook using the chemical environment for depiction.
* `SMIRNOFF99Frosst` - under-development manual conversion of amber-parm99+parm@Frosst to SMIRNOFF format via an intermediate SMIRNOFFishFrcmod format. The hope is to have an initial guess SMIRNOFF for small molecules.
* `mixedFF_structure` - example of how to use SMIRNOFF forcefields for small molecules in combination with more conventional force fields for proteins and other components of your system.
* `partial_bondorder/` - example of using partial bond orders for parameterizing molecules, and showing how the SMIRNOFF format extends; see https://github.com/open-forcefield-group/smarty/issues/53
* `host_guest_simulation/` - example of setting up a basic MD simulation for a host-guest system using the SMIRNOFF force fields.
