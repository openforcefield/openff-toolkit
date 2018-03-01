# SMIRNOFF FFXML files for openforcefield.typing.engines.smirnoff.forcefield

## Manifest
- `Frosst_AlkEthOH.offxml`: Original OFFXML file for the AlkEthOH region of chemical space/set as hand created by Christopher Bayly
- `Frosst_AlkEthOH_parmAtFrosst.offxml`: Same as the above but with additional "bugs" introduced (see comments in the FFXML) in order to reproduce problems introduced by atom typing in the AMBER forcefields, typically missing torsions which resulted in specialized torsions instead being given the default. In other words, this forcefield file can reproduce erroneous AMBER energies.
- `Frosst_AlkEthOH_MDL.offxml`: Same as Frosst_AlkEthOH.offxml but uses the MDL aromaticity model for bond perception. No functional difference for this region of chemical space, but illustrates how to invoke that model.
- `Frosst_AlkEthOH_withIDs`: Adds parameter IDs to Frosst_AlkEthOH.offxml for use by the forcefield labeler to show which parameters are used in which places. Otherwise identical.
- `benzene_minimal.offxml`: Minimal offxml file for benzene constructed from smirnoff99Frosst 9/22/16 (dx.doi.org/10.5281/zenodo.154235), intended for use in testing impropers.
- `smirnoff99Frosst.offxml`: SMIRNOFF99Frosst force field as generated in `utilities/convert_frosst` from `smirnoffishFrcmod.parm99Frosst.txt` as per the README.md there. Also maintained as a separate package on GitHub at https://github.com/openforcefield/smirnoff99frosst .
