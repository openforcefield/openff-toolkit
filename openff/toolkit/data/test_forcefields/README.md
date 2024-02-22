# SMIRNOFF FFXML files for openff.toolkit.typing.engines.smirnoff.forcefield

## Manifest
- `Frosst_AlkEthOH.offxml`: Original OFFXML file for the AlkEthOH region of chemical space/set as hand created by Christopher Bayly
- `Frosst_AlkEthOH_parmAtFrosst.offxml`: Same as the above but with additional "bugs" introduced (see comments in the FFXML) in order to reproduce problems introduced by atom typing in the AMBER force fields, typically missing torsions which resulted in specialized torsions instead being given the default. In other words, this force field file can reproduce erroneous AMBER energies.
- `Frosst_AlkEthOH_MDL.offxml`: Same as Frosst_AlkEthOH.offxml but uses the MDL aromaticity model for bond perception. No functional difference for this region of chemical space, but illustrates how to invoke that model.
- `Frosst_AlkEthOH_withIDs`: Adds parameter IDs to Frosst_AlkEthOH.offxml for use by the force field labeler to show which parameters are used in which places. Otherwise identical.
- `benzene_minimal.offxml`: Minimal offxml file for benzene constructed from smirnoff99Frosst 9/22/16 (dx.doi.org/10.5281/zenodo.154235), intended for use in testing impropers.
- `smirnoff99Frosst.offxml`: (DEPRECATED) SMIRNOFF99Frosst force field as generated in `utilities/convert_frosst` from `smirnoffishFrcmod.parm99Frosst.txt` as per the README.md there. 
  Also maintained as a separate package on GitHub at https://github.com/openforcefield/smirnoff99frosst. 
  **NOTE:** There are many versions of the smirnoff99Frosst force field, and this file **SHOULD NOT** be assumed to correspond to **ANY** released version.
  This file exists only as a basis for automated testing and may be updated at any time.
  To find properly-versioned smirnoff99Frosst files suitable for use in publication, please see https://github.com/openforcefield/smirnoff99Frosst/
- `smirnoff99frosst_experimental.offxml`: Temporarily added to address hierarchy problems, redundancies, SMIRKS pattern typos etc., as documented in issue #367. 
  Will ultimately be propagated to an updated force field in the openforcefield/smirnoff99frosst repo.`
- `smirff99Frosst_reference_0_1_spec.offxml`: a SMIRNOFF 0.1 spec file enclosed by the legacy SMIRFF tag. This file is used in backwards-compatibility testing.
- `smirnoff99Frosst_reference_0_1_spec.offxml`: a SMIRNOFF 0.1 spec file enclosed by the legacy SMIRFF tag. This file is used in backwards-compatibility testing.
- `smirnoff99Frosst_reference_0_2_spec.offxml`: a SMIRNOFF 0.2 spec file enclosed by the legacy SMIRFF tag. This file is used in backwards-compatibility testing.
- `smirnoff99Frosst_reference_0_1_spec.offxml`: a SMIRNOFF 0.1 spec file enclosed by the legacy SMIRFF tag. This file is used in backwards-compatibility testing.
- `GBSA_HCT-1.0.offxml`: An experimental force field used in validation tests against OpenMM's implementation of the [Hawkins-Cramer-Truhlar](http://docs.openmm.org/latest/userguide/zbibliography.html#hawkins1995) GBSA model (corresponding to igb=1 in AMBER) 
- `GBSA_OBC1-1.0.offxml`: An experimental force field used in validation tests against OpenMM's implementation of the [Onufriev-Bashford-Case](http://docs.openmm.org/latest/userguide/zbibliography.html#onufriev2004) using the GB(OBC)I model (corresponding to igb=2 in AMBER) 
- `GBSA_OBC2-1.0.offxml`: An experimental force field used in validation tests against OpenMM's implementation of the [Onufriev-Bashford-Case](http://docs.openmm.org/latest/userguide/zbibliography.html#onufriev2004) using the GB(OBC)II model (corresponding to igb=5 in AMBER) 
