# Previously released SMIRNOFF force field

These files are mostly used for regression tests. They have been modified by hand and do not necessarily represent the parameters in any version of the smirnoff99Frosst force field. Each FF is stored in three versions of the SMIRNOFF spec, and tests are run with all versions to ensure equivalence of SMIRKS-based parameter assignment.

## Manifest
- `test_ff_0_0_2_spec_0_X.offxml`: Derived from a version of the smirnoff99Frosst force field released with `openforcefield 0.0.2`.
- `test_ff_0_0_4_spec_0_X.offxml`: Derived from a version of the smirnoff99Frosst force field released with `openforcefield 0.0.4`.
- `test_ff_0_0_4_fixed_spec_0_X.offxml`: Derived from a version of the smirnoff99Frosst force field released with `openforcefield 0.0.4` with a bugfix correction for the issue described in #179.
- `hbonds.offxml`: Separate file adding only hydrogen bond constraints to the parameterization.
