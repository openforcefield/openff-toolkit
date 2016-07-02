# Example application of SMARTY atom type sampler to recover parm99 typing of alkanes, ethers, and alcohols

In this example, the SMARTY `AtomTypeSampler` is used to attempt to recover SMARTS atom types that recapitulate the typing rules from a referenced set of typed molecules.

## Manifest
* `smarty.py` - example command-line driver
* `atomtypes/` - input atom type sample specification files
* `molecules/` - typed molecule datasets
* `scripts/` - useful conversion scripts

## Usage

Usage

Example:
```
smarty --basetypes=atomtypes/basetypes-elemental.smarts --initialtypes=atomtypes/basetypes-elemental.smarts --decorators=atomtypes/decorators.smarts --substitutions=atomtypes/substitutions.smarts \
    --molecules=molecules/test_filt1_tripos.mol2 --reference=molecules/test_filt1_ff.mol2 \
    --iterations 500 --temperature=0.1 >| smartyAlkEtOH.iter500.log
```

The results below are extracted from the log file smartyAlkEtOH.iter500.log.
Initially, the base atom types are added to the pool of current atom types, and the number of atoms and molecules matched by each atom type are shown:
```
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS
    1 :        464         42 |                                                       c_hydrogen                             [#1]
    2 :        232         42 |                                                         c_carbon                             [#6]
    3 :          0          0 |                                                       c_nitrogen                             [#7]
    4 :        107         42 |                                                         c_oxygen                             [#8]
    5 :          0          0 |                                                       c_fluorine                             [#9]
    6 :          0          0 |                                                    c_phosphorous                            [#15]
    7 :          0          0 |                                                         c_sulfur                            [#16]
    8 :          0          0 |                                                       c_chlorine                            [#17]
    9 :          0          0 |                                                        c_bromine                            [#35]
   10 :          0          0 |                                                         c_iodine                            [#53]
TOTAL :        803         42
```
After a number of iterations, the pool of current atom types will have diverged, with a child having been added to the set and unused atom types removed from the original set.
```
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS REF TYPE        FRACTION OF REF TYPED MOLECULES MATCHED
    1 :        464         42 |                                                       c_hydrogen                             [#1]       HC              244 /              244 (100.000%)
    2 :          0          0 |                                                         c_carbon                             [#6]
    3 :        232         42 |                                                 c_carbon neutral                          [#6&+0]       CT              232 /              232 (100.000%)
    4 :        107         42 |                                                         c_oxygen                             [#8]       OH               68 /               68 (100.000%)
TOTAL :        803         42 |                                                                                                         544 /      803 match (67.746 %)
```
