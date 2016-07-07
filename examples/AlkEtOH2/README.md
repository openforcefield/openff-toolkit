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
smarty --basetypes=atomtypes/basetypes-elemental.smarts --decorators=atomtypes/decorators.smarts --substitutions=atomtypes/substitutions.smarts \
    --molecules=molecules/test_filt1_tripos.mol2 --reference=molecules/test_filt1_ff.mol2 \
    --iterations 500 --temperature=0.1 >| smartyAlkEtOH.iter500.log
```

The results below are extracted from the log file smartyAlkEtOH.iter500.log.
Initially, the base atom types are added to the pool of current atom types, and the number of atoms and molecules matched by each atom type are shown:
```
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS
    1 :        464         42 |                                                         hydrogen                             [#1]
    2 :        232         42 |                                                           carbon                             [#6]
    3 :          0          0 |                                                         nitrogen                             [#7]
    4 :        107         42 |                                                           oxygen                             [#8]
    5 :          0          0 |                                                         fluorine                             [#9]
    6 :          0          0 |                                                      phosphorous                            [#15]
    7 :          0          0 |                                                           sulfur                            [#16]
    8 :          0          0 |                                                         chlorine                            [#17]
    9 :          0          0 |                                                          bromine                            [#35]
   10 :          0          0 |                                                           iodine                            [#53]
TOTAL :        803         42
```
After a number of iterations, the pool of current atom types will have diverged, with a child having been added to the set and unused atom types removed from the original set.
```
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS REF TYPE        FRACTION OF REF TYPED MOLECULES MATCHED
    1 :        464         42 |                                                         hydrogen                             [#1]       HC              244 /              244 (100.000%)
    2 :        232         42 |                                                           carbon                             [#6]       CT              232 /              232 (100.000%)
    3 :         39         30 |                                                           oxygen                             [#8]       OS               39 /               39 (100.000%)
    4 :         68         42 |                                           oxygen total-h-count-1                          [#8&H1]       OH               68 /               68 (100.000%)
TOTAL :        803         42 |                                                                                                         583 /      803 match (72.603 %)
```
