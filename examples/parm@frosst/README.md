# Example application of SMARTY atom type sampler to recover parm@frosst typing

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
    --molecules=molecules/zinc-subset-tripos.mol2.gz --reference=molecules/zinc-subset-parm@frosst.mol2.gz --iterations 1000 --temperature=0.1
```

Initially, the base atom types are added to the pool of current atom types, and the number of atoms and molecules matched by each atom type are shown:
```
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS REF TYPE        FRACTION OF REF TYPED MOLECULES MATCHED
    1 :      88148       7487 |                                                         hydrogen                             [#1]       HA            28720 /            28720 (100.000%)
    2 :      90146       7505 |                                                           carbon                             [#6]       CA            37143 /            37143 (100.000%)
    3 :      20838       6806 |                                                         nitrogen                             [#7]       NB             7612 /             7612 (100.000%)
    4 :      12829       5946 |                                                           oxygen                             [#8]        O             4876 /             4876 (100.000%)
    5 :       1001        444 |                                                         fluorine                             [#9]        F             1001 /             1001 (100.000%)
    6 :          5          5 |                                                      phosphorous                            [#15]        P                5 /                5 (100.000%)
    7 :       3171       2593 |                                                           sulfur                            [#16]        S             2544 /             2544 (100.000%)
    8 :        574        463 |                                                         chlorine                            [#17]       CL              574 /              574 (100.000%)
    9 :         84         73 |                                                          bromine                            [#35]       BR               84 /               84 (100.000%)
   10 :          8          8 |                                                           iodine                            [#53]        I                8 /                8 (100.000%)
TOTAL :     216804       7505 |                                                                                                       82567 /   216804 match (38.084 %)
```
After a few iterations, the pool of current atom types will have diverged, with some children having been added to the set or atom types removed from the original set.
```

```
