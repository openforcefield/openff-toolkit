# Example application of SMARTY atom type sampler to recover parm@frosst typing

In this example, the SMARTY `AtomTypeSampler` is used to attempt to recover SMARTS atom types that recapitulate the typing rules from a referenced set of typed molecules.

## Manifest
* `smarty.py` - example command-line driver
* `atomtypes/` - input atom type sample specification files
* `molecules/` - typed molecule datasets
* `convert-atom-names-to-tripos.py` - utility to convert atom names to Tripos in mol2 files

## Usage

Usage

Example:
```
python sample.py --basetypes=atomtypes/basetypes-elemental.smarts --decorators=atomtypes/decorators.smarts --substitutions=atomtypes/substitutions.smarts \
    --molecules=molecules/zinc-subset-tripos.mol2.gz --reference=molecules/zinc-subset-parm@frosst.mol2.gz --iterations 1000 --temperature=0.1
```

Initially, the base atom types are added to the pool of current atom types, and the number of atoms and molecules matched by each atom type are shown:
```
INDEX        ATOMS  MOLECULES                                          TYPE NAME                                           SMARTS
    1 :       4371        483 |                                         hydrogen                                             [#1]
    2 :       2717        488 |                                           carbon                                             [#6]
    3 :        117        102 |                                         nitrogen                                             [#7]
    4 :        338        216 |                                           oxygen                                             [#8]
    5 :         69         24 |                                         fluorine                                             [#9]
    6 :          2          2 |                                      phosphorous                                            [#15]
    7 :         20         18 |                                           sulfur                                            [#16]
    8 :        116         62 |                                         chlorine                                            [#17]
    9 :         28         23 |                                          bromine                                            [#35]
   10 :         12         11 |                                           iodine                                            [#53]
TOTAL         7790        490
```
After many iterations, the pool of current atom types will have diverged, with some children having been added to the set.  (Base atom types can never be deleted.)
```
Iteration 999 / 1000: True
INDEX        ATOMS  MOLECULES                                          TYPE NAME                                           SMARTS
    1 :       4366        483 |                                         hydrogen                                             [#1]
    2 :          5          5 |                         hydrogen sulfur-adjacent                                  [#1&$(*~[#16])]
    3 :       2602        474 |                                           carbon                                             [#6]
    4 :         22         18 |                         carbon fluorine-adjacent                                   [#6&$(*~[#9])]
    5 :          7          6 |       carbon fluorine-adjacent hydrogen-adjacent                         [#6&$(*~[#9])&$(*~[#1])]
    6 :         25         23 |                          carbon bromine-adjacent                                  [#6&$(*~[#35])]
    7 :         61         33 |         carbon total-h-count-1 nitrogen-adjacent                                [#6&H1&$(*~[#7])]
    8 :        105         92 |                                         nitrogen                                             [#7]
    9 :         12         12 |                           nitrogen triple-bonded                                  [#7&$([*]#[*])]
   10 :        338        216 |                                           oxygen                                             [#8]
   11 :         69         24 |                                         fluorine                                             [#9]
   12 :          2          2 |                                      phosphorous                                            [#15]
   13 :         16         14 |                                           sulfur                                            [#16]
   14 :          4          4 |                                 sulfur valence-6                                         [#16&v6]
   15 :        116         62 |                                         chlorine                                            [#17]
   16 :         28         23 |                                          bromine                                            [#35]
   17 :         12         11 |                                           iodine                                            [#53]
TOTAL         7790        490
```
