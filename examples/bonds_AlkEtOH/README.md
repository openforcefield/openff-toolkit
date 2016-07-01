# Example application of SMARTY atom type sampler to recover parm99 typing of alkanes, ethers, and alcohols

In this example, the SMARTY `AtomTypeSampler` is used to attempt to recover SMARTS atom types that recapitulate the typing rules from a referenced set of typed molecules.

## Manifest
* `smarty.py` - example command-line driver
* `atomtypes/` - input atom type sample specification files
* `bondtypes/` - input bond type sample specification files 
* `molecules/` - typed molecule datasets
* `scripts/` - useful conversion scripts

## Usage

Usage

Example:
```
smarty --basetypes=bondtypes/basetype-bonds.smarts --initialtypes=bondtypes/basetype-bonds.smarts --decorators=atomtypes/decorators.smarts --substitutions=bondtypes/substitutions.smarts --molecules=molecules/test_filt1_tripos.mol2 --reference=molecules/test_filt1_ff.mol2 --iterations 500 --temperature=0.1 >|smartyBonds.iter500.log
```

The results below are extracted from the log file smartyBonds.iter500.log.
Initially, the base atom types are added to the pool of current atom types, and the number of atoms and molecules matched by each atom type are shown:
```
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS
    1 :          0          0 |                                                               cc                        [#6]~[#6]
    2 :        518         42 |                                                               ch                        [#6]~[#1]
    3 :        149         42 |                                                               co                        [#6]~[#8]
    4 :        136         42 |                                                               oh                        [#8]~[#1]
    5 :          0          0 |                                                               hh                        [#1]~[#1]
    6 :          0          0 |                                                               oo                        [#8]~[#8]
TOTAL :        803         42
Creating graph matching current atom types with reference atom types...
Graph creation took 0.006 s
Computing maximum weight match...
Maximum weight match took 0.001 s
PROPOSED:
Atom type matches:
cc                                                                       no match
ch                                                               matches       HC :      244 atoms matched
co                                                               matches       CT :      110 atoms matched
oh                                                               matches       OH :       68 atoms matched
hh                                                                       no match
oo                                                                       no match
422 / 803 total atoms match (52.553 %)
Iteration 0 / 500
```

After a number of iterations, the pool of current atom types will have diverged, with a child having been added to the set and unused atom types removed from the original set.
```
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS REF TYPE        FRACTION OF REF TYPED MOLECULES MATCHED
    1 :        518         42 |                                                               ch                        [#6]~[#1]       HC              244 /              244 (100.000%)
    2 :        149         42 |                                                               co                        [#6]~[#8]       CT              110 /              232 ( 47.414%)
    3 :        136         42 |                                                               oh                        [#8]~[#1]       OH               68 /               68 (100.000%)
TOTAL :        803         42 |                                                                                                         422 /      803 match (52.553 %)
```

