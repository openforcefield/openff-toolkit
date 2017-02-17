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
    --iterations 500 --temperature=0 >| smartyAlkEtOH.iter500.log
```

At T=0, this has roughly a 50% success rate (it's simply a greedy algorithm in that case); at T=0.01, a long enough run will typically achieve 100%, and at T=0.1, success is relatively rare.

The results below are extracted from the log file smartyAlkEtOH.iter500.log from T=0.
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
    1 :         68         42 |                                                       c_hydrogen                             [#1]       HO               68 /               68 (100.000%)
    2 :        232         42 |                                                         c_carbon                             [#6]       CT              232 /              232 (100.000%)
    3 :         39         30 |                                                         c_oxygen                             [#8]       OS               39 /               39 (100.000%)
    4 :         68         42 |                                       c_oxygen hydrogen-adjacent                    [#8$(*~[#1])]       OH               68 /               68 (100.000%)
    5 :        244         38 |                                      c_hydrogen simply c_carbon                     [#1$(*-[#6])]       HC              244 /              244 (100.000%)
    6 :        125         37 |                       c_hydrogen simply c_carbon simply c_oxygen               [#1$(*-[#6]-[#8])]       H1              116 /              116 (100.000%)
    7 :         27         21 | c_hydrogen simply c_carbon simply c_oxygen (any c_carbon) (simply c_oxygen) [#1$(*-[#6](-[#8])(~[#6])-[#8])]       H2               27 /               33 ( 81.818%)
TOTAL :        803         42 |                                                                                                         794 /      803 match (98.879 %)
```

In this case, by the end of the run, SMARTY has found and maintained SMARTS patterns which match all of the AlkEthOH atom types, giving it a 100% success rate:

```
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS REF TYPE        FRACTION OF REF TYPED MOLECULES MATCHED
    1 :         68         42 |                                                       c_hydrogen                             [#1]       HO               68 /               68 (100.000%)
    2 :        232         42 |                                                         c_carbon                             [#6]       CT              232 /              232 (100.000%)
    3 :         39         30 |                                                         c_oxygen                             [#8]       OS               39 /               39 (100.000%)
    4 :         68         42 |                                       c_oxygen hydrogen-adjacent                    [#8$(*~[#1])]       OH               68 /               68 (100.000%)
    5 :        244         38 |                                      c_hydrogen simply c_carbon                     [#1$(*-[#6])]       HC              244 /              244 (100.000%)
    6 :        116         37 |                       c_hydrogen simply c_carbon simply c_oxygen               [#1$(*-[#6]-[#8])]       H1              116 /              116 (100.000%)
    7 :         33         22 |     c_hydrogen simply c_carbon simply c_oxygen (simply c_oxygen)        [#1$(*-[#6](-[#8])-[#8])]       H2               33 /               33 (100.000%)
    8 :          3          3 | c_hydrogen simply c_carbon simply c_oxygen (simply c_oxygen) (simply c_oxygen) [#1$(*-[#6](-[#8])(-[#8])-[#8])]       H3                3 /                3 (100.000%)
TOTAL :        803         42 |                                                                                                         803 /      803 match (100.000 %)
```

The end of the log file concludes with printing a hierarchy of the atom types which smarty discovered:
```
Atom type hierarchy:
    [#1]
        [#1$(*-[#6])]
            [#1$(*-[#6]-[#8])]
                [#1$(*-[#6](-[#8])-[#8])]
                    [#1$(*-[#6](-[#8])(-[#8])-[#8])]
    [#8]
        [#8$(*~[#1])]
    [#6]
```

This is in keeping with what is expected for AlkEthOH -- hydrogen has OH, HC, H1, H2, and H3, where HC, H1, H2, and H3 differ by the number of electron withdrawing groups attached to the carbon adjacent to the hydrogen.
