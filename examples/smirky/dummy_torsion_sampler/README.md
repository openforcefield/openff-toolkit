# Torsion sampling with dumb decorators

This example uses the default decorators and compares sampled results
to the reference SMIRFF for AlkEthOH. 

Called from the command line:
```
smirky --molecules AlkEthOH_test_filt1_ff.mol2 --typetag Torsion --verbose \
True --smirff forcefield/Frosst_AlkEtOH.ffxml --iterations 1000
```

Printed to command line:

```
Loading molecules from '/Users/bannanc/anaconda/lib/python2.7/site-packages/smarty/data/molecules/AlkEthOH_test_filt1_ff.mol2'...
42 molecules read
0.007 s elapsed
Attempting to parse file 'odds_files/atom_OR_bases.smarts'
Attempting to parse file 'odds_files/atom_decorators.smarts'
Attempting to parse file 'odds_files/atom_decorators.smarts'
Attempting to parse file 'odds_files/bond_OR_bases.smarts'
Attempting to parse file 'odds_files/bond_AND_decorators.smarts'
Attempting to parse file 'odds_files/atom_index_odds.smarts'
Attempting to parse file 'odds_files/bond_index_odds.smarts'
Creating Torsion sampler took 1.464 s
Iteration 0 / 1000
Iteration 1 / 1000
...
Iteration 998 / 1000
Iteration 999 / 1000
1000 iterations took 532.072 s (0.532 s / iteration)
```

Also created two files
* Torsion_1.00e-1.log - log of all details for each iteration 
* Torsion_1.00e-1.csv - 'trajectory' storing summary information for each iteration

Final score was 1282 / 2175 
Created Torsion                 Matched Reference (number of torsions)
[*:1]~[*:2]~[*:3]~[*:4]         [#1:1]-[#6X4:2]-[#6X4:3]-[#1:4] (574)
[*:1]~;@[*:2]~[*:3]~[*:4]       [#1:1]-[#6X4:2]-[#6X4:3]-[#6X4:4] (552)
[*:1]~;@[$ewg1&A:2]~[*:3]~[*:4] [a,A:1]-[#6X4:2]-[#8X2:3]-[!#1:4] (156)
