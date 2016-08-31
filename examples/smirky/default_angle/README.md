### SMIRKY Basic Example

This is a demonstration of the default behavior of smirky
that is only the required information was given. 
A set of molecules to type and a typetag. 
Files for the decorators and odds_files being used are available at
smart/data/odds_files. 

From command line called:
```
smirky -molecules AlkEthOH_test_filt1_ff.mol2 --typetag Angle --verbose True
```

Output to command line:
```
Loading molecules from '/Users/bannanc/anaconda/lib/python2.7/site-packages/smarty/data/molecules/AlkEthOH_test_filt1_ff.mol2'...
42 molecules read
0.009 s elapsed
Attempting to parse file 'odds_files/atom_OR_bases.smarts'
Attempting to parse file 'odds_files/atom_decorators.smarts'
Attempting to parse file 'odds_files/atom_decorators.smarts'
Attempting to parse file 'odds_files/bond_OR_bases.smarts'
Attempting to parse file 'odds_files/bond_AND_decorators.smarts'
Attempting to parse file 'odds_files/atom_index_odds.smarts'
Attempting to parse file 'odds_files/bond_index_odds.smarts'
Creating Angle sampler took 0.315 s
Iteration 0 / 150
Iteration 1 / 150
Iteration 2 / 150
...
Iteration 148 / 150
Iteration 149 / 150
150 iterations took 37.055 s (0.247 s / iteration)
```

Also created output files
* Angle_1.00e-1.csv - 'trajectory' file with summary of changes at each iteration
* Angle_1.00e-1.log - detailed information with each step
