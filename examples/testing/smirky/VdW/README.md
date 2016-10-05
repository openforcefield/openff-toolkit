# Example of sampling VdW (atomtypes)

This example uses the the decorator files in this directory 
instead of the defaults in smarty/data/odds_files. 

Also, the script run_smirky.bash was created so that repeated experiments 
with the odds in these decorator files can be done more easily. 
It was called from the command line with 
`./run_smirky.bash 10000 0.0001`
The large number of iterations was chosen to see what the sampler was capable of at this point

Final score was 95.517%
Created Type                                    Matched Reference
* 0: [*:1]                                      n0003: [$([#1]-C-[#7,#8,F,#16,Cl,Br]):1]
* 1932: [#1!X2:1]-[#6!H1](-[#6!H0])-[#6!X2]     n0002: [$([#1]-C):1]
* 6057: [#1!X2:1]-[#8!X4]-[#6!H0]               n0006: [#1$(*-[#8]):1]
* 6190: [#8H1:1]                                n0011: [#8X2+0$(*-[#1]):1]
* 5034: [#6!H0:1]                               n0008: [#6X4:1]
* 4177: [#8H0:1]                                n0010: [#8X2:1]
