# Example of sampling Angles 

This example uses the the decorator files in this directory 
instead of the defaults in smarty/data/odds_files. 

Also, the script run_smirky.bash was created so that repeated experiments 
with the odds in these decorator files can be done more easily. 
It was called from the command line with 
`./run_smirky.bash 10000 0.0001`
The large number of iterations was chosen to see what the sampler was capable of at this point

Final score was 82.655%
Created Type                                        Matched Reference
* 0: [*:1]~[*:2]~[*:3]                              a0002: [#1:1]-[#6X4:2]-[#1:3]
* 7394: [*:1]~[*:2]~[*:3]-[#1!H3]                   a0001: [a,A:1]-[#6X4:2]-[a,A:3]
* 2381: [#6H2:1]~[*:2]~[*:3]-[#1!H3]                a0003: [#6X4:1]-[#6X4:2]-[#6X4:3]
* 6359: [*:1]~[#8!H3:2]~[*:3]-[#1!H3]               a0005: [#6X4:1]-[#8X2:2]-[#1:3]
* 6022: [*:1](-[#1!H2])~[#8!H3:2]~[*:3]-[#1!H3]     a0006: [#6X4:1]-[#8X2:2]-[#6X4:3]
