# Example of sampling Torsion 

This example uses the the decorator files in this directory 
instead of the defaults in smarty/data/odds_files. 

Also, the script run_smirky.bash was created so that repeated experiments 
with the odds in these decorator files can be done more easily. 
It was called from the command line with 
`./run_smirky.bash 10000 0.0001`
The large number of iterations was chosen to see what the sampler was capable of at this point

Final score was 92.460%
Created Type                                            Matched Reference
* 0: [*:1]~[*:2]~[*:3]~[*:4]                            t0004: [#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]
* 5669: [*:1](-[#1!X4])~[*:2]~[*:3]~[*:4]               t0003: [a,A:1]-[#6X4:2]-[#8X2:3]-[!#1:4]
* 4147: [#8!X4:1]~[*:2]~[#6X4:3]~[#6!X2:4]              t0001: [a,A:1]-[#6X4:2]-[#6X4:3]-[a,A:4]
* 4706: [#6H2:1]~[*:2]~[#6X4:3]~[#6!X2:4]               t0008: [#6X4:1]-[#6X4:2]-[#8X2:3]-[#6X4:4]
* 3363: [#1!H3:1]~[*:2](-[#1!H3])~[#6X4:3]~[#6!X2:4]    t0005: [#1:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]
* 1797: [*:1](-[#1!X2])~[*:2]~[#8H1:3]~[*:4]            t0006: [#6X4:1]-[#6X4:2]-[#8X2:3]-[#1:4]
* 5384: [#1!X4:1]~[*:2]~[#8H1:3]~[*:4]                  t0002: [a,A:1]-[#6X4:2]-[#8X2:3]-[#1:4]
* 6449: [#8!X4:1]~[*:2]~[*:3]~[#8!X4:4]                 t0010: [#8X2:1]-[#6X4:2]-[#6X4:3]-[#8X2:4]
* 7550: [#8!X4:1]~[*:2]~[#6X4:3]~[#1!H3:4]              t0012: [#1:1]-[#6X4:2]-[#6X4:3]-[OX2:4]
* 1207: [#6H2:1]~[*:2](-[#1!X4])~[#6X4:3]~[#6!X2:4]     t0007: [#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]
