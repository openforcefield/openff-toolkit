# Example of sampling Bonds 

This example uses the the decorator files in this directory 
instead of the defaults in smarty/data/odds_files. 

Also, the script run_smirky.bash was created so that repeated experiments 
with the odds in these decorator files can be done more easily. 
It was called from the command line with 
`./run_smirky.bash 10000 0.0001`
The large number of iterations was chosen to see what the sampler was capable of at this point

Final score was 100.00%
Created Type                        Matched Reference
* 0: [*:1]~[*:2]                    b0001: [#6X4:1]-[#6X4:2]
* 8101: [#1!X2:1]~[#6!H3:2]         b0002: [#6X4:1]-[#1:2]
* 1628: [#1!X2:1]~[#8H1:2]          b0006: [#8X2:1]-[#1:2]
* 8455: [#8!H0:1]~[#6!H3:2]         b0004: [#6X4:1]-[O&X2&H1:2]
* 1275: [*:1]~[#8H0:2]              b0005: [#6X4:1]-[O&X2&H0:2]
