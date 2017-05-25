# Example/test showing how a SMIRNOFF file with partial bond order information can be used to have a single line give bond parameters relevant for benzene but also singly and doubly-bonded trivalent carbon.


## Manifest
* `Frosst_AlkEthOH_extracarbons.ffxml`: SMIRNOFF FFXML file adding extra carbon parameters to cover benzene (adapted from work by Christopher Bayly, and using partial bond orders to compress the [#6X3]-[#6X3], [#6X3]:[#6X3] and [#6X3]=[#6X3] bond parameters to a single line.
* `test_partialbondorder.ipynb`: Apply these parameters to benzene and print info about the bonded parameters which are assigned.
