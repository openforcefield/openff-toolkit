## Using partial bond orders

In this example, we demonstrate how a SMIRNOFF force field can use partial bond order information to assign bond parameters relevant for benzene as well as singly and doubly-bonded trivalent carbon, all with a single bond entry.

The file `Frosst_AlkEthOH_extracarbons.offxml` is a SMIRNOFF FFXML file adding extra carbon parameters to cover benzene (adapted from work by Christopher Bayly) and uses partial bond orders to compress the `[#6X3]-[#6X3]` , `[#6X3]:[#6X3]` , and `[#6X3]=[#6X3]` bond parameters to a single line.

The `<Bonds>` section contains
```XML
<Bonds potential="harmonic" length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2" fractional_bondorder_method="Wiberg" fractional_bondorder_interpolation="linear">
   ...
   <Bond smirks="[#6X3:1]!#[#6X3:2]" k_bondorder1="820.0" k_bondorder2="1098" length_bondorder1="1.45" length_bondorder2="1.35"/> <!-- Christopher Bayly from parm99, Aug 2016 -->
   <Bond smirks="[#6X3:1]-[#1:2]" k="734.0" length="1.080"/> <!-- Christopher Bayly from parm99, Aug 2016 -->
</Bonds>
```
The jupyter notebook [test_partialbondorder.ipynb](https://github.com/openforcefield/openforcefield/blob/master/examples/partial_bondorder/test_partialbondorder.ipynb) illustrates the application of these parameters to benzene, and inspection of the resulting bonded parameters.
