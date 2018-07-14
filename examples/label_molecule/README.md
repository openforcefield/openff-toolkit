# Inspecting which SMIRNOFF parameters are assigned to specific molecules

This example illustrates how to inspect which SMIRNOFF parameters are assigned to specific atoms for a test molecule, using the AlkEthOH parameter set.

## Example files

* [label_molecule.py](https://github.com/openforcefield/openforcefield/blob/master/examples/label_molecule/label_molecule.py): Creates a simple example molecule and uses the labeler to indicate what force terms are applied to which atoms
* [get_parameter_statistics.py](https://github.com/openforcefield/openforcefield/blob/master/examples/label_molecule/get_parameter_statistics.py): Pulls details on which parameters are used in each molecule, and which molecules each parameter occurs in
