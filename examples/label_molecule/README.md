# Usage example for forcefield_labeler

Shows how to use forcefield_labeler to see what SMIRFF parameter numbers are assigned to specific atoms for a test molecule (which must be covered by the AlkEthOH parameter set).

## Manifest
* label_molecule.py: Creates a simple example molecule and uses the labeler to indicate what force terms are applied to which atoms.
* get_parameter_statistics.py: Pulls details on which parameters are used in each molecule, and which molecules each parameter occurs in.
