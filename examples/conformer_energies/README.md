## Computing small molecule energies using a SMIRNOFF force field

These examples illustrate how to compute small molecule vacuum energies (for example, for different conformers of the same molecule) using a SMIRNOFF force field.

* [`conformer_energies.ipynb`](conformer_energies.ipynb) contains an example illustrating the computation of conformer energies using a SMIRNOFF force field
* [`conformer_energies.py`](conformer_energies.py) is a CLI utility for computing conformer energies using the "Sage" (openff-2.1.0) SMIRNOFF force field.
  It can be run on example data with `python conformer_energies.py -f ruxolitinib_conformers.sdf`
