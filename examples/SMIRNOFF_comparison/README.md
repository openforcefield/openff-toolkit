## Compare small molecule parm@frosst energies between SMIRNOFF and AMBER `prmtop/inpcrd`

This example demonstrates how to compare parm@frosst energies of molecules from the AlkEthOH set between systems generated from a SMIRNOFF `.offxml` file and systems generated from AMBER `prmtop/inpcrd` files.

* [compare_molecule_energies.py](https://github.com/openforcefield/openforcefield/blob/master/examples/SMIRNOFF_comparison/compare_molecule_energies.py): compare energies for a single molecule
* [compare_set_energies.py](https://github.com/openforcefield/openforcefield/blob/master/examples/SMIRNOFF_comparison/compare_set_energies.py): compare energies across a set of molecules

The key idea here is to use the utility `compare_system_energies` to compare the molecule parameterized by `ForceField` with a corresponding system constructed via an AMBER `prmtop/inpcrd` pair:

```python
# Load molecule
from openmmtools.topology import Molecule
molecule = Molecule.from_file(mol_filename)

# Load forcefield
from openforcefield.typing.engines.smirnoff import ForceField
forcefield = ForceField(offxml_filename)

# Compare energies
from openforcefield.tests.utils import compare_molecule_energies
energies = compare_amber_smirnoff(prmtop_filename, inpcrd_filename, forcefield, molecule)
```
