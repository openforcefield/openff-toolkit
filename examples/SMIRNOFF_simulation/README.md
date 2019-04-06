## `SMIRNOFF_simulation`: Simulate a molecule in the gas phase using SMIRNOFF forcefield format

This example code ([run_molecule.py](https://github.com/openforcefield/openforcefield/blob/master/examples/SMIRNOFF_simulation/run_molecule.py)) demonstrates how to simulate a molecule in the gas phase using the SMIRNOFF forcefield format.

The essential idea is to create an openforcefield `Topology` object (as well as positions, if present) from a mol2 file:
```python
# Load molecule
from openforcefield.topology import Molecule
molecule = Molecule.from_file(mol_filename)

# Create an openforcefield Topology
from openforcefield.topology import Topology
topology = Topology.from_molecules(molecule)

# Get positions in OpenMM-compatible format
positions = molecule.positions
```
We then load a SMIRNOFF forcefield and create an OpenMM `System` object to simulate:
```python
# Load a SMIRNOFF small molecule forcefield for alkanes, ethers, and alcohols
forcefield = ForceField(offxml_filename)

# Create the OpenMM system
system = forcefield.create_system(topology)
```
