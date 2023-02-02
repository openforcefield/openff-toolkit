# Release History

Releases follow the `major.minor.micro` scheme recommended by [PEP440](https://www.python.org/dev/peps/pep-0440/#final-releases), where

* `major` increments denote a change that may break API compatibility with previous `major` releases
* `minor` increments add features but do not break API compatibility
* `micro` increments represent bugfix releases or improvements in documentation

## Current development

### New features
- [PR #1502](https://github.com/openforcefield/openff-toolkit/pull/1502): Adds Gasteiger charge computation using the RDKit backend.
- [PR #1498](https://github.com/openforcefield/openff-toolkit/pull/1498): `Molecule.remap()` now supports partial mappings with the `partial` argument.
- [PR #1528](https://github.com/openforcefield/openff-toolkit/pull/1528): `Topology.box_vectors` are can now be set with `openmm.unit.Quantity`s, which are internally converted.

### Behavior changes
- [PR #1498](https://github.com/openforcefield/openff-toolkit/pull/1498): New, more complete, and more descriptive errors for `Molecule.remap()`.
- [PR #1525](https://github.com/openforcefield/openff-toolkit/pull/1525): Some unrelesed force fields previously accesible from `"openff/toolkit/data/test_forcefields/"` are no longer implicitly available to the `ForceField` constructor.

### Improved documentation and warnings
- [PR #1498](https://github.com/openforcefield/openff-toolkit/pull/1498): Improved documentation for `Molecule.remap()`, `Molecule.from_smiles()`, and `Molecule.from_mapped_smiles()`, emphasizing the relationships between these methods. In particular, the documentation now clearly states that `from_smiles()` will not reorder atoms based on SMILES atom mapping.
- [PR #1525](https://github.com/openforcefield/openff-toolkit/pull/1525): Improves reporting failures when loading force fields.
- [PR #1513](https://github.com/openforcefield/openff-toolkit/pull/1513): Improves error messages and documentation around supported aromaticity models (currently only "OEAroModel_MDL").


## 0.12.0

### New features
- [PR #1484](https://github.com/openforcefield/openff-toolkit/pull/1484): A `positions` argument has been added to `Topology.from_openmm()` and `Topology.from_mdtraj()`, which allows the topology's positions to be set more conveniently.
- [PR #1468](https://github.com/openforcefield/openff-toolkit/pull/1468): Track which `ParameterHandler`s are loaded as plugins.

### Behavior changes
- [PR #1481](https://github.com/openforcefield/openff-toolkit/pull/1481):
  Removes `compute_partial_charges_am1bcc`, which was deprecated in 0.11.0.
- [PR #1466](https://github.com/openforcefield/openff-toolkit/pull/1466):
  Replaces the use of `collections.OrderedDict` throughout the toolkit with
  the built-in `dict`.
  `attach_units`, `detach_units`, and `extract_serialized_units_from_dict` have been removed from
  `openff.toolkit.utils.utils`.
- [PR #1472](https://github.com/openforcefield/openff-toolkit/pull/1472):
  Removes [`ParameterHandler._VALENCE_TYPE`] and the same attribute of its subclasses, which were
  previously not used. Also deprecates `ChemicalEnvironment` and, by extension, the
  `openff.toolkit.typing.chemistry` submodule.

[`ChemicalEnvironment`]: ChemicalEnvironment
[`ParameterHandler._VALENCE_TYPE`]: ParameterHandler._VALENCE_TYPE


### Bugfixes
- [PR #1476](https://github.com/openforcefield/openff-toolkit/pull/1476): Fixes
  [#1475](https://github.com/openforcefield/openff-toolkit/issues/1475) by also registering
  a `ParameterHandler`'s class when calling `ForceField.register_parameter_handler`.
- [PR #1480](https://github.com/openforcefield/openff-toolkit/pull/1480): Fixes
  [#1479](https://github.com/openforcefield/openff-toolkit/issues/1479) by requiring that `Atom.atomic_number` is a positive integer.
- [PR #1494](https://github.com/openforcefield/openff-toolkit/pull/1494): Fixes
  [#1493](https://github.com/openforcefield/openff-toolkit/issues/1493) in which some OFFXML file
  contents were parsed to `unit.Quantity` objects despite not representing physical quantities.

[`Atom.atomic_number`]: Atom.atomic_number

### Improved documentation and warnings
- [PR #1484](https://github.com/openforcefield/openff-toolkit/pull/1484): The docstrings for `Topology.from_openmm()` and `Topology.from_mdtraj()` have been improved.
- [PR #1483](https://github.com/openforcefield/openff-toolkit/pull/1483): Simplified and clarified errors and warnings related to undefined stereochemistry with RDKit.

## 0.11.4 Bugfix release

### Behavior changes
- [PR #1462](https://github.com/openforcefield/openff-toolkit/pull/1462): Makes residue
  numbers added by `Molecule.perceive_residues` strings (previously they were ints), to 
  match the behavior of `Topology.from_openmm` and other hierarchy info-setting methods. 

### Bugfixes
- [PR #1459](https://github.com/openforcefield/openff-toolkit/pull/1459): Fixes
  [#1430](https://github.com/openforcefield/openff-toolkit/issues/1430), where 
  `Topology.from_openmm` would mis-assign atom names (and probably also 
  hierarchy metadata as well).
- [PR #1462](https://github.com/openforcefield/openff-toolkit/pull/1462): Fixes
  [#1461](https://github.com/openforcefield/openff-toolkit/issues/1461), where the 
  default `Molecule.residues` iterator wouldn't sort by residue number correctly 
  when residue information was added by `Molecule.perceive_residues`.


## 0.11.3 Bugfix release


- [PR #1460](https://github.com/openforcefield/openff-toolkit/pull/1460): Disables error causing 
  [Issue #1432](https://github.com/openforcefield/openff-toolkit/issues/1432), where 
  `Molecule.from_polymer_pdb` would sometimes issue stereochemistry errors when reading valid 
  PDBs using the RDKit backend.  


### Bugfixes
- [PR #1436](https://github.com/openforcefield/openff-toolkit/pull/1436): Fix a small bug introduced in 0.11.2, where running with OpenEye installed but not licensed could lead to a crash.
- [PR #1444](https://github.com/openforcefield/openff-toolkit/pull/1444): Update for pint 0.20.

### Examples updates
- [PR #1447](https://github.com/openforcefield/openff-toolkit/pull/1447): Fixed units of tolerance used in OpenMM minimization in Toolkit Showcase example notebook (from @ziyuanzhao2000)

### Improved documentation and warnings
- [PR #1442](https://github.com/openforcefield/openff-toolkit/pull/1442): Doctests added to CI, leading to numerous fixed docstrings and examples therein.

### Miscellaneous
- [PR #1413](https://github.com/openforcefield/openff-toolkit/pull/1413): Remove some large and unused data files from the test suite.
- [PR #1434](https://github.com/openforcefield/openff-toolkit/pull/1434): Remove dependency on `typing_extensions`.

## 0.11.2 Bugfix release

### Behavior changes
- [PR #1421](https://github.com/openforcefield/openff-toolkit/pull/1421): Allow `Molecule.from_rdkit()` to load D- and F- block radicals, which cannot have implicit hydrogens.

### Bug fixes
- [PR #1417](https://github.com/openforcefield/openff-toolkit/pull/1417): Ensure the properties dict is copied when a `Molecule` is.

### Improved documentation and warnings
- [PR #1426](https://github.com/openforcefield/openff-toolkit/pull/1426): A warning about OpenEye Toolkits being unavailable is only emitted when they are installed but the license file is not found.

## 0.11.1 Minor release forbidding loading radicals

### Behavior changes
- [PR #1398](https://github.com/openforcefield/openff-toolkit/pull/1398): Updates the [`Bond.bond_order`] setter to only accept int values.
- [PR #1236](https://github.com/openforcefield/openff-toolkit/pull/1236): [`from_rdkit`] and [`from_openeye`] now 
  raise an `RadicalsNotSupportedError` when loading radicals. It's not clear that the OpenFF Toolkit was ever safely 
  handling radicals - they appear to be the root cause of many instances of unintended hydrogen addition and other 
  connection table changes. If this change affects a workflow that was previously working correctly, please let us 
  know on [this issue](https://github.com/openforcefield/openff-toolkit/issues/1075) so we can refine this behavior. 

### Examples changed
- [PR #1236](https://github.com/openforcefield/openff-toolkit/pull/1236): `examples/check_dataset_parameter_coverage` has
  been deprecated. 

[`Bond.bond_order`]: Bond.bond_order
[`from_rdkit`]: Molecule.from_rdkit
[`from_openeye`]: Molecule.from_openeye

### Bug fixes
- [PR #1400](https://github.com/openforcefield/openff-toolkit/pull/1400): Fixes a bug where `Molecule.from_pdb_and_smiles` could incorrectly order coordinates.
- [PR #1404](https://github.com/openforcefield/openff-toolkit/pull/1404): Support default hierarchy schemes in outputs of `Molecule.from_pdb_and_smiles()` and `Topology.from_openmm()`

## 0.11.0 Major release adding support for proteins and refactoring the Topology class.

## Migration guide

### New [`Molecule.from_polymer_pdb()`] method for loading proteins from PDB files

The Toolkit now supports loading protein PDB files through the [`Molecule.from_polymer_pdb()`] class method. For now, PDB files must consist of only a single protein molecule composed only of the 20 standard amino acids, their common protonated and deprotonated conjugates, and the N-methyl and acetyl caps.

[`Molecule.from_polymer_pdb()`]: Molecule.from_polymer_pdb

### Important API points re-exported from `openff.toolkit`

A number of commonly used API points have been re-exported from the package root. This should make using the Toolkit simpler for most people. The previous API points remain available. These API points are lazy-loaded so that parts of the toolkit can still be loaded without loading the entire code-base.

The most important of these are the `ForceField`, `Molecule`, and `Topology` classes:

```diff
- from openff.toolkit.typing.engines.smirnoff import ForceField
- from openff.toolkit.topology import Molecule, Topology
+ from openff.toolkit import ForceField, Molecule, Topology
```

A number of other useful API points are also available through this mechanism:

```diff
- from openff.toolkit.typing.engines.smirnoff import get_available_force_fields
- from openff.toolkit.utils.toolkits import (
-     GLOBAL_TOOLKIT_REGISTRY,
-     AmberToolsToolkitWrapper,
-     BuiltInToolkitWrapper,
-     OpenEyeToolkitWrapper,
-     RDKitToolkitWrapper,
-     ToolkitRegistry,
- )
+ from openff.toolkit import (
+     get_available_force_fields,
+     GLOBAL_TOOLKIT_REGISTRY,
+     AmberToolsToolkitWrapper,
+     BuiltInToolkitWrapper,
+     OpenEyeToolkitWrapper,
+     RDKitToolkitWrapper,
+     ToolkitRegistry,
+ )
```

The `topology`, `typing`, and `utils` modules can also be lazy loaded after importing only the top-level module:

```diff
- import openff.toolkit.topology
+ import openff.toolkit
  atom = openff.toolkit.topology.Atom()
```

### Units

The use of OpenMM units has been replaced by the new OpenFF Units package, based on Pint.

Import the [unit registry](https://pint.readthedocs.io/en/stable/developers_reference.html?highlight=unitregistry#pint.UnitRegistry) provided by `openff-units`:

```
from openff.units import unit
```

Create a `unit.Quantity` object:
```
value = unit.Quantity(1.0, unit.nanometer)  # or 1.0 * unit.nanometer
```

Inspect the value and unit of this quantity:
```
print(value.magnitude)  # or value.m
# 1.0
print(value.units)
# <Unit('nanometer')>
```

Convert to compatible units:
```
converted = value.to(unit.angstrom)
print(converted)
# 10.0 <Unit('angstrom')>
```

Report the value in compatible units:
```
print(value.m_as(unit.angstrom))  # Note that value.magnitude_as() does not exist
# 10.0 <Unit('angstrom')>
```

Convert to and from OpenMM quantities:
```
from openff.units.openmm import to_openmm, from_openmm
value_openmm = to_openmm(value)
print(value_openmm)
# Quantity(value=1.0, unit=nanometer)
print(type(value_openmm))
# 1.0 <Unit('nanometer')>
value_roundtrip = from_openmm(value_openmm)
print(value_roundtrip)
# 1.0 <Unit('nanometer')>
```

#### Breaking change: Removal of `openff.toolkit.utils.check_units_are_compatible()`

The `openff.toolkit.utils.check_units_are_compatible()` function has been removed. Use [`openff.units.Quantity.is_compatible_with()`] and [`openff.units.openmm.from_openmm()`] instead:

```diff
- check_units_are_compatible("length", length, openmm.unit.angstrom)
+ from_openmm(length).is_compatible_with(openff.units.unit.angstrom)
```

[`openff.units.Quantity.is_compatible_with()`]: openff.units.Quantity.is_compatible_with
[`openff.units.openmm.from_openmm()`]: openff.units.openmm.from_openmm

### Breaking change: Interchange now responsible for system parametrization

Code for applying parameters to topologies has been removed from the Toolkit. This is now the responsibility of [OpenFF Interchange]. This change improves support for working with parametrized systems (through the [`Interchange`] class), and adds support for working with simulation engines other than OpenMM.

The [`ForceField.create_interchange()`] method has been added, and the [`ForceField.create_openmm_system()`] method now uses Interchange under the hood.

As part of this change, the `UnsupportedKeywordArgumentsError` has been removed;
passing unknown arguments to `create_openmm_system` now raises a `TypeError`, as is normal in Python.

The following classes and methods have been **removed** from `openff.toolkit.typing.engines.smirnoff.parameters`:
- `NonintegralMoleculeChargeException`
- `NonbondedMethod`
- `ParameterHandler.assign_parameters()`
- `ParameterHandler.postprocess_system()`
- `ParameterHandler.check_partial_bond_orders_from_molecules_duplicates()`
- `ParameterHandler.assign_partial_bond_orders_from_molecules()`

In addition, the `ParameterHandler.create_force()` method has been deprecated and its functionality has been removed. It will be removed in a future release. 

The `return_topology` argument of `create_openmm_system` has also been deprecated, and will be removed in 0.12.0. To create an OpenMM topology, use `Interchange`:

```diff
- omm_sys, off_top = force_field.create_openmm_system(
-     topology,
-     return_topology=True,
- )
- omm_top = off_top.to_openmm()
+ interchange = force_field.create_interchange(topology)
+ omm_sys = interchange.to_openmm(combine_nonbonded_forces=True)
+ omm_top = interchange.to_openmm_topology()
```

If you need access to the modified OpenFF topology for some other reason, create an `Interchange` and retrieve it there:

```diff
- omm_sys, off_top = force_field.create_openmm_system(
-     topology,
-     return_topology=True,
- )
+ interchange = force_field.create_interchange(topology)
+ off_top = interchange.topology
+ omm_sys = interchange.to_openmm(combine_nonbonded_forces=True)
```

[`Interchange`]: openff.interchange.Interchange
[`ForceField.create_interchange()`]: openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.create_interchange
[`ForceField.create_openmm_system()`]: openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system
[OpenFF Interchange]: https://docs.openforcefield.org/interchange

### Breaking change: `Topology` molecule representation

`Topology` objects now store complete copies of their constituent `Molecule` objects, rather than using simplified classes specific to `Topology`. This dramatically simplifies the code base and allows the use of the full `Molecule` API on molecules inside topologies.

The following classes have been **removed**:
- `TopologyAtom` (use [`Atom`](Atom) instead)
- `TopologyBond` (use [`Bond`](Bond) instead)
- `TopologyMolecule` (use [`Molecule`](Molecule) instead)

The following properties have been **deprecated** and will be removed in a future release:
- `Topology.n_topology_atoms` (use [`Topology.n_atoms`](Topology.n_atoms) instead)
- `Topology.topology_atoms` (use [`Topology.atoms`](Topology.atoms) instead)
- `Topology.n_topology_bonds` (use [`Topology.n_bonds`](Topology.n_bonds) instead)
- `Topology.topology_bonds` (use [`Topology.bonds`](Topology.bonds) instead)
- `Topology.n_topology_particles` (use [`Topology.n_particles`](Topology.n_particles) instead)
- `Topology.topology_particles` (use [`Topology.particles`](Topology.particles) instead)
- `Topology.reference_molecules` (use [`Topology.unique_molecules`](Topology.unique_molecules) instead)
- `Topology.n_reference_molecules` (use [`Topology.n_unique_molecules`](Topology.n_unique_molecules) instead)
- `Topology.n_topology_molecules` (use [`Topology.n_molecules`](Topology.n_molecules) instead)
- `Topology.topology_molecules` (use [`Topology.molecules`](Topology.molecules) instead)
- `Topology.n_particles` (use [`Topology.n_atoms`](Topology.n_atoms) instead)
- `Topology.particles` (use [`Topology.atoms`](Topology.atoms) instead)
- `Topology.particle_index` (use [`Topology.atom_index`](Topology.atom_index) instead)

In addition, the [`Topology.identical_molecule_groups`] property has been added, to facilitate iterating over copies of isomorphic molecules in a `Topology`.

[`Topology.identical_molecule_groups`]: Topology.identical_molecule_groups

### Breaking change: Removed virtual site handling from topologies

To maintain a clear distinction between a model and the chemistry it represents, virtual site handling has been removed from the Toolkit's `Topology` and `Molecule` classes. Virtual site support remains in the force field side of the toolkit, but creating virtual sites for particular molecules is now the responsibility of OpenFF Interchange. This allows the same `Topology` to be used for force fields that use different virtual sites; for example, a topology representing a solvated protein might be parametrized separately with a 3-point and 4-point water model.

As part of this change, the distinction between `Atom` and `Particle` is deprecated. The `Particle` class will be removed in a future release.

The following classes have been **removed**:
- `BondChargeVirtualSite`
- `DivalentLonePairVirtualSite`
- `MonovalentLonePairVirtualSite`
- `TrivalentLonePairVirtualSite`
- `VirtualParticle`
- `VirtualSite`
- `TopologyVirtualParticle`
- `TopologyVirtualSite`

The following methods and properties have been **removed**:
- `Atom.add_virtual_site()`
- `Atom.virtual_sites`
- `FrozenMolecule.compute_virtual_site_positions_from_conformer()`
- `FrozenMolecule.compute_virtual_site_positions_from_atom_positions()`
- `FrozenMolecule.n_virtual_sites`
- `FrozenMolecule.n_virtual_particles`
- `FrozenMolecule.virtual_sites()`
- `Molecule.add_bond_charge_virtual_site()`
- `Molecule.add_monovalent_lone_pair_virtual_site()`
- `Molecule.add_divalent_lone_pair_virtual_site()`
- `Molecule.add_trivalent_lone_pair_virtual_site()`
- `Molecule.add_bond_charge_virtual_site()`
- `Topology.n_topology_virtual_sites`
- `Topology.topology_virtual_sites`
- `Topology.virtual_site()`
- `Topology.add_particle()`

The following properties have been **deprecated** and will be removed in a future release:

- `Molecule.n_particles` (use [`Molecule.n_atoms`](Molecule.n_atoms) instead)
- `Molecule.particles` (use [`Molecule.atoms`](Molecule.atoms) instead)
- `Molecule.particle` (use [`Molecule.atom`](Molecule.atom) instead)
- `Molecule.particle_index` (use [`Topology.n_atoms`](Topology.n_atoms) instead)

### Atom metadata and hierarchy schemes for iterating over residues, chains, etc.

The new `Atom.metadata` attribute is a dictionary that can store arbitrary metadata. Atom metadata commonly includes residue names, residue sequence numbers, chain identifiers, and other metadata that is not essential to the functioning of the Toolkit. Metadata can then be passed on when a `Molecule` is converted to another package; see [](users/molecule_conversion).

Metadata can also support iteration through the [`HierarchyScheme`](openff.toolkit.topology.HierarchyScheme) class. A hierarchy scheme is defined by some uniqueness criteria. Iterating over the scheme iterates over groups of atoms that have identical metadata values for the defined uniqueness criteria. For more information, see the API docs for [`HierarchyScheme`](openff.toolkit.topology.HierarchyScheme) and its related methods.

### Breaking change: Removed `Topology.charge_model` and `Topology.fractional_bond_order_model`

Due to flaws in previous versions of the OFF Toolkit, these properties never had an effect on the assigned parameters. To resolve this bug and maintain a clear distinction between a model and the chemistry it represents, the `Topology.charge_model` and `Topology.fractional_bond_order_model` properties have been removed. Charge models and FBOs are now the responsibility of the ForceField.

### Breaking change: Removed `Atom.element`

`Atom.element` has been removed to reduce our dependency on OpenMM for core functions:

```diff
- atomic_number = atom.element.atomic_number
+ atomic_number = atom.atomic_number
- atom_mass = atom.element.mass
+ atom_mass = atom.mass
- atom_elem_symbol = atom.element.symbol
+ atom_elem_symbol = atom.symbol
```

### `Topology.to_file()`

The [`Topology.to_file()`](openff.toolkit.topology.Topology.to_file) method has been significantly revised, including three breaking changes.

#### Breaking change: `filename` argument renamed `file`

The `filename` argument has been renamed `file`, and now supports file-like objects in addition to file names:

```diff
-  topology.to_file(filename="out.pdb", positions=xyz)
+  topology.to_file(file="out.pdb", positions=xyz)
```

#### Breaking change: Atom names guaranteed unique per residue by default

The default behavior is now to ensure that atom names are unique within a residue, rather than within a molecule. The `ensure_unique_atom_names` argument has been added to control this behavior. The previous behavior can be achieved by passing `True` to `ensure_unique_atom_names`:

```diff
- topology.to_file("out.pdb", xyz)
+ topology.to_file("out.pdb", xyz, ensure_unique_atom_names=True)
```

The `ensure_unique_atom_names` argument can also take the name of a `HierarchyScheme`, in which case atom names will be unique within the elements of that scheme (instead of within the atoms of a molecule). If the scheme is missing from a molecule, atom names will be unique within that molecule. The default value of this argument is `"residues"` to preserve atom names from the PDB.

#### Breaking change: `keepIds` argument renamed `keep_ids`

The `keepIds` argument has been renamed to the more Pythonic `keep_ids`. Its behavior and position in the argument list has not changed.

```diff
- topology.to_file("out.pdb", xyz, keepIds=True)
+ topology.to_file("out.pdb", xyz, keep_ids=True)
```

#### Non-breaking changes

In addition to these breaking changes, the `positions` argument is now optional. If it is not provided, positions will be taken from the first conformer of each molecule in the topology. If any molecule has no conformers, an error will be raised.

### Positions in topologies

The [`Topology.get_positions()`] and [`Topology.set_positions()`] methods have been added to facilitate working with coordinates in topologies. A topology's positions are defined by the zeroth conformer of each molecule. If any molecule lacks conformers, the entire topology has no positions.

[`Topology.get_positions()`]: Topology.get_positions
[`Topology.set_positions()`]: Topology.set_positions

### Parameter types moved out of handler classes

To facilitate their discovery and documentation, re-exports for the `ParameterType` classes have been added to the `openff.toolkit.typing.engines.smirnoff.parameters` module. Previously, they were accessible only within their associated `ParameterHandler` classes. This is not a breaking change.

```diff
- from openff.toolkit.typing.engines.smirnoff.parameters import BondHandler
- BondType = BondHandler.BondType
+ from openff.toolkit.typing.engines.smirnoff.parameters import BondType
```

### Breaking change: `MissingDependencyError` renamed `MissingPackageError`

The `MissingDependencyError` exception has been renamed [`MissingPackageError`] to better reflect its purpose.

```diff
  try:
      ...
- except MissingDependencyError:
+ except MissingPackageError:
      pass
```

[`MissingPackageError`]: openff.toolkit.utils.exceptions.MissingPackageError

### `compute_partial_charges_am1bcc()` deprecated

The `compute_partial_charges_am1bcc()` methods of the `Molecule`, `AmberToolsToolkitWrapper` and `OpenEyeToolkitWrapper` classes have been deprecated and will be removed in a future release. Their functionality has been incorporated into [`assign_partial_charges()`] for more consistency with other charge generation methods:

```diff
- mol.compute_partial_charges_am1bcc()
+ mol.assign_partial_charges(partial_charge_method='am1bcc')
```

[`assign_partial_charges()`]: Molecule.assign_partial_charges


### Additional changes and bugfixes

- [PR #1105](https://github.com/openforcefield/openff-toolkit/pull/1105), [PR #1195](https://github.com/openforcefield/openff-toolkit/pull/1195), [PR #1301](https://github.com/openforcefield/openff-toolkit/pull/1301), [PR #1331](https://github.com/openforcefield/openff-toolkit/pull/1331), [PR #1322](https://github.com/openforcefield/openff-toolkit/pull/1322), [PR #1372](https://github.com/openforcefield/openff-toolkit/pull/1372): Add `Molecule.from_polymer_pdb`
- [PR #1377](https://github.com/openforcefield/openff-toolkit/pull/1377): Adds 
  `Topology.unique_molecules`, which largely replaces `Topology.reference_molecules`. 
- [PR #1313](https://github.com/openforcefield/openff-toolkit/pull/1313): Fixes 
  [Issue #1287](https://github.com/openforcefield/openff-toolkit/issues/1287), where  
  `OpenEyeToolkitWrapper.assign_partial_charges` didn't request symmetrized charges when
  the charge model was set to `AM1-Mulliken`.
- [PR #1348](https://github.com/openforcefield/openff-toolkit/pull/1348): Allows
  `pathlib.Path` objects to be passed to
  [`Molecule.from_file`](openff.toolkit.topology.Molecule.from_file).
- [PR #1276](https://github.com/openforcefield/openff-toolkit/pull/1276): Removes the
  `use_interchange` argument to
  [`create_openmm_system`](openff.toolkit.typing.engines.smirnoff.ForceField.create_openmm_system).
  Deletes the `create_force` and `postprocess_system` methods of `ParameterHandler`
  [`ParameterHandler.create_force`](openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler.create_force),
  [`ParameterHandler.postprocess_system`](openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler.postprocess_system) and other methods related to creating OpenMM systems and forces. This is now handled in Interchange.
- [PR #1303](https://github.com/openforcefield/openff-toolkit/pull/1303): Deprecates `Topology.particles`,
  `Topology.n_particles`, `Topology.particle_index` as `Molecule` objects do not store virtual sites,
  only atoms.
- [PR #1297](https://github.com/openforcefield/openff-toolkit/pull/1297): Drops support
  for Python 3.7, following [NEP-29](https://numpy.org/neps/nep-0029-deprecation_policy.html).
- [PR #1194](https://github.com/openforcefield/openforcefield/pull/1194): Adds
  [`Topology.__add__`](openff.toolkit.topology.Topology.__add__), allowing `Topology` objects to be
  added together, including added in-place, using the `+` operator.
- [PR #1277](https://github.com/openforcefield/openff-toolkit/pull/1277): Adds support for version
  0.4 of the `<Electrostatics>` section of the SMIRNOFF specification.
- [PR #1279](https://github.com/openforcefield/openforcefield/pull/1279):
  [`ParameterHandler.version`](openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler.version)
  and the ``.version`` attribute of its subclasses is now a
  [``Version``](https://packaging.pypa.io/en/latest/version.html#packaging.version.Version)
  object. Previously it was a string, which is not safe for
  [PEP440](https://www.python.org/dev/peps/pep-0440/#final-releases)-style versioning.
- [PR #1250](https://github.com/openforcefield/openff-toolkit/pull/1250): Adds support for
  `return_topology` in the Interchange code path in
  [`create_openmm_system`](openff.toolkit.typing.engines.smirnoff.ForceField.create_openmm_system).
- [PR #964](https://github.com/openforcefield/openff-toolkit/pull/964): Adds initial implementation
  of atom metadata dictionaries.
- [PR #1097](https://github.com/openforcefield/openff-toolkit/pull/1097): Deprecates TopologyMolecule.
- [PR #1097](https://github.com/openforcefield/openff-toolkit/pull/1097): Topology.from_openmm
  is no longer guaranteed to maintain the ordering of bonds, but now explicitly guarantees that it maintains
  the order of atoms (Neither of these ordering guarantees were explicitly documented before, but this may be a
  change from the previous behavior).
- [PR #1165](https://github.com/openforcefield/openforcefield/pull/1165): Adds the boolean argument
  `use_interchange` to
  [`create_openmm_system`](openff.toolkit.typing.engines.smirnoff.ForceField.create_openmm_system)
  with a default value of False. Setting it to True routes `openmm.System` creation through
  Interchange.
- [PR #1192](https://github.com/openforcefield/openforcefield/pull/1192): Add re-exports for core classes to the new
  `openff.toolkit.app` module and re-exports for parameter types to the new `openff.toolkit.topology.parametertypes` module.
  This does not affect existing paths and gives some new, easier to remember paths to core objects.
- [PR #1198](https://github.com/openforcefield/openforcefield/pull/1198): Ensure the vdW switch
  width is correctly set and enabled.
- [PR #1213](https://github.com/openforcefield/openff-toolkit/pull/1213): Removes
  `Topology.charge_model` and `Topology.fractional_bond_order_model`.
- [PR #1140](https://github.com/openforcefield/openff-toolkit/pull/1140): Adds the `Topology.identical_molecule_groups` property, which provides a way of grouping the instances of a specific chemical species in the topology.
- [PR #1200](https://github.com/openforcefield/openforcefield/pull/1200): Fixes a bug
  ([Issue #1199](https://github.com/openforcefield/openff-toolkit/issues/428)) in which library
  charges were ignored in some force fields, including `openff-2.0.0` code name "Sage." This resulted in
  the TIP3P partial charges included Sage not being applied correctly in versions 0.10.1 and 0.10.2 of the
  OpenFF Toolkit. This regression was not present in version 0.10.0 and earlier and therefore is not
  believed to have any impact on the fitting or benchmarking of the first release of Sage (version
  2.0.0). The change causing regression only affected library charges and therefore no other
  parameter types are believed to be affected.
- [PR #1346](https://github.com/openforcefield/openff-toolkit/pull/1346): Conformer generation with RDKit
  will use `useRandomCoords=True` on a second attempt if the first attempt fails, which sometimes
  happens with large molecules.
- [PR #1277](https://github.com/openforcefield/openff-toolkit/pull/1277): Version 0.3 `<Electrostatics>`
  sections of OFFXML files will automatically be up-converted (in memory) to version 0.4 according
  to the recomendations provided in
  [OFF-EP 0005](https://openforcefield.github.io/standards/enhancement-proposals/off-ep-0005/). Note
  this means the `method` attribute is replaced by `periodic_potential`, `nonperiodic_potential`,
  and `exception_potential`.
- [PR #1277](https://github.com/openforcefield/openff-toolkit/pull/1277): Fixes a bug in which
  attempting to convert
  [`ElectrostaticsHandler.switch_width`](openff.toolkit.typing.engines.smirnoff.parameters.ElectrostaticsHandler)
  did nothing.
- [PR #1130](https://github.com/openforcefield/openforcefield/pull/1130): Running unit tests will
  no longer generate force field files in the local directory.
- [PR #1182](https://github.com/openforcefield/openforcefield/pull/1182): Removes `Atom.element`,
  thereby also removing `Atom.element.symbol`, `Atom.element.mass` and `Atom.element.atomic_number`.
  These are replaced with corresponding properties directly on the
  [`Atom`](openff.toolkit.topology.molecule.Atom) class:
  [`Atom.symbol`](openff.toolkit.topology.molecule.Atom.symbol),
  [`Atom.mass`](openff.toolkit.topology.molecule.Atom.mass), and
  [`Atom.atomic_number`](openff.toolkit.topology.molecule.Atom.atomic_number).
- [PR #1209](https://github.com/openforcefield/openforcefield/pull/1209): Fixes
  [Issue #1073](https://github.com/openforcefield/openff-toolkit/issues/1073), where the
  `fractional_bondorder_method` kwarg to the 
  [`BondHandler`](openff.toolkit.typing.engines.smirnoff.parameters.BondHandler) initializer 
  was being ignored.
- [PR #1214](https://github.com/openforcefield/openforcefield/pull/1214): A long overdue fix
  for [Issue #837](https://github.com/openforcefield/openff-toolkit/issues/837)! If OpenEye is
  available, the `ToolkitAM1BCCHandler` will use the ELF10 method to select conformers for AM1BCC
  charge assignment. 
- [PR #1160](https://github.com/openforcefield/openforcefield/pull/1160): Fixes the bug identified in
  [Issue #1159](https://github.com/openforcefield/openff-toolkit/issues/1159), in which the order of 
  atoms defining a `BondChargeVirtualSite` (and possibly other virtual sites types too) might be reversed 
  if the `match` attribute of the virtual site has a value of `"once"`.
- [PR #1231](https://github.com/openforcefield/openforcefield/pull/1231): Fixes
  [Issue #1181](https://github.com/openforcefield/openff-toolkit/issues/1181) and
  [Issue #1190](https://github.com/openforcefield/openff-toolkit/issues/1190), where in rare cases 
  double bond stereo would cause `to_rdkit` to raise an error. The transfer of double bond stereochemistry
  from OpenFF's E/Z representation to RDKit's local representation is now handled as a constraint
  satisfaction problem.
- [PR #1368](https://github.com/openforcefield/openff-toolkit/pull/1368): Adds the `Topology.get_positions()` and `Topology.set_positions()` methods for working with topology positions. Positions are represented as the first conformer of each molecule in the topology.
- [PR #1368](https://github.com/openforcefield/openff-toolkit/pull/1368): Allows setting the `ensure_unique_atom_names` argument of `Topology.to_openmm()`to the name of a hierarchy scheme, in which case atom names are guaranteed unique per element of that scheme rather than per molecule. Changes the default value to `"residues"`.
- [PR #1368](https://github.com/openforcefield/openff-toolkit/pull/1368): Adds the `ensure_unique_atom_names` argument to the `Topology.to_file()`, which mimics the same argument in `Topology.to_openmm()`. Renames the `keepIds` argument to `keep_ids`. Renames the `filename` argument to `file` and allows a file-like object to be passed instead of a filename. Makes the `positions` argument optional; if it is not given, positions are take from the first conformer of each molecule in the Topology.
- [PR #1290](https://github.com/openforcefield/openforcefield/pull/1290): Fixes
  [Issue #1216](https://github.com/openforcefield/openff-toolkit/issues/1216) by adding internal logic to handle
  the possibility that multiple vsites share the same parent atom, and makes the return value of 
  `VirtualSiteHandler.find_matches` be closer to the base class.

### Examples added

- [PR #1113](https://github.com/openforcefield/openff-toolkit/pull/1113): Updates the Amber/GROMACS
  example to use Interchange.

- [PR #1368](https://github.com/openforcefield/openff-toolkit/pull/1368): Updates the Toolkit showcase with the new polymer handling and Interchange support

### Tests updated

- [PR #1188](https://github.com/openforcefield/openff-toolkit/pull/1188): Add an `<Electrostatics>`
  section to the TIP3P force field file used in testing (`test_forcefields/tip3p.offxml`)



## 0.10.5 Bugfix release

- [PR #1252](https://github.com/openforcefield/openforcefield/pull/1252): Refactors virtual 
  site support, resolving
  [Issue #1235](https://github.com/openforcefield/openff-toolkit/issues/1235), 
  [Issue #1233](https://github.com/openforcefield/openff-toolkit/issues/1233), 
  [Issue #1222](https://github.com/openforcefield/openff-toolkit/issues/1222),
  [Issue #1221](https://github.com/openforcefield/openff-toolkit/issues/1221), and
  [Issue #1206](https://github.com/openforcefield/openff-toolkit/issues/1206).
  
  - Attempts to make virtual site handler more resilient through code simplification.
  - Virtual sites are now associated with a particular 'parent' atom, rather than with a set of atoms. In particular, when checking if a v-site has been assigned we now only check the main 'parent' atom associated with the v-site, rather than all additional orientation atoms. As an example, if a force field contained a bond-charge v-site that matches [O:1]=[C:2] and a monovalent lone pair that matches [O:1]=[C:2]-[*:3] in that order, then only the monovalent lone pair will be assigned to formaldehyde as the oxygen is the main atom that would be associated with both v-sites, and the monovalent lone pair appears later in the hierarchy. This constitutes a behaviour change over previous versions.
  - All v-site exclusion policies have been removed except for 'parents' which has been updated to match [OFF-EP 0006](https://openforcefield.github.io/standards/enhancement-proposals/off-ep-0006/).
  - checks have been added to enforce that the 'match' keyword complies with the SMIRNOFF spec.
  - Molecule virtual site classes no longer store FF data such as epsilon and sigma.
  - Sanity checks have been added when matching chemical environments for v-sites that ensure the environment looks like one of our expected test cases.
  - Fixes di- and trivalent lone pairs mixing the `:1` and `:2` indices.
  - Fixes trivalent v-site positioning.
  - Correctly splits `TopologyVirtualSite` and `TopologyVirtualParticle` so that virtual particles no longer have attributes such as `particles`, and ensure that indexing methods now work correctly.

## 0.10.4 Bugfix release

### Critical bugfixes

- [PR #1242](https://github.com/openforcefield/openforcefield/pull/1242): Fixes
  [Issue #837](https://github.com/openforcefield/openff-toolkit/issues/837).
  If OpenEye Toolkits are available,
  [`ToolkitAM1BCCHandler`](openff.toolkit.typing.engines.smirnoff.parameters.ToolkitAM1BCCHandler)
  will use the ELF10 method to select conformers for AM1-BCC charge assignment.
- [PR #1184](https://github.com/openforcefield/openforcefield/pull/1184): Fixes
  [Issue #1181](https://github.com/openforcefield/openff-toolkit/issues/1181) and
  [Issue #1190](https://github.com/openforcefield/openff-toolkit/issues/1190), where in rare cases
  double bond stereochemistry would cause
  [`Molecule.to_rdkit`](openff.toolkit.topology.Molecule.to_rdkit) to raise an error. The transfer
  of double bond stereochemistry from OpenFF's E/Z representation to RDKit's local representation is
  now handled as a constraint satisfaction problem.

## 0.10.3 Bugfix release

### Critical bugfixes

- [PR #1200](https://github.com/openforcefield/openforcefield/pull/1200): Fixes a bug
  ([Issue #1199](https://github.com/openforcefield/openff-toolkit/issues/428)) in which library
  charges were ignored in some force fields, including `openff-2.0.0` code name "Sage." This resulted in
  the TIP3P partial charges included Sage not being applied correctly in versions 0.10.1 and 0.10.2 of the
  OpenFF Toolkit. This regression was not present in version 0.10.0 and earlier and therefore is not
  believed to have any impact on the fitting or benchmarking of the first release of Sage (version
  2.0.0). The change causing the regression only affected library charges and therefore no other
  parameter types are believed to be affected.

### API breaking changes
- [PR #855](https://github.com/openforcefield/openff-toolkit/pull/855): In earlier
  versions of the toolkit, we had mistakenly made the assumption that cheminformatics
  toolkits agreed on the number and membership of rings. However we later learned that this
  was not true. This PR removes
  [`Molecule.rings`](openff.toolkit.topology.Molecule.rings) and
  [`Molecule.n_rings`](openff.toolkit.topology.Molecule.n_rings). To find rings in
  a molecule, directly use a cheminformatics toolkit after using
  [`Molecule.to_rdkit`](openff.toolkit.topology.Molecule.to_rdkit) or
  [`Molecule.to_openeye`](openff.toolkit.topology.Molecule.to_openeye).
  [`Atom.is_in_ring`](openff.toolkit.topology.Atom.is_in_ring) and
  [`Bond.is_in_ring`](openff.toolkit.topology.Bond.is_in_ring) are now methods, not properties.

### Behaviors changed and bugfixes

- [PR #1171](https://github.com/openforcefield/openforcefield/pull/1171): Failure of
  [`Molecule.apply_elf_conformer_selection()`] due to excluding all available conformations ([Issue #428](https://github.com/openforcefield/openff-toolkit/issues/428))
  now provides a better error. The `make_carboxylic_acids_cis` argument (`False` by default) has been added to
  [`Molecule.generate_conformers()`] to mitigate a common cause of this error. By setting this argument to `True` in internal use of this method, trans carboxylic
  acids are no longer generated in [`Molecule.assign_partial_charges()`] and
  [`Molecule.assign_fractional_bond_orders()`] methods (though users may still pass trans conformers in, they'll just be pruned by ELF methods). This should work around most instances
  of the OpenEye Omega bug where trans carboxylic acids are more common than they should be.

[`Molecule.apply_elf_conformer_selection()`]: openff.toolkit.topology.Molecule.apply_elf_conformer_selection
[`Molecule.generate_conformers()`]: openff.toolkit.topology.Molecule.generate_conformers
[`Molecule.assign_partial_charges()`]: openff.toolkit.topology.Molecule.assign_partial_charges
[`Molecule.assign_fractional_bond_orders()`]: openff.toolkit.topology.Molecule.assign_fractional_bond_orders

### Behaviors changed and bugfixes
- [PR #1185](https://github.com/openforcefield/openff-toolkit/pull/1185):
  Removed length check in ValenceDict and fixed checking the permutations of dihedrals

### Improved documentation and warnings
- [PR #1172](https://github.com/openforcefield/openff-toolkit/pull/1172): Adding
  discussion about constraints to the FAQ
- [PR #1173](https://github.com/openforcefield/openforcefield/pull/1173): Expand
  on the SMIRNOFF section of the toolkit docs
- [PR #855](https://github.com/openforcefield/openff-toolkit/pull/855): Refactors
  [`Atom.is_in_ring`](openff.toolkit.topology.Atom.is_in_ring) and
  [`Bond.is_in_ring`](openff.toolkit.topology.Bond.is_in_ring) to use corresponding
  functionality in OpenEye and RDKit wrappers.

### API breaking changes
- [PR #855](https://github.com/openforcefield/openff-toolkit/pull/855): Removes
  [`Molecule.rings`](openff.toolkit.topology.Molecule.rings) and
  [`Molecule.n_rings`](openff.toolkit.topology.Molecule.n_rings). To find rings in
  a molecule, directly use a cheminformatics toolkit after using
  [`Molecule.to_rdkit`](openff.toolkit.topology.Molecule.to_rdkit) or
  [`Molecule.to_openeye`](openff.toolkit.topology.Molecule.to_openeye).
  [`Atom.is_in_ring`](openff.toolkit.topology.Atom.is_in_ring) and
  [`Bond.is_in_ring`](openff.toolkit.topology.Bond.is_in_ring) are now methods, not properties.


## 0.10.2 Bugfix release

### API-breaking changes

- [PR #1118](https://github.com/openforcefield/openforcefield/pull/1118):
  [`Molecule.to_hill_formula`](openff.toolkit.topology.Molecule.to_hill_formula) is now a class method
  and no longer accepts input of NetworkX graphs.

### Behaviors changed and bugfixes

- [PR #1160](https://github.com/openforcefield/openforcefield/pull/1160): Fixes a major bug identified in
  [Issue #1159](https://github.com/openforcefield/openff-toolkit/issues/1159), in which the order of
  atoms defining a `BondChargeVirtualSite` (and possibly other virtual sites types too) might be reversed
  if the `match` attribute of the virtual site has a value of `"once"`.
- [PR #1130](https://github.com/openforcefield/openforcefield/pull/1130): Running unit tests will
  no longer generate force field files in the local directory.
- [PR #1148](https://github.com/openforcefield/openforcefield/pull/1148): Adds a new exception
  [`UnsupportedFileTypeError`](openff.toolkit.utils.exceptions.UnsupportedFileTypeError) and
  descriptive error message when attempting to use
  [`Molecule.from_file`](openff.toolkit.topology.Molecule.from_file) to parse XYZ/`.xyz` files.
- [PR #1153](https://github.com/openforcefield/openforcefield/pull/1153): Fixes
  [Issue #1152](https://github.com/openforcefield/openff-toolkit/issues/1052) in which running
  [`Molecule.generate_conformers`](openff.toolkit.topology.Molecule.generate_conformers)
  using the OpenEye backend would use the stereochemistry from an existing conformer instead 
  of the stereochemistry from the molecular graph, leading to undefined behavior if the molecule had a 2D conformer. 
- [PR #1158](https://github.com/openforcefield/openff-toolkit/pull/1158): Fixes the default
  representation of [`Molecule`](openff.toolkit.topology.Molecule) failing in Jupyter notebooks when
  NGLview is not installed.
- [PR #1151](https://github.com/openforcefield/openforcefield/pull/1151): Fixes 
  [Issue #1150](https://github.com/openforcefield/openff-toolkit/issues/1150), in which calling 
  [`Molecule.assign_fractional_bond_orders`](openff.toolkit.topology.Molecule.assign_fractional_bond_orders)
  with all default arguments would lead to an error as a result of trying to lowercase `None`.
- [PR #1149](https://github.com/openforcefield/openforcefield/pull/1149):
  [`TopologyAtom`](openff.toolkit.topology.TopologyAtom),
  [`TopologyBond`](openff.toolkit.topology.TopologyBond), and
  [`TopologyVirtualSite`](openff.toolkit.topology.TopologyVirtualSite) now properly reference their
  reference molecule from their `.molecule` attribute.
- [PR #1155](https://github.com/openforcefield/openforcefield/pull/1155): Ensures big-endian
  byte order of NumPy arrays when serialized to dictionaries or files formats except JSON.
- [PR #1163](https://github.com/openforcefield/openforcefield/pull/1163): Fixes the bug identified in
  [Issue #1161](https://github.com/openforcefield/openff-toolkit/issues/1161), which was caused by the use
  of the deprecated `pkg_resources` package. Now the recommended `importlib_metadata` package is used instead.


### Breaking changes
- [PR #1118](https://github.com/openforcefield/openforcefield/pull/1118):
  [`Molecule.to_hill_formula`](openff.toolkit.topology.Molecule.to_hill_formula) is now a class method
  and no longer accepts input of NetworkX graphs.
- [PR #1156](https://github.com/openforcefield/openforcefield/pull/1156): Removes `ParseError` and
  `MessageException`, which has been deprecated since version 0.10.0.

### Examples added

- [PR #1113](https://github.com/openforcefield/openff-toolkit/pull/1113): Updates the Amber/GROMACS
  example to use Interchange.

## 0.10.1 Minor feature and bugfix release

### Behaviors changed and bugfixes

- [PR #1096](https://github.com/openforcefield/openforcefield/pull/1096): Atom names generated by
  [`Molecule.generate_unique_atom_names`](openff.toolkit.topology.Molecule.generate_unique_atom_names)
  are now appended with an `"x"`. See the linked issue for more details.
- [PR #1050](https://github.com/openforcefield/openforcefield/pull/1050): In
  [`Molecule.generate_conformers`](openff.toolkit.topology.Molecule.generate_conformers), a single
  toolkit wrapper failing to generate conformers is no longer fatal, but if all wrappers in a registry
  fail, then a `ValueError` will be raised. This mirrors the behavior of
  [`Molecule.assign_partial_charges`](openff.toolkit.topology.Molecule.assign_partial_charges).
- [PR #1050](https://github.com/openforcefield/openforcefield/pull/1050): Conformer generation
  failures in
  [`OpenEyeToolkitWrapper.generate_conformers`](openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.generate_conformers), and
  [`RDKitToolkitWrapper.generate_conformers`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.generate_conformers)
  now each raise
  [`openff.toolkit.utils.exceptions.ConformerGenerationError`](openff.toolkit.utils.exceptions.ConformerGenerationError)
  if conformer generation fails. The same behavior occurs in
  [`Molecule.generate_conformers`](openff.toolkit.topology.Molecule.generate_conformers), but only
  when the ``toolkit_registry`` argument is a
  [`ToolkitWrapper`](openff.toolkit.utils.toolkits.ToolkitWrapper), not when it is a
  [`ToolkitRegistry`](openff.toolkit.utils.toolkits.ToolkitRegistry). 
- [PR #1046](https://github.com/openforcefield/openforcefield/pull/1046): Changes OFFXML output to
  replace tabs with 4 spaces to standardize representation in different text viewers. 
- [PR #1001](https://github.com/openforcefield/openff-toolkit/pull/1001): RDKit `Mol` objects 
  created through the [`Molecule.to_rdkit()`](openff.toolkit.topology.Molecule.to_rdkit)
  method have the `NoImplicit` property set to `True` on all atoms. This prevents RDKit from
  incorrectly adding hydrogen atoms to to molecule. 
- [PR #1058](https://github.com/openforcefield/openforcefield/pull/1058): Removes the unimplemented methods
  [`ForceField.create_parmed_structure`](openff.toolkit.typing.engines.smirnoff.ForceField.create_parmed_structure),
  [`Topology.to_parmed`](openff.toolkit.topology.Topology.to_parmed), and
  [`Topology.from_parmed`](openff.toolkit.topology.Topology.from_parmed).
- [PR #1065](https://github.com/openforcefield/openforcefield/pull/1065): The example `conformer_energies.py` script
  now uses the Sage 2.0.0 force field.
- [PR #1036](https://github.com/openforcefield/openforcefield/pull/1036): SMARTS matching
  logic for library charges was updated to use only one unique match instead of
  enumerating all possible matches. This results in faster matching, particularly
  with larger molecules. No adverse side effects
  were found in testing, but bad behavior may possibly exist in some unknown cases.
  Note that the default behavior for other parameter handlers was not updated.
- [PR #1001](https://github.com/openforcefield/openff-toolkit/pull/1001): Revamped the 
  [`Molecule.visualize()`](openff.toolkit.topology.Molecule.visualize) method's `rdkit` 
  backend for more pleasing and idiomatic 2D visualization by default.
- [PR #1087](https://github.com/openforcefield/openff-toolkit/pull/1087): Fixes
  [Issue #1073](https://github.com/openforcefield/openff-toolkit/issues/1073) in which
  [`Molecule.__repr__`](openff.toolkit.topology.Molecule.__repr__) fails if the molecule can not be represented as 
  a SMILES pattern. Now, if SMILES generation fails, the molecule will be described by its Hill formula.
- [PR #1052](https://github.com/openforcefield/openff-toolkit/pull/1052): Fixes
  [Issue #986](https://github.com/openforcefield/openff-toolkit/issues/986)
  by raising a subclass of `AttributeError` in
  `_ParameterAttributeHandler.__getattr__`
- [PR #1030](https://github.com/openforcefield/openforcefield/pull/1030): Fixes a bug
  in which the expectations for capitalization for values of `bond_order_model` attributes and 
  keywords are inconsistent.
- [PR #1101](https://github.com/openforcefield/openff-toolkit/pull/1101): Fixes a bug
  in which calling `to_qcschema` on a molecule with no connectivity feeds
  `QCElemental.Molecule` an empty list for the `connectivity` field; now feeds `None`.

### Tests updated
- [PR #1017](https://github.com/openforcefield/openforcefield/pull/1017): Ensures that OpenEye-only CI builds really
  do lack both AmberTools and RDKit.  

### Improved documentation and warnings
 - [PR #1065](https://github.com/openforcefield/openforcefield/pull/1017): Example notebooks were updated to use the
   Sage Open Force Field
 - [PR #1062](https://github.com/openforcefield/openforcefield/pull/1062): 
   Rewrote installation guide for clarity and comprehensiveness.

## 0.10.0 Improvements for force field fitting

### Behaviors changed

- [PR #1021](https://github.com/openforcefield/openforcefield/pull/1021): Renames
  [`openff.toolkit.utils.exceptions.ParseError`](openff.toolkit.utils.exceptions.ParseError) to
  [`openff.toolkit.utils.exceptions.SMILESParseError`](openff.toolkit.utils.exceptions.SMILESParseError) to
  avoid a conflict with an identically-named exception in the SMIRNOFF XML parsing code.
- [PR #1021](https://github.com/openforcefield/openforcefield/pull/1021): Renames and moves
  [`openff.toolkit.typing.engines.smirnoff.forcefield.ParseError`](openff.toolkit.typing.engines.smirnoff.forcefield.ParseError) to
  [`openff.toolkit.utils.exceptions.SMIRNOFFParseError`](openff.toolkit.utils.exceptions.SMIRNOFFParseError).
  This `ParseError` is deprecated and will be removed in a future release.

### New features and behaviors changed

- [PR #1027](https://github.com/openforcefield/openforcefield/pull/1027): Corrects interconversion of Molecule objects 
  with OEMol objects by ensuring atom names are correctly accessible via the `OEAtomBase.GetName()` and 
  `OEAtomBase.SetName()` methods, rather that the non-standard `OEAtomBase.GetData("name")` and 
  `OEAtomBase.SetData("name", name)`.
- [PR #1007](https://github.com/openforcefield/openforcefield/pull/1007): Resolves
  [Issue #456](https://github.com/openforcefield/openff-toolkit/issues/456) by adding the 
  `normalize_partial_charges` (default is `True`) keyword argument to 
  [`Molecule.assign_partial_charges`](openff.toolkit.topology.Molecule.assign_partial_charges),
  [`AmberToolsToolkitWrapper.assign_partial_charges`](openff.toolkit.utils.toolkits.AmberToolsToolkitWrapper.assign_partial_charges),  
  [`OpenEyeToolkitWrapper.assign_partial_charges`](openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.assign_partial_charges), 
  [`RDKitToolkitWrapper.assign_partial_charges`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.assign_partial_charges), and
  [`BuiltInToolkitWrapper.assign_partial_charges`](openff.toolkit.utils.toolkits.BuiltInToolkitWrapper.assign_partial_charges).
  This adds an offset to each atom's partial charge to ensure that their sum is equal to the net charge on the molecule
  (to the limit of a python float's precision, generally less than 1e-6 electron charge). **Note that, because this new 
  behavior is ON by default, it may slightly affect the partial charges and energies of systems generated by running
  [`create_openmm_system`](openff.toolkit.typing.engines.smirnoff.ForceField.create_openmm_system).**
- [PR #954](https://github.com/openforcefield/openforcefield/pull/954): Adds
  [`LibraryChargeType.from_molecule`](openff.toolkit.typing.engines.smirnoff.parameters.LibraryChargeHandler.LibraryChargeType.from_molecule)
  which returns a 
  [`LibraryChargeType`](openff.toolkit.typing.engines.smirnoff.parameters.LibraryChargeHandler.LibraryChargeType)
   object that will match the full molecule being parameterized, and assign
  it the same partial charges as are set on the input molecule.
- [PR #923](https://github.com/openforcefield/openforcefield/pull/923): Adds
  [`Molecule.nth_degree_neighbors`](openff.toolkit.topology.Molecule.nth_degree_neighbors),
  [`Topology.nth_degree_neighbors`](openff.toolkit.topology.Topology.nth_degree_neighbors),
  [`TopologyMolecule.nth_degree_neighbors`](openff.toolkit.topology.TopologyMolecule.nth_degree_neighbors),
  which returns pairs of atoms that are separated in a molecule or topology by _exactly_ N atoms.
- [PR #917](https://github.com/openforcefield/openforcefield/pull/917):
  [`ForceField.create_openmm_system`](openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system) 
  now ensures that the cutoff of the `NonbondedForce` is set to
  the cutoff of the `vdWHandler` when it and a `Electrostatics` handler are present in the force field.
- [PR #850](https://github.com/openforcefield/openforcefield/pull/850):
  [`OpenEyeToolkitWrapper.is_available`](openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.is_available)
  now returns `True` if _any_ OpenEye tools are licensed (and installed). This allows, i.e,
  use of functionality that requires `OEChem` without having an `OEOmega` license.
- [PR #909](https://github.com/openforcefield/openforcefield/pull/909): Virtual site positions can now
  be computed directly in the toolkit. This functionality is accessed through
  - [`FrozenMolecule.compute_virtual_site_positions_from_conformer`](openff.toolkit.topology.FrozenMolecule.compute_virtual_site_positions_from_conformer)
  - [`VirtualSite.compute_positions_from_conformer`](openff.toolkit.topology.VirtualSite.compute_positions_from_conformer)
  - [`VirtualParticle.compute_position_from_conformer`](openff.toolkit.topology.VirtualParticle.compute_position_from_conformer)
  - [`FrozenMolecule.compute_virtual_site_positions_from_atom_positions`](openff.toolkit.topology.FrozenMolecule.compute_virtual_site_positions_from_atom_positions)
  - [`VirtualSite.compute_positions_from_atom_positions`](openff.toolkit.topology.VirtualSite.compute_positions_from_atom_positions)
  - [`VirtualParticle.compute_position_from_atom_positions`](openff.toolkit.topology.VirtualParticle.compute_position_from_atom_positions)
    where the positions can be computed from a stored conformer, or an input vector of atom positions.
  - Tests have been added (`TestMolecule.test_*_virtual_site_position`) to check for sane behavior. The tests do
    not directly compare OpenMM position equivalence, but offline tests show that they are equivalent.
  - The helper method 
    [`VirtualSiteHandler.create_openff_virtual_sites`](openff.toolkit.typing.engines.smirnoff.parameters.VirtualSiteHandler.create_openff_virtual_sites) 
    is now public, which returns a modified topology with virtual sites added.
  - Virtual sites now expose the parameters used to create its local frame via the read-only properties
    - [`VirtualSite.local_frame_weights`](openff.toolkit.topology.VirtualSite.local_frame_weights)
    - [`VirtualSite.local_frame_position`](openff.toolkit.topology.VirtualSite.local_frame_position)
  - Adding virtual sites via the `Molecule` API now have defaults for `sigma`, `epsilon`, and `charge_increment`
    set to 0 with appropriate units, rather than `None`
- [PR #956](https://github.com/openforcefield/openforcefield/pull/956): Added 
  [`ForceField.get_partial_charges()`](openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.get_partial_charges)
  to more easily compute the partial charges assigned by a force field for a molecule.
- [PR  #1006](https://github.com/openforcefield/openff-toolkit/pull/1006):
  Two behavior changes in the SMILES output for `to_file()` and `to_file_obj()`:
  - The RDKit and OpenEye wrappers now output the same SMILES as `to_smiles()`.
   This uses explicit hydrogens rather than the toolkit's default of implicit hydrogens.
  - The RDKit wrapper no longer includes a header line. This improves
  the consistency between the OpenEye and RDKit outputs.

### Bugfixes

- [PR #1024](https://github.com/openforcefield/openforcefield/pull/1024): Small changes
  for compatibility with OpenMM 7.6.
- [PR #1003](https://github.com/openforcefield/openforcefield/pull/1003): Fixes
  [Issue #1000](https://github.com/openforcefield/openff-toolkit/issues/1000), where a stereochemistry
  warning is sometimes erroneously emitted when loading a stereogenic molecule using
  [`Molecule.from_pdb_and_smiles`](openff.toolkit.topology.Molecule.from_pdb_and_smiles)
- [PR #1002](https://github.com/openforcefield/openforcefield/pull/1002): Fixes a bug in which OFFXML files could
  inadvertently be loaded from subdirectories.
- [PR #969](https://github.com/openforcefield/openforcefield/pull/969): Fixes a bug in which the cutoff distance
  of the `NonbondedForce` generated by
  [`ForceField.create_openmm_system`](openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system)
  was not set to the value specified by the vdW and Electrostatics handlers.
- [PR #909](https://github.com/openforcefield/openforcefield/pull/909): Fixed several bugs related to creating an
  OpenMM system with virtual sites created via the `Molecule` virtual site API
- [PR  #1006](https://github.com/openforcefield/openff-toolkit/pull/1006):
  Many small fixes to the toolkit wrapper I/O for better error
  handling, improved consistency between reading from a file vs. file
  object, and improved consistency between the RDKit and OEChem
  toolkit wrappers. For the full list see
  [Issue #1005](https://github.com/openforcefield/openff-toolkit/issues/1005). Some
  of the more significant fixes are:
  - [`RDKitToolkitWrapper.from_file_obj()`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.from_file_obj) now uses the same
    structure normaliation as `from_file()`.
  - `from_smiles()` now raises an `openff.toolkit.utils.exceptions.SMILESParsingError` if
  the SMILES could not be parsed.
  - OEChem input and output files now raise an OSError if the file
  could not be opened.
  - All input file object readers now support file objects open in binary mode.

### Examples added

- [PR #763](https://github.com/openforcefield/openff-toolkit/pull/763):
  Adds an introductory example showcasing the toolkit parameterizing a protein-ligand simulation.
- [PR #955](https://github.com/openforcefield/openff-toolkit/pull/955): Refreshed the force field modification example
- [PR #934](https://github.com/openforcefield/openff-toolkit/pull/934)
  and [conda-forge/openff-toolkit-feedstock#9](https://github.com/conda-forge/openff-toolkit-feedstock/pull/9):
  Added `openff-toolkit-examples` Conda package for easy installation of examples and their
  dependencies. Simply `conda install -c conda-forge openff-toolkit-examples` and then run
  the `openff-toolkit-examples` script to copy the examples suite to a convenient place to
  run them!

### Tests updated

- [PR #963](https://github.com/openforcefield/openff-toolkit/pull/963):
  Several tests modules used functions from test_forcefield.py that created an OpenFF Molecule
  without a toolkit. These functions are now in their own module so they can be imported directly,
  without the overhead of going through test_forcefield.
- [PR #997](https://github.com/openforcefield/openff-toolkit/pull/997):
  Several XML snippets in `test_forcefield.py` that were scattered around inside of classes and
  functions are now moved to the module level.
  
  
## 0.9.2 Minor feature and bugfix release

### New features and behaviors changed

- [PR #762](https://github.com/openforcefield/openforcefield/pull/762):
  [`Molecule.from_rdkit`](openff.toolkit.topology.Molecule.from_rdkit) now converts
  implicit hydrogens into explicit hydrogens by default. This change may affect
  [`RDKitToolkitWrapper/Molecule.from_smiles`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.from_smiles),
  [`from_mapped_smiles`](openff.toolkit.topology.Molecule.from_mapped_smiles),
  [`from_file`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.from_file),
  [`from_file_obj`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.from_file_obj),
  [`from_inchi`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.from_inchi), and
  [`from_qcschema`](openff.toolkit.topology.Molecule.from_qcschema).
  This new behavior can be disabled using the
  `hydrogens_are_explicit=True` keyword argument to
  [`from_smiles`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.from_smiles),
  or loading the molecule into the desired explicit protonation state in RDKit, then calling
  [`from_rdkit`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.from_rdkit) on the RDKit molecule with
  `hydrogens_are_explicit=True`.
- [PR #894](https://github.com/openforcefield/openforcefield/pull/894): Calls to
  [`Molecule.from_openeye`](openff.toolkit.topology.Molecule.from_openeye),
  [`Molecule.from_rdkit`](openff.toolkit.topology.Molecule.from_rdkit),
  [`Molecule.from_smiles`](openff.toolkit.topology.Molecule.from_smiles),
  [`OpenEyeToolkitWrapper.from_smiles`](openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.from_smiles), and
  [`RDKitToolkitWrapper.from_smiles`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.from_smiles)
  will now load atom maps into the the resulting
  `Molecule's` `offmol.properties['atom_map']` field, even if not all atoms have map indices assigned.
- [PR #904](https://github.com/openforcefield/openforcefield/pull/904):
  [`TopologyAtom`](openff.toolkit.topology.TopologyAtom.element) objects now have
  an element getter [`TopologyAtom.element`](openff.toolkit.topology.TopologyAtom.element).

### Bugfixes

- [PR #891](https://github.com/openforcefield/openforcefield/pull/891): Calls to
  [`Molecule/OpenEyeToolkitWrapper.from_openeye`](openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.from_openeye)
  no longer mutate the input OE molecule.
- [PR #897](https://github.com/openforcefield/openforcefield/pull/897): Fixes enumeration of stereoisomers for
  molecules with already defined stereochemistry using
  [`RDKitToolkitWrapper.enumerate_stereoisomers`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.enumerate_stereoisomers).
- [PR #859](https://github.com/openforcefield/openforcefield/pull/859): Makes
  [`RDKitToolkitWrapper.enumerate_tautomers`](openff.toolkit.utils.toolkits.RDKitToolkitWrapper.enumerate_tautomers)
  actually use the `max_states` keyword argument during tautomer generation, which will reduce resource use in some
  cases.

### Improved documentation and warnings
- [PR #862](https://github.com/openforcefield/openforcefield/pull/862): Clarify that `System` objects produced by the
  toolkit are OpenMM `System`s in anticipation of forthcoming OpenFF `System`s. Fixes
  [Issue #618](https://github.com/openforcefield/openforcefield/issues/618).
- [PR #863](https://github.com/openforcefield/openff-toolkit/pull/863): Documented how to build the docs in the
  developers guide.
- [PR #870](https://github.com/openforcefield/openff-toolkit/pull/870): Reorganised documentation to improve
  discoverability and allow future additions.
- [PR #871](https://github.com/openforcefield/openff-toolkit/pull/871): Changed Markdown parser from m2r2 to MyST for
  improved documentation rendering.
- [PR #880](https://github.com/openforcefield/openff-toolkit/pull/880): Cleanup and partial rewrite of the developer's
  guide.
- [PR #906](https://github.com/openforcefield/openff-toolkit/pull/906): Cleaner instructions on how to setup
  development environment.

:::{TODO}
- Translate previous release history to MyST markdown
:::

## Earlier releases

:::{eval-rst}

0.9.1 - Minor feature and bugfix release
----------------------------------------

New features
""""""""""""
- `PR #839 <https://github.com/openforcefield/openforcefield/pull/839>`_: Add support for computing WBOs from multiple
  conformers using the AmberTools and OpenEye toolkits, and from ELF10 conformers using the OpenEye toolkit wrapper.
- `PR #832 <https://github.com/openforcefield/openforcefield/pull/832>`_: Expose ELF conformer selection through the
  ``Molecule`` API via a new ``apply_elf_conformer_selection`` function.
- `PR #831 <https://github.com/openforcefield/openff-toolkit/pull/831>`_: Expose ELF conformer selection through the
  OpenEye wrapper.
- `PR #790 <https://github.com/openforcefield/openforcefield/pull/790>`_: Fixes `Issue #720
  <https://github.com/openforcefield/openforcefield/issues/720>`_ where qcschema roundtrip to/from results
  in an error due to missing cmiles entry in attributes.
- `PR #793 <https://github.com/openforcefield/openff-toolkit/pull/793>`_: Add an initial ELF conformer selection
  implementation which uses RDKit.
- `PR #799 <https://github.com/openforcefield/openff-toolkit/pull/799>`_: Closes
  `Issue #746 <https://github.com/openforcefield/openff-toolkit/issues/746>`_ by adding
  :py:meth:`Molecule.smirnoff_impropers <openff.toolkit.topology.FrozenMolecule.smirnoff_impropers>`,
  :py:meth:`Molecule.amber_impropers <openff.toolkit.topology.FrozenMolecule.amber_impropers>`,
  :py:meth:`TopologyMolecule.smirnoff_impropers <openff.toolkit.topology.TopologyMolecule.smirnoff_impropers>`,
  :py:meth:`TopologyMolecule.amber_impropers <openff.toolkit.topology.TopologyMolecule.amber_impropers>`,
  :py:meth:`Topology.smirnoff_impropers <openff.toolkit.topology.Topology.smirnoff_impropers>`, and
  :py:meth:`Topology.amber_impropers <openff.toolkit.topology.Topology.amber_impropers>`.
- `PR #847 <https://github.com/openforcefield/openforcefield/pull/847>`_: Instances of
  :py:class:`ParameterAttribute <openff.toolkit.typing.engines.smirnoff.parameters.ParameterAttribute>`
  documentation can now specify their docstrings with the optional ``docstring`` argument to the
  ``__init__()`` method.
- `PR #827 <https://github.com/openforcefield/openff-toolkit/pull/827>`_: The
  setter for :py:class:`Topology.box_vectors <openff.toolkit.topology.Topology>` now infers box vectors
  when box lengths are pass as a list of length 3.

Behavior changed
""""""""""""""""
- `PR #802 <https://github.com/openforcefield/openforcefield/pull/802>`_: Fixes
  `Issue #408 <https://github.com/openforcefield/openforcefield/issues/408>`_. The 1-4 scaling
  factor for electrostatic interactions is now properly set by the value specified in the force
  field. Previously it fell back to a default value of 0.83333. The toolkit may now produce
  slightly different energies as a result of this change.
- `PR #839 <https://github.com/openforcefield/openforcefield/pull/839>`_: The average WBO will now be returned when
  multiple conformers are provided to ``assign_fractional_bond_orders`` using ``use_conformers``.
- `PR #816 <https://github.com/openforcefield/openforcefield/pull/816>`_: Force field file paths
  are now loaded in a case-insensitive manner.

Bugfixes
""""""""
- `PR #849 <https://github.com/openforcefield/openforcefield/pull/849>`_: Changes
  :py:meth:`create_openmm_system <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system>` so
  that it no longer uses the conformers on existing reference molecules (if present) to calculate Wiberg
  bond orders. Instead, new conformers are always generated during parameterization.

Improved documentation and warnings
"""""""""""""""""""""""""""""""""""
- `PR #838 <https://github.com/openforcefield/openforcefield/pull/838>`_: Corrects spacing of "forcefield" to "force
  field" throughout documentation. Fixes `Issue #112 <https://github.com/openforcefield/openforcefield/issues/112>`_.
- `PR #846 <https://github.com/openforcefield/openff-toolkit/pull/846>`_: Corrects dead links throughout release history.
  Fixes `Issue #835 <https://github.com/openforcefield/openff-toolkit/issues/835>`_.
- `PR #847 <https://github.com/openforcefield/openforcefield/pull/847>`_: Documentation now compiles
  with far fewer warnings, and in many cases more correctly. Additionally, :py:class:`ParameterAttribute
  <openff.toolkit.typing.engines.smirnoff.parameters.ParameterAttribute>` documentation no longer
  appears incorrectly in classes where it is used. Fixes `Issue #397
  <https://github.com/openforcefield/openforcefield/issues/397>`_.

0.9.0 - Namespace Migration
---------------------------

This release marks the transition from the old ``openforcefield`` branding over to its new
identity as ``openff-toolkit``. This change has been made to better represent the role of the
toolkit, and highlight its place in the larger Open Force Field (OpenFF) ecosystem.

From version ``0.9.0`` onwards the toolkit will need to be imported as ``import openff.toolkit.XXX`` and
``from openff.toolkit import XXX``.

API-breaking changes
""""""""""""""""""""
- `PR #803 <https://github.com/openforcefield/openff-toolkit/pull/803>`_: Migrates ``openforcefield``
  imports to ``openff.toolkit``.


0.8.4 - Minor feature and bugfix release
----------------------------------------

**This release is intended to be functionally identical to 0.9.1.
The only difference is that it uses the "openforcefield" namespace.**

This release is a final patch for the ``0.8.X`` series of releases of the toolkit, and also marks the last
version of the toolkit which will be imported as ``import openforcefield.XXX`` / ``from openforcefield import XXX``.
From version ``0.9.0`` onwards the toolkit will be importable only as ``import openff.toolkit.XXX`` /
``from openff.toolkit import XXX``.

**Note** This change will also be accompanied by a renaming of the package from ``openforcefield`` to ``openff-toolkit``,
so users need not worry about accidentally pulling in a version with changed imports. Users will have to explicitly
choose to install the ``openff-toolkit`` package once released which will contain the breaking import changes.


0.8.3 - Major bugfix release
----------------------------

This release fixes a critical bug in van der Waals parameter assignment.

This release is also a final patch for the ``0.8.X`` series of releases of the toolkit, and also marks the last
version of the toolkit which will be imported as ``import openforcefield.XXX`` / ``from openforcefield import XXX``.
From version ``0.9.0`` onwards the toolkit will be importable only as ``import openff.toolkit.XXX`` /
``from openff.toolkit import XXX``.

**Note** This change will also be accompanied by a renaming of the package from ``openforcefield`` to ``openff-toolkit``,
so users need not worry about accidentally pulling in a version with changed imports. Users will have to explicitly
choose to install the ``openff-toolkit`` package once released which will contain the breaking import changes.

Bugfixes
""""""""
- `PR #808 <https://github.com/openforcefield/openff-toolkit/pull/808>`_: Fixes
  `Issue #807 <https://github.com/openforcefield/openff-toolkit/issues/807>`_,
  which tracks a major bug in the interconversion between a vdW ``sigma``
  and ``rmin_half`` parameter.


New features
""""""""""""
- `PR #794 <https://github.com/openforcefield/openff-toolkit/pull/794>`_: Adds a decorator
  ``@requires_package`` that denotes a function requires an optional dependency.
- `PR #805 <https://github.com/openforcefield/openff-toolkit/pull/805>`_: Adds a deprecation warning for the up-coming
  release of the ``openff-toolkit`` package and its import breaking changes.

0.8.2 - Bugfix release
----------------------

**WARNING: This release was later found to contain a major bug,**
`Issue #807 <https://github.com/openforcefield/openff-toolkit/issues/807>`_,
**and produces incorrect energies.**

Bugfixes
""""""""
- `PR #786 <https://github.com/openforcefield/openff-toolkit/pull/786>`_: Fixes `Issue #785
  <https://github.com/openforcefield/openff-toolkit/issues/785>`_ where RDKitToolkitWrapper would
  sometimes expect stereochemistry to be defined for non-stereogenic bonds when loading from
  SDF.
- `PR #786 <https://github.com/openforcefield/openff-toolkit/pull/786>`_: Fixes an issue where
  using the :py:class:`Molecule <openff.toolkit.topology.Molecule>` copy constructor
  (``newmol = Molecule(oldmol)``) would result
  in the copy sharing the same ``.properties`` dict as the original (as in, changes to the
  ``.properties`` dict of the copy would be reflected in the original).
- `PR #789 <https://github.com/openforcefield/openff-toolkit/pull/789>`_: Fixes a regression noted in
  `Issue #788 <https://github.com/openforcefield/openff-toolkit/issues/788>`_
  where creating
  :py:class:`vdWHandler.vdWType <openff.toolkit.typing.engines.smirnoff.parameters.vdWHandler.vdWType>`
  or setting ``sigma`` or ``rmin_half`` using Quantities represented as strings resulted in an error.


0.8.1 - Bugfix and minor feature release
----------------------------------------

**WARNING: This release was later found to contain a major bug,**
`Issue #807 <https://github.com/openforcefield/openff-toolkit/issues/807>`_,
**and produces incorrect energies.**

API-breaking changes
""""""""""""""""""""
- `PR #757 <https://github.com/openforcefield/openff-toolkit/pull/757>`_: Renames
  ``test_forcefields/smirnoff99Frosst.offxml`` to ``test_forcefields/test_forcefield.offxml``
  to avoid confusion with any of the ACTUAL released FFs in the
  `smirnoff99Frosst line <https://github.com/openforcefield/smirnoff99Frosst/>`_
- `PR #751 <https://github.com/openforcefield/openff-toolkit/pull/751>`_: Removes the
  optional ``oetools=("oechem", "oequacpac", "oeiupac", "oeomega")`` keyword argument from
  :py:meth:`OpenEyeToolkitWrapper.is_available <openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.is_available>`, as
  there are no special behaviors that are accessed in the case of partially-licensed OpenEye backends. The
  new behavior of this method is the same as if the default value above is always provided.

Behavior Changed
""""""""""""""""
- `PR #583 <https://github.com/openforcefield/openff-toolkit/pull/583>`_: Methods
  such as :py:meth:`Molecule.from_rdkit <openff.toolkit.topology.Molecule.from_rdkit>`
  and :py:meth:`Molecule.from_openeye <openff.toolkit.topology.Molecule.from_openeye>`,
  which delegate their internal logic to :py:class:`ToolkitRegistry <openff.toolkit.utils.toolkits.ToolkitRegistry>`
  functions, now guarantee that they will return an object of the correct type when being called on ``Molecule``-derived classes. Previously,
  running these constructors using subclasses of :py:class:`FrozenMolecule <openff.toolkit.topology.Molecule>`
  would not return an instance of that subclass, but rather just an instance of a
  :py:class:`Molecule <openff.toolkit.topology.Molecule>`.
- `PR #753 <https://github.com/openforcefield/openff-toolkit/pull/753>`_: ``ParameterLookupError``
  is now raised when passing to
  :py:meth:`ParameterList.index <openff.toolkit.typing.engines.smirnoff.parameters.ParameterList>`
  a SMIRKS pattern not found in the parameter list.

New features
""""""""""""
- `PR #751 <https://github.com/openforcefield/openff-toolkit/pull/751>`_: Adds
  ``LicenseError``, a subclass of ``ToolkitUnavailableException`` which is raised when attempting to
  add a cheminformatics :py:class:`ToolkitWrapper <openff.toolkit.utils.toolkits.ToolkitWrapper>` for
  a toolkit that is installed but unlicensed.
- `PR #678 <https://github.com/openforcefield/openff-toolkit/pull/678>`_: Adds
  :py:meth:`ForceField.deregister_parameter_handler <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.deregister_parameter_handler>`.
- `PR #730 <https://github.com/openforcefield/openff-toolkit/pull/730>`_: Adds
  :py:class:`Topology.is_periodic <openff.toolkit.topology.Topology>`.
- `PR #753 <https://github.com/openforcefield/openff-toolkit/pull/753>`_: Adds
  :py:meth:`ParameterHandler.__getitem__ <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>`
  to look up individual :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>`
  objects.

Bugfixes
""""""""
- `PR #745 <https://github.com/openforcefield/openff-toolkit/pull/745>`_: Fixes bug when
  serializing molecule with conformers to JSON.
- `PR #750 <https://github.com/openforcefield/openff-toolkit/pull/750>`_: Fixes a bug causing either
  ``sigma`` or ``rmin_half`` to sometimes be missing on
  :py:class:`vdWHandler.vdWType <openff.toolkit.typing.engines.smirnoff.parameters.vdWHandler.vdWType>`
  objects.
- `PR #756 <https://github.com/openforcefield/openff-toolkit/pull/756>`_: Fixes bug when running
  :py:meth:`vdWHandler.create_force <openff.toolkit.typing.engines.smirnoff.parameters.vdWHandler>`
  using a ``vdWHandler`` that was initialized using the API.
- `PR #776 <https://github.com/openforcefield/openff-toolkit/pull/776>`_: Fixes a bug in which
  the :py:meth:`Topology.from_openmm <openff.toolkit.topology.Topology.from_openmm>` and
  :py:meth:`Topology.from_mdtraj <openff.toolkit.topology.Topology.from_mdtraj>` methods would
  dangerously allow ``unique_molecules=None``.
- `PR #777 <https://github.com/openforcefield/openff-toolkit/pull/777>`_:
  :py:class:`RDKitToolkitWrapper <openff.toolkit.utils.toolkits.RDKitToolkitWrapper>`
  now outputs the full warning message when ``allow_undefined_stereo=True`` (previously the
  description of which stereo was undefined was squelched)


0.8.0 - Virtual Sites
---------------------

**Major Feature: Support for the SMIRNOFF VirtualSite tag**

This release implements the SMIRNOFF virtual site specification. The implementation enables support
for models using off-site charges, including 4- and 5-point water models, in addition to lone pair
modeling on various functional groups. The primary focus was on the ability to parameterize a
system using virtual sites, and generating an OpenMM system with all virtual sites present and
ready for evaluation. Support for formats other than OpenMM has not be implemented in this release,
but may come with the appearance of the OpenFF system object. In addition to implementing the
specification, the toolkit :py:class:`Molecule <openff.toolkit.topology.Molecule>` objects now
allow the creation and manipulation of virtual sites.

This change is documented in the `Virtual sites page <virtualsites.html>`_ of the user guide.


**Minor Feature: Support for the 0.4 ChargeIncrementModel tag**

To allow for more convenient fitting of ``ChargeIncrement`` parameters, it is now possible to specify one less
``charge_increment`` value than there are tagged atoms in a ``ChargeIncrement``'s ``smirks``. The missing
``charge_increment`` value will be calculated at parameterization-time to make the sum of
the charge contributions from a ``ChargeIncrement`` parameter equal to zero.
Since this change allows for force fields that are incompatible with
the previous specification, this new style of ``ChargeIncrement`` must specify a ``ChargeIncrementModel``
section version of ``0.4``. All ``0.3``-compatible ``ChargeIncrement`` parameters are compatible with
the ``0.4`` ``ChargeIncrementModel`` specification.

More details and examples of this change are available in `The ChargeIncrementModel tag in the SMIRNOFF specification <https://openforcefield.github.io/standards/standards/smirnoff/#chargeincrementmodel-small-molecule-and-fragment-charges>`_


New features
""""""""""""
- `PR #726 <https://github.com/openforcefield/openff-toolkit/pull/726>`_: Adds support for the 0.4
  ChargeIncrementModel spec, allowing for the specification of one fewer ``charge_increment`` values
  than there are tagged atoms in the ``smirks``, and automatically assigning the final atom an offsetting charge.
- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Adds support for the ``VirtualSites`` tag in the SMIRNOFF specification

- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Adds ``replace`` and ``all_permutations`` kwarg to

  - :py:meth:`Molecule.add_bond_charge_virtual_site <openff.toolkit.topology.Molecule.add_bond_charge_virtual_site>`
  - :py:meth:`Molecule.add_monovalent_lone_pair_virtual_site <openff.toolkit.topology.Molecule.add_monovalent_lone_pair_virtual_site>`
  - :py:meth:`Molecule.add_divalent_lone_pair_virtual_site <openff.toolkit.topology.Molecule.add_divalent_lone_pair_virtual_site>`
  - :py:meth:`Molecule.add_trivalent_lone_pair_virtual_site <openff.toolkit.topology.Molecule.add_trivalent_lone_pair_virtual_site>`

- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Adds ``orientations`` to

  - :py:class:`BondChargeVirtualSite <openff.toolkit.topology.BondChargeVirtualSite>`
  - :py:class:`MonovalentLonePairVirtualSite <openff.toolkit.topology.MonovalentLonePairVirtualSite>`
  - :py:class:`DivalentLonePairVirtualSite <openff.toolkit.topology.DivalentLonePairVirtualSite>`
  - :py:class:`TrivalentLonePairVirtualSite <openff.toolkit.topology.TrivalentLonePairVirtualSite>`

- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Adds

  - :py:class:`VirtualParticle <openff.toolkit.topology.VirtualParticle>`
  - :py:class:`TopologyVirtualParticle <openff.toolkit.topology.TopologyVirtualParticle>`
  - :py:meth:`BondChargeVirtualSite.get_openmm_virtual_site <openff.toolkit.topology.BondChargeVirtualSite.get_openmm_virtual_site>`
  - :py:meth:`MonovalentLonePairVirtualSite.get_openmm_virtual_site <openff.toolkit.topology.MonovalentLonePairVirtualSite.get_openmm_virtual_site>`
  - :py:meth:`DivalentLonePairVirtualSite.get_openmm_virtual_site <openff.toolkit.topology.DivalentLonePairVirtualSite.get_openmm_virtual_site>`
  - :py:meth:`TrivalentLonePairVirtualSite.get_openmm_virtual_site <openff.toolkit.topology.TrivalentLonePairVirtualSite.get_openmm_virtual_site>`
  - :py:meth:`ValenceDict.key_transform <openff.toolkit.topology.ValenceDict.key_transform>`
  - :py:meth:`ValenceDict.index_of <openff.toolkit.topology.ValenceDict.index_of>`
  - :py:meth:`ImproperDict.key_transform <openff.toolkit.topology.ImproperDict.key_transform>`
  - :py:meth:`ImproperDict.index_of <openff.toolkit.topology.ImproperDict.index_of>`

- `PR #705 <https://github.com/openforcefield/openff-toolkit/pull/705>`_: Adds interpolation
  based on fractional bond orders for harmonic bonds. This includes interpolation for both
  the force constant ``k`` and/or equilibrium bond distance ``length``. This is accompanied by a
  bump in the ``<Bonds>`` section of the SMIRNOFF spec (but not the entire spec).
- `PR #718 <https://github.com/openforcefield/openff-toolkit/pull/718>`_: Adds ``.rings`` and
  ``.n_rings`` to :py:class:`Molecule <openff.toolkit.topology.Molecule>` and ``.is_in_ring``
  to :py:class:`Atom <openff.toolkit.topology.Atom>` and
  :py:class:`Bond <openff.toolkit.topology.Bond>`

Bugfixes
"""""""""
- `PR #682 <https://github.com/openforcefield/openff-toolkit/pull/682>`_: Catches failures in
  :py:meth:`Molecule.from_iupac <openff.toolkit.topology.Molecule.from_iupac>` instead of silently
  failing.
- `PR #743 <https://github.com/openforcefield/openff-toolkit/pull/743>`_: Prevents the non-bonded
  (vdW) cutoff from silently falling back to the OpenMM default of 1 nm in
  :py:meth:`Forcefield.create_openmm_system
  <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system>` and instead
  sets its to the value specified by the force field.
- `PR #737 <https://github.com/openforcefield/openff-toolkit/pull/737>`_: Prevents OpenEye from
  incidentally being used in the conformer generation step of
  :py:class:`AmberToolsToolkitWrapper.assign_fractional_bond_orders
  <openff.toolkit.utils.toolkits.AmberToolsToolkitWrapper.assign_fractional_bond_orders>`.

Behavior changed
""""""""""""""""
- `PR #705 <https://github.com/openforcefield/openff-toolkit/pull/705>`_: Changes the default values
  in the ``<Bonds>`` section of the SMIRNOFF spec to ``fractional_bondorder_method="AM1-Wiberg"``
  and ``potential="(k/2)*(r-length)^2"``, which is backwards-compatible with and equivalent to
  ``potential="harmonic"``.

Examples added
""""""""""""""
- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Adds a virtual site example notebook to run
  an OpenMM simulation with virtual sites, and compares positions and potential energy of TIP5P water between OpenFF
  and OpenMM force fields.

API-breaking changes
""""""""""""""""""""
- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Methods

  - :py:meth:`Molecule.add_bond_charge_virtual_site <openff.toolkit.topology.Molecule.add_bond_charge_virtual_site>`
  - :py:meth:`Molecule.add_monovalent_lone_pair_virtual_site <openff.toolkit.topology.Molecule.add_monovalent_lone_pair_virtual_site>`
  - :py:meth:`Molecule.add_divalent_lone_pair_virtual_site <openff.toolkit.topology.Molecule.add_divalent_lone_pair_virtual_site>`
  - :py:meth:`Molecule.add_trivalent_lone_pair_virtual_site <openff.toolkit.topology.Molecule.add_trivalent_lone_pair_virtual_site>`
    now only accept a list of atoms, not a list of integers, to define to parent atoms

- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Removes
  :py:meth:`VirtualParticle.molecule_particle_index <openff.toolkit.topology.VirtualParticle.molecule_particle_index>`

- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Removes ``outOfPlaneAngle`` from

  - :py:class:`DivalentLonePairVirtualSite <openff.toolkit.topology.DivalentLonePairVirtualSite>`
  - :py:class:`TrivalentLonePairVirtualSite <openff.toolkit.topology.TrivalentLonePairVirtualSite>`

- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Removes ``inPlaneAngle`` from
  :py:class:`TrivalentLonePairVirtualSite <openff.toolkit.topology.TrivalentLonePairVirtualSite>`

- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Removes ``weights`` from

  - :py:class:`BondChargeVirtualSite <openff.toolkit.topology.BondChargeVirtualSite>`
  - :py:class:`MonovalentLonePairVirtualSite <openff.toolkit.topology.MonovalentLonePairVirtualSite>`
  - :py:class:`DivalentLonePairVirtualSite <openff.toolkit.topology.DivalentLonePairVirtualSite>`
  - :py:class:`TrivalentLonePairVirtualSite <openff.toolkit.topology.TrivalentLonePairVirtualSite>`

Tests added
"""""""""""

- `PR #548 <https://github.com/openforcefield/openff-toolkit/pull/548>`_: Adds test for

  - The virtual site parameter handler
  - TIP5P water dimer energy and positions
  - Adds tests to for virtual site/particle indexing/counting


0.7.2 - Bugfix and minor feature release
----------------------------------------

New features
""""""""""""
- `PR #662 <https://github.com/openforcefield/openff-toolkit/pull/662>`_: Adds ``.aromaticity_model``
  of :py:class:`ForceField <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>` and ``.TAGNAME``
  of :py:class:`ParameterHandler <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>` as
  public attributes.
- `PR #667 <https://github.com/openforcefield/openff-toolkit/pull/667>`_ and
  `PR #681 <https://github.com/openforcefield/openff-toolkit/pull/681>`_ linted the codebase with
  ``black`` and ``isort``, respectively.
- `PR #675 <https://github.com/openforcefield/openff-toolkit/pull/675>`_ adds
  ``.toolkit_version`` to
  :py:class:`ToolkitWrapper <openff.toolkit.utils.toolkits.ToolkitWrapper>` and
  ``.registered_toolkit_versions`` to
  :py:class:`ToolkitRegistry <openff.toolkit.utils.toolkits.ToolkitRegistry>`.
- `PR #696 <https://github.com/openforcefield/openff-toolkit/pull/696>`_ Exposes a setter for
  :py:class:`ForceField.aromaticity_model <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>`
- `PR #685 <https://github.com/openforcefield/openff-toolkit/pull/685>`_ Adds a custom ``__hash__``
  function to
  :py:class:`ForceField <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>`


Behavior changed
""""""""""""""""
- `PR #684 <https://github.com/openforcefield/openff-toolkit/pull/684>`_: Changes
  :py:class:`ToolkitRegistry <openff.toolkit.utils.toolkits.ToolkitRegistry>` to return an empty
  registry when initialized with no arguments, i.e. ``ToolkitRegistry()`` and makes the
  ``register_imported_toolkit_wrappers`` argument private.
- `PR #711 <https://github.com/openforcefield/openff-toolkit/pull/711>`_: The
  setter for :py:class:`Topology.box_vectors <openff.toolkit.topology.Topology>`
  now infers box vectors (a 3x3 matrix) when box lengths
  (a 3x1 array) are passed, assuming an orthogonal box.
- `PR #649 <https://github.com/openforcefield/openff-toolkit/pull/648>`_: Makes SMARTS
  searches stereochemistry-specific (if stereo is specified in the SMARTS) for both OpenEye
  and RDKit backends. Also ensures molecule
  aromaticity is re-perceived according to the ForceField's specified
  aromaticity model, which may overwrite user-specified aromaticity on the ``Molecule``
- `PR #648 <https://github.com/openforcefield/openff-toolkit/pull/648>`_: Removes the
  ``utils.structure`` module, which was deprecated in 0.2.0.
- `PR #670 <https://github.com/openforcefield/openff-toolkit/pull/670>`_: Makes the
  :py:class:`Topology <openff.toolkit.topology.Topology>` returned by ``create_openmm_system``
  contain the partial charges and partial bond orders (if any) assigned during parameterization.
- `PR #675 <https://github.com/openforcefield/openff-toolkit/pull/675>`_ changes the
  exception raised when no ``antechamber`` executable is found from ``IOError`` to
  ``AntechamberNotFoundError``
- `PR #696 <https://github.com/openforcefield/openff-toolkit/pull/696>`_ Adds an
  ``aromaticity_model`` keyword argument to the
  :py:class:`ForceField <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>`
  constructor, which defaults to ``DEFAULT_AROMATICITY_MODEL``.

Bugfixes
"""""""""
- `PR #715 <https://github.com/openforcefield/openff-toolkit/pull/715>`_: Closes issue `Issue #475
  <https://github.com/openforcefield/openff-toolkit/issues/475>`_ writing a "PDB" file using OE backend rearranges
  the order of the atoms by pushing the hydrogens to the bottom.
- `PR #649 <https://github.com/openforcefield/openff-toolkit/pull/648>`_: Prevents 2020 OE
  toolkit from issuing a warning caused by doing stereo-specific smarts searches on certain
  structures.
- `PR #724 <https://github.com/openforcefield/openff-toolkit/pull/724>`_: Closes issue `Issue #502
  <https://github.com/openforcefield/openff-toolkit/issues/502>`_ Adding a utility function Topology.to_file() to
  write topology and positions to a "PDB" file using openmm backend for pdb file write.

Tests added
"""""""""""
- `PR #694 <https://github.com/openforcefield/openff-toolkit/pull/694>`_: Adds automated testing
  to code snippets in docs.
- `PR #715 <https://github.com/openforcefield/openff-toolkit/pull/715>`_: Adds tests for pdb file writes using OE
  backend.
- `PR #724 <https://github.com/openforcefield/openff-toolkit/pull/724>`_: Adds tests for the utility function Topology.to_file().


0.7.1 - OETK2020 Compatibility and Minor Update
-----------------------------------------------

This is the first of our patch releases on our new planned monthly release schedule.

Detailed release notes are below, but the major new features of this release are updates for
compatibility with the new 2020 OpenEye Toolkits release, the
``get_available_force_fields`` function, and the disregarding of pyrimidal nitrogen stereochemistry
in molecule isomorphism checks.

Behavior changed
""""""""""""""""
- `PR #646 <https://github.com/openforcefield/openff-toolkit/pull/646>`_: Checking for
  :py:class:`Molecule <openff.toolkit.topology.Molecule>`
  equality using the ``==`` operator now disregards all pyrimidal nitrogen stereochemistry
  by default. To re-enable, use
  :py:class:`Molecule.{is|are}_isomorphic <openff.toolkit.topology.Molecule>`
  with the ``strip_pyrimidal_n_atom_stereo=False`` keyword argument.
- `PR #646 <https://github.com/openforcefield/openff-toolkit/pull/646>`_: Adds
  an optional ``toolkit_registry`` keyword argument to
  :py:class:`Molecule.are_isomorphic <openff.toolkit.topology.Molecule>`,
  which identifies the toolkit that should be used to search for pyrimidal nitrogens.


Bugfixes
""""""""
- `PR #647 <https://github.com/openforcefield/openff-toolkit/pull/647>`_: Updates
  :py:class:`OpenEyeToolkitWrapper <openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper>`
  for 2020.0.4 OpenEye Toolkit behavior/API changes.
- `PR #646 <https://github.com/openforcefield/openff-toolkit/pull/646>`_: Fixes a bug where
  :py:class:`Molecule.chemical_environment_matches <openff.toolkit.topology.Molecule>`
  was not able to accept a :py:class:`ChemicalEnvironment <openff.toolkit.typing.chemistry.ChemicalEnvironment>` object
  as a query.
- `PR #634 <https://github.com/openforcefield/openff-toolkit/pull/634>`_: Fixes a bug in which calling
  :py:class:`RDKitToolkitWrapper.from_file <openff.toolkit.utils.toolkits.RDKitToolkitWrapper>` directly
  would not load files correctly if passed lowercase ``file_format``. Note that this bug did not occur when calling
  :py:class:`Molecule.from_file <openff.toolkit.topology.Molecule>`.
- `PR #631 <https://github.com/openforcefield/openff-toolkit/pull/631>`_: Fixes a bug in which calling
  :py:class:`unit_to_string <openff.toolkit.utils.utils.unit_to_string>` returned
  ``None`` when the unit is dimensionless. Now ``"dimensionless"`` is returned.
- `PR #630 <https://github.com/openforcefield/openff-toolkit/pull/630>`_: Closes issue `Issue #629
  <https://github.com/openforcefield/openff-toolkit/issues/629>`_ in which the wrong exception is raised when
  attempting to instantiate a :py:class:`ForceField <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>`
  from an unparsable string.

New features
""""""""""""
- `PR #632 <https://github.com/openforcefield/openff-toolkit/pull/632>`_: Adds
  :py:class:`ForceField.registered_parameter_handlers <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>`
- `PR #614 <https://github.com/openforcefield/openff-toolkit/pull/614>`_: Adds
  :py:class:`ToolkitRegistry.deregister_toolkit <openff.toolkit.utils.toolkits.ToolkitRegistry>`
  to de-register registered toolkits, which can include toolkit wrappers loaded into ``GLOBAL_TOOLKIT_REGISTRY``
  by default.
- `PR #656 <https://github.com/openforcefield/openff-toolkit/pull/656>`_: Adds
  a new allowed ``am1elf10`` option to the OpenEye implementation of
  :py:class:`assign_partial_charges <openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper>` which
  calculates the average partial charges at the AM1 level of theory using conformers selected using the ELF10 method.
- `PR #643 <https://github.com/openforcefield/openff-toolkit/pull/643>`_: Adds
  :py:class:`openforcefield.typing.engines.smirnoff.forcefield.get_available_force_fields <openff.toolkit.typing.engines.smirnoff.forcefield.get_available_force_fields>`,
  which returns paths to the files of force fields available through entry point plugins.


0.7.0 - Charge Increment Model, Proper Torsion interpolation, and new Molecule methods
--------------------------------------------------------------------------------------

This is a relatively large release, motivated by the idea that changing existing functionality is bad
so we shouldn't do it too often, but when we do change things we should do it all at once.

Here's a brief rundown of what changed, migration tips, and how to find more details in the full release notes below:

* To provide more consistent partial charges for a given molecule, existing conformers are now disregarded by default
  by ``Molecule.assign_partial_charges``. Instead, new conformers are generated for use in semiempirical calculations.
  Search for ``use_conformers``.
* Formal charges are now always returned as ``simtk.unit.Quantity`` objects, with units of elementary charge.
  To convert them to integers, use ``from simtk import unit`` and
  ``atom.formal_charge.value_in_unit(unit.elementary_charge)`` or
  ``mol.total_charge.value_in_unit(unit.elementary_charge)``.
  Search ``atom.formal_charge``.
* The OpenFF Toolkit now automatically reads and writes partial charges in SDF files. Search for
  ``atom.dprop.PartialCharges``.
* The OpenFF Toolkit now has different behavior for handling multi-molecule and multi-conformer SDF files. Search
  ``multi-conformer``.
* The OpenFF Toolkit now distinguishes between partial charges that are all-zero and partial charges that are unknown.
  Search ``partial_charges = None``.
* ``Topology.to_openmm`` now assigns unique atoms names by default. Search ``ensure_unique_atom_names``.
* Molecule equality checks are now done by graph comparison instead of SMILES comparison.
  Search ``Molecule.are_isomorphic``.
* The ``ChemicalEnvironment`` module was almost entirely removed, as it is an outdated duplicate of some Chemper
  functionality. Search ``ChemicalEnvironment``.
* ``TopologyMolecule.topology_particle_start_index`` has been removed from the ``TopologyMolecule`` API, since atoms
  and virtualsites are no longer contiguous in the ``Topology`` particle indexing system. Search
  ``topology_particle_start_index``.
* ``compute_wiberg_bond_orders`` has been renamed to ``assign_fractional_bond_orders``.

There are also a number of new features, such as:

* Support for ``ChargeIncrementModel`` sections in force fields.
* Support for ``ProperTorsion`` ``k`` interpolation in force fields using fractional bond orders.
* Support for AM1-Mulliken, Gasteiger, and other charge methods using the new ``assign_partial_charges`` methods.
* Support for AM1-Wiberg bond order calculation using either the OpenEye or RDKit/AmberTools backends and the
  ``assign_fractional_bond_orders`` methods.
* Initial (limited) interoperability with QCArchive, via ``Molecule.to_qcschema`` and ``from_qcschema``.
* A ``Molecule.visualize`` method.
* Several additional ``Molecule`` methods, including state enumeration and mapped SMILES creation.

**Major Feature: Support for the SMIRNOFF ChargeIncrementModel tag**

`The ChargeIncrementModel tag in the SMIRNOFF specification <https://openforcefield.github.io/standards/standards/smirnoff/#chargeincrementmodel-small-molecule-and-fragment-charges>`_
provides analagous functionality to AM1-BCC, except that instead of AM1-Mulliken charges, a number of different charge
methods can be called, and instead of a fixed library of two-atom charge corrections, an arbitrary number of
SMIRKS-based, N-atom charge corrections can be defined in the SMIRNOFF format.

The initial implementation of the SMIRNOFF ``ChargeIncrementModel`` tag accepts keywords for ``version``,
``partial_charge_method``, and ``number_of_conformers``. ``partial_charge_method`` can be any string, and it is
up to the ``ToolkitWrapper``'s ``compute_partial_charges`` methods to understand what they mean. For
geometry-independent ``partial_charge_method`` choices, ``number_of_conformers`` should be set to zero.

SMIRKS-based parameter application for ``ChargeIncrement`` parameters is different than other SMIRNOFF sections.
The initial implementation of ``ChargeIncrementModelHandler`` follows these rules:

* an atom can be subject to many ``ChargeIncrement`` parameters, which combine additively.
* a ``ChargeIncrement`` that matches a set of atoms is overwritten only if another ``ChargeIncrement``
  matches the same group of atoms, regardless of order. This overriding follows the normal SMIRNOFF hierarchy.

To give a concise example, what if a molecule ``A-B(-C)-D`` were being parametrized, and the force field
defined ``ChargeIncrement`` SMIRKS in the following order?

1) ``[A:1]-[B:2]``
2) ``[B:1]-[A:2]``
3) ``[A:1]-[B:2]-[C:3]``
4) ``[*:1]-[B:2](-[*:3])-[*:4]``
5) ``[D:1]-[B:2](-[*:3])-[*:4]``

In the case above, the ChargeIncrement from parameters 1 and 4 would NOT be applied to the molecule,
since another parameter matching the same set of atoms is specified further down in the parameter hierarchy
(despite those subsequent matches being in a different order).

Ultimately, the ChargeIncrement contributions from parameters 2, 3, and 5 would be summed and applied.

It's also important to identify a behavior that these rules were written to *avoid*: if not for the
"regardless of order" clause in the second rule, parameters 4 and 5 could actually have been applied six and two times,
respectively (due to symmetry in the SMIRKS and the use of wildcards). This situation could also arise as a result
of molecular symmetry. For example, a methyl group could match the SMIRKS ``[C:1]([H:2])([H:3])([H:4])`` six ways
(with different orderings of the three hydrogen atoms), but the user would almost certainly not intend for the charge
increments to be applied six times. The "regardless of order" clause was added specifically to address this.

In short, the first time a group of atoms becomes involved in a ``ChargeIncrement`` together, the OpenMM ``System`` gains a new
parameter "slot". Only another ``ChargeIncrement`` which applies to the exact same group of atoms (in any order) can
take over the "slot", pushing the original ``ChargeIncrement`` out.

**Major Feature: Support for ProperTorsion k value interpolation**

`Chaya Stern's work <https://chayast.github.io/fragmenter-manuscript/>`_
showed that we may be able to produce higher-quality proper torsion parameters by taking into
account the "partial bond order" of the torsion's central bond. We now have the machinery to compute AM1-Wiberg
partial bond orders for entire molecules using the ``assign_fractional_bond_orders`` methods of either  ``OpenEyeToolkitWrapper`` or ``AmberToolsToolkitWrapper``. The thought is that, if some simple electron population analysis shows
that a certain aromatic bond's order is 1.53, maybe rotations about that bond can be described well by interpolating
53% of the way between the single and double bond k values.

Full details of how to define a torsion-interpolating SMIRNOFF force fields are available in
`the ProperTorsions section of the SMIRNOFF specification <https://openforcefield.github.io/standards/standards/smirnoff/#fractional-torsion-bond-orders>`_.

Behavior changed
""""""""""""""""
- `PR #508 <https://github.com/openforcefield/openff-toolkit/pull/508>`_:
  In order to provide the same results for the same chemical species, regardless of input
  conformation,
  :py:class:`Molecule <openff.toolkit.topology.Molecule>`
  ``assign_partial_charges``, ``compute_partial_charges_am1bcc``, and
  ``assign_fractional_bond_orders`` methods now default to ignore input conformers
  and generate new conformer(s) of the molecule before running semiempirical calculations.
  Users can override this behavior by specifying the keyword argument
  ``use_conformers=molecule.conformers``.
- `PR #281 <https://github.com/openforcefield/openff-toolkit/pull/281>`_: Closes
  `Issue #250 <https://github.com/openforcefield/openff-toolkit/issues/250>`_
  by adding support for partial charge I/O in SDF. The partial charges are stored as a property in the
  SDF molecule block under the tag ``<atom.dprop.PartialCharge>``.
- `PR #281 <https://github.com/openforcefield/openff-toolkit/pull/281>`_: If a
  :py:class:`Molecule <openff.toolkit.topology.Molecule>`'s
  ``partial_charges`` attribute is set to ``None`` (the default value), calling ``to_openeye`` will
  now produce a OE molecule with partial charges set to ``nan``. This would previously produce an OE
  molecule with partial charges of 0.0, which was a loss of information, since it wouldn't be clear
  whether the original OFFMol's partial charges were REALLY all-zero as opposed to ``None``. OpenEye toolkit
  wrapper methods such as ``from_smiles`` and ``from_file`` now produce OFFMols with
  ``partial_charges = None`` when appropriate (previously these would produce OFFMols with
  all-zero charges, for the same reasoning as above).
- `PR #281 <https://github.com/openforcefield/openff-toolkit/pull/281>`_:
  :py:class:`Molecule <openff.toolkit.topology.Molecule>`
  ``to_rdkit``
  now sets partial charges on the RDAtom's ``PartialCharges`` property (this was previously set
  on the ``partial_charges`` property). If the
  :py:class:`Molecule <openff.toolkit.topology.Molecule>`'s partial_charges attribute is ``None``, this property
  will not be defined on the RDAtoms.
- `PR #281 <https://github.com/openforcefield/openff-toolkit/pull/281>`_:
  Enforce the behavior during SDF I/O that a SDF may contain multiple
  `molecules`, but that the OFF Toolkit
  does not assume that it contains multiple `conformers of the same molecule`. This is an
  important distinction, since otherwise there is ambiguity around whether properties of one
  entry in a SDF are shared among several molecule blocks or not, or how to resolve conflicts if properties
  are defined differently for several "conformers" of chemically-identical species (More info
  `here <https://docs.eyesopen.com/toolkits/python/oechemtk/oemol.html#dude-where-s-my-sd-data>`_).
  If the user requests the OFF Toolkit to write a multi-conformer
  :py:class:`Molecule <openff.toolkit.topology.Molecule>` to SDF, only the first conformer will be written.
  For more fine-grained control of writing properties, conformers, and partial charges, consider
  using ``Molecule.to_rdkit`` or ``Molecule.to_openeye`` and using the functionality offered by
  those packages.
- `PR #281 <https://github.com/openforcefield/openff-toolkit/pull/281>`_: Due to different
  constraints placed on the data types allowed by external toolkits, we make our best effort to
  preserve :py:class:`Molecule <openff.toolkit.topology.Molecule>`
  ``properties`` when converting molecules to other packages, but users should be aware that
  no guarantee of data integrity is made. The only data format for keys and values in the property dict that
  we will try to support through a roundtrip to another toolkit's Molecule object is ``string``.
- `PR #574 <https://github.com/openforcefield/openff-toolkit/pull/574>`_: Removed check that all
  partial charges are zero after assignment by ``quacpac`` when AM1BCC used for charge assignment.
  This check fails erroneously for cases in which the partial charge assignments are correctly all zero,
  such as for ``N#N``. It is also an unnecessary check given that ``quacpac`` will reliably indicate when
  it has failed to assign charges.
- `PR #597 <https://github.com/openforcefield/openff-toolkit/pull/597>`_: Energy-minimized sample systems
  with Parsley 1.1.0.
- `PR #558 <https://github.com/openforcefield/openff-toolkit/pull/558>`_: The
  :py:class:`Topology <openff.toolkit.topology.Topology>`
  particle indexing system now orders :py:class:`TopologyVirtualSites <openff.toolkit.topology.TopologyVirtualSite>`
  after all atoms.
- `PR #469 <https://github.com/openforcefield/openff-toolkit/pull/469>`_:
  When running :py:meth:`Topology.to_openmm <openff.toolkit.topology.Topology.to_openmm>`, unique atom names
  are generated if the provided atom names are not unique (overriding any existing atom names). This
  uniqueness extends only to atoms in the same molecule. To disable this behavior, set the kwarg
  ``ensure_unique_atom_names=False``.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_:
  :py:meth:`Molecule.__eq__ <openff.toolkit.topology.Molecule>` now uses the new
  :py:meth:`Molecule.are_isomorphic <openff.toolkit.topology.Molecule.are_isomorphic>` to perform the
  similarity checking.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_:
  The :py:meth:`Topology.from_openmm <openff.toolkit.topology.Topology.from_openmm>` and
  :py:meth:`Topology.add_molecule <openff.toolkit.topology.Topology.add_molecule>` methods now use the
  :py:meth:`Molecule.are_isomorphic <openff.toolkit.topology.Molecule.are_isomorphic>` method to match
  molecules.
- `PR #551 <https://github.com/openforcefield/openff-toolkit/pull/551>`_: Implemented the
  :py:meth:`ParameterHandler.get_parameter <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler.get_parameter>`
  function (would previously return ``None``).

API-breaking changes
""""""""""""""""""""
- `PR #471 <https://github.com/openforcefield/openff-toolkit/pull/471>`_: Closes
  `Issue #465 <https://github.com/openforcefield/openff-toolkit/issues/465>`_.
  ``atom.formal_charge`` and ``molecule.total_charge`` now return ``simtk.unit.Quantity`` objects
  instead of integers. To preserve backward compatibility, the setter for ``atom.formal_charge``
  can accept either a ``simtk.unit.Quantity`` or an integer.
- `PR #601 <https://github.com/openforcefield/openff-toolkit/pull/601>`_: Removes
  almost all of the previous
  :py:class:`ChemicalEnvironment <openff.toolkit.typing.chemistry.ChemicalEnvironment>`
  API, since this entire module was simply copied from
  `Chemper <https://github.com/MobleyLab/chemper>`_ several years ago and has fallen behind on updates.
  Currently only
  :py:meth:`ChemicalEnvironment.get_type <openff.toolkit.typing.chemistry.ChemicalEnvironment.get_type>`,
  :py:meth:`ChemicalEnvironment.validate <openff.toolkit.typing.chemistry.ChemicalEnvironment.validate>`,
  and an equivalent classmethod
  :py:meth:`ChemicalEnvironment.validate_smirks <openff.toolkit.typing.chemistry.ChemicalEnvironment.validate_smirks>`
  remain. Also, please comment on
  `this GitHub issue <https://github.com/MobleyLab/chemper/issues/90>`_ if you HAVE been using
  the previous extra functionality in this module and would like us to prioritize creation of a Chemper
  conda package.
- `PR #558 <https://github.com/openforcefield/openff-toolkit/pull/558>`_: Removes
  ``TopologyMolecule.topology_particle_start_index``, since the :py:class:`Topology <openff.toolkit.topology.Topology>`
  particle indexing system now orders :py:class:`TopologyVirtualSites <openff.toolkit.topology.TopologyVirtualSite>`
  after all atoms.
  :py:meth:`TopologyMolecule.atom_start_topology_index <openff.toolkit.topology.TopologyMolecule.atom_start_topology_index>`
  and
  :py:meth:`TopologyMolecule.virtual_particle_start_topology_index <openff.toolkit.topology.TopologyMolecule.virtual_particle_start_topology_index>`
  are still available to access the appropriate values in the respective topology indexing systems.
- `PR #508 <https://github.com/openforcefield/openff-toolkit/pull/508>`_:
  ``OpenEyeToolkitWrapper.compute_wiberg_bond_orders`` is now
  :py:meth:`OpenEyeToolkitWrapper.assign_fractional_bond_orders <openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.assign_fractional_bond_orders>`.
  The ``charge_model`` keyword is now ``bond_order_model``. The allowed values of this keyword have
  changed from ``am1`` and ``pm3`` to ``am1-wiberg`` and ``pm3-wiberg``, respectively.
- `PR #508 <https://github.com/openforcefield/openff-toolkit/pull/508>`_:
  ``Molecule.compute_wiberg_bond_orders`` is now
  :py:meth:`Molecule.assign_fractional_bond_orders <openff.toolkit.topology.Molecule.assign_fractional_bond_orders>`.
- `PR #595 <https://github.com/openforcefield/openff-toolkit/pull/595>`_: Removed functions
  ``openforcefield.utils.utils.temporary_directory`` and
  ``openforcefield.utils.utils.temporary_cd`` and replaced their behavior with
  ``tempfile.TemporaryDirectory()``.

New features
""""""""""""
- `PR #471 <https://github.com/openforcefield/openff-toolkit/pull/471>`_: Closes
  `Issue #208 <https://github.com/openforcefield/openff-toolkit/issues/208>`_
  by implementing support for the
  ``ChargeIncrementModel`` tag in the `SMIRNOFF specification <https://openforcefield.github.io/standards/standards/smirnoff/#chargeincrementmodel-small-molecule-and-fragment-charges>`_.
- `PR #471 <https://github.com/openforcefield/openff-toolkit/pull/471>`_: Implements
  ``Molecule.assign_partial_charges``, which calls one of the newly-implemented
  ``OpenEyeToolkitWrapper.assign_partial_charges``, and
  ``AmberToolsToolkitWrapper.assign_partial_charges``. ``strict_n_conformers`` is a
  optional boolean keyword argument indicating whether an ``IncorrectNumConformersError`` should be raised if an invalid
  number of conformers is supplied during partial charge calculation. For example, if two conformers are
  supplied, but ``partial_charge_method="AM1BCC"`` is also set, then there is no clear use for
  the second conformer. The previous behavior in this case was to raise a warning, and to preserve that
  behavior, ``strict_n_conformers`` defaults to a value of ``False``.
- `PR #471 <https://github.com/openforcefield/openff-toolkit/pull/471>`_: Adds
  keyword argument ``raise_exception_types`` (default: ``[Exception]``) to
  :py:meth:`ToolkitRegistry.call <openff.toolkit.utils.toolkits.ToolkitRegistry.call>`.
  The default value will provide the previous OpenFF Toolkit behavior, which is that the first ToolkitWrapper
  that can provide the requested method is called, and it either returns on success or raises an exception. This new
  keyword argument allows the ToolkitRegistry to *ignore* certain exceptions, but treat others as fatal.
  If ``raise_exception_types = []``, the ToolkitRegistry will attempt to call each ToolkitWrapper that provides the
  requested method and if none succeeds, a single ``ValueError`` will be raised, with text listing the
  errors that were raised by each ToolkitWrapper.
- `PR #601 <https://github.com/openforcefield/openff-toolkit/pull/601>`_: Adds
  :py:meth:`RDKitToolkitWrapper.get_tagged_smarts_connectivity <openff.toolkit.utils.toolkits.RDKitToolkitWrapper.get_tagged_smarts_connectivity>`
  and
  :py:meth:`OpenEyeToolkitWrapper.get_tagged_smarts_connectivity <openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.get_tagged_smarts_connectivity>`,
  which allow the use of either toolkit for smirks/tagged smarts validation.
- `PR #600 <https://github.com/openforcefield/openff-toolkit/pull/600>`_:
  Adds :py:meth:`ForceField.__getitem__ <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>`
  to look up ``ParameterHandler`` objects based on their string names.
- `PR #508 <https://github.com/openforcefield/openff-toolkit/pull/508>`_:
  Adds :py:meth:`AmberToolsToolkitWrapper.assign_fractional_bond_orders <openff.toolkit.utils.toolkits.AmberToolsToolkitWrapper.assign_fractional_bond_orders>`.
- `PR #469 <https://github.com/openforcefield/openff-toolkit/pull/469>`_: The
  :py:class:`Molecule <openff.toolkit.topology.Molecule>` class adds
  :py:meth:`Molecule.has_unique_atom_names <openff.toolkit.topology.Molecule.has_unique_atom_names>`
  and :py:meth:`Molecule.has_unique_atom_names <openff.toolkit.topology.Molecule.generate_unique_atom_names>`.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_:
  Adds to the :py:class:`Molecule <openff.toolkit.topology.Molecule>` class
  :py:meth:`Molecule.are_isomorphic <openff.toolkit.topology.Molecule.are_isomorphic>`
  and :py:meth:`Molecule.is_isomorphic_with <openff.toolkit.topology.Molecule.is_isomorphic_with>`
  and :py:meth:`Molecule.hill_formula <openff.toolkit.topology.Molecule.hill_formula>`
  and :py:meth:`Molecule.to_hill_formula <openff.toolkit.topology.Molecule.to_hill_formula>`
  and :py:meth:`Molecule.to_qcschema <openff.toolkit.topology.Molecule.to_qcschema>`
  and :py:meth:`Molecule.from_qcschema <openff.toolkit.topology.Molecule.from_qcschema>`
  and :py:meth:`Molecule.from_mapped_smiles <openff.toolkit.topology.Molecule.from_mapped_smiles>`
  and :py:meth:`Molecule.from_pdb_and_smiles <openff.toolkit.topology.Molecule.from_pdb_and_smiles>`
  and :py:meth:`Molecule.canonical_order_atoms <openff.toolkit.topology.Molecule.canonical_order_atoms>`
  and :py:meth:`Molecule.remap <openff.toolkit.topology.Molecule.remap>`

    .. note::
       The to_qcschema method accepts an extras dictionary which is passed into the validated qcelemental.models.Molecule
       object.

- `PR #506 <https://github.com/openforcefield/openff-toolkit/pull/506>`_:
  The :py:class:`Molecule <openff.toolkit.topology.Molecule>` class adds
  :py:meth:`Molecule.find_rotatable_bonds <openff.toolkit.topology.Molecule.find_rotatable_bonds>`
- `PR #521 <https://github.com/openforcefield/openff-toolkit/pull/521>`_:
  Adds :py:meth:`Molecule.to_inchi <openff.toolkit.topology.Molecule.to_inchi>`
  and :py:meth:`Molecule.to_inchikey <openff.toolkit.topology.Molecule.to_inchikey>`
  and :py:meth:`Molecule.from_inchi <openff.toolkit.topology.Molecule.from_inchi>`

    .. warning::
       InChI was not designed as an molecule interchange format and using it as one is not recommended. Many round trip
       tests will fail when using this format due to a loss of information. We have also added support for fixed
       hydrogen layer nonstandard InChI which can help in the case of tautomers, but overall creating molecules from InChI should be
       avoided.

- `PR #529 <https://github.com/openforcefield/openff-toolkit/pull/529>`_: Adds the ability to write out to XYZ files via
  :py:meth:`Molecule.to_file <openff.toolkit.topology.Molecule.to_file>` Both single frame and multiframe XYZ files are supported.
  Note reading from XYZ files will not be supported due to the lack of connectivity information.
- `PR #535 <https://github.com/openforcefield/openff-toolkit/pull/535>`_: Extends the the API for the
  :py:meth:`Molecule.to_smiles <openff.toolkit.topology.Molecule.to_smiles>` to allow for the creation of cmiles
  identifiers through combinations of isomeric, explicit hydrogen and mapped smiles, the default settings will return
  isomeric explicit hydrogen smiles as expected.

        .. warning::
           Atom maps can be supplied to the properties dictionary to modify which atoms have their map index included,
           if no map is supplied all atoms will be mapped in the order they appear in the
           :py:class:`Molecule <openff.toolkit.topology.Molecule>`.

- `PR #563 <https://github.com/openforcefield/openff-toolkit/pull/563>`_:
  Adds ``test_forcefields/ion_charges.offxml``, giving ``LibraryCharges`` for monatomic ions.
- `PR #543 <https://github.com/openforcefield/openff-toolkit/pull/543>`_:
  Adds 3 new methods to the :py:class:`Molecule <openff.toolkit.topology.Molecule>` class which allow the enumeration of molecule
  states. These are :py:meth:`Molecule.enumerate_tautomers <openff.toolkit.topology.Molecule.enumerate_tautomers>`,
  :py:meth:`Molecule.enumerate_stereoisomers <openff.toolkit.topology.Molecule.enumerate_stereoisomers>`,
  :py:meth:`Molecule.enumerate_protomers <openff.toolkit.topology.Molecule.enumerate_protomers>`

      .. warning::
         Enumerate protomers is currently only available through the OpenEye toolkit.

- `PR #573 <https://github.com/openforcefield/openff-toolkit/pull/573>`_:
  Adds ``quacpac`` error output to ``quacpac`` failure in ``Molecule.compute_partial_charges_am1bcc``.
- `PR #560 <https://github.com/openforcefield/openff-toolkit/issues/560>`_: Added visualization method to the the Molecule class.
- `PR #620 <https://github.com/openforcefield/openff-toolkit/pull/620>`_: Added the ability to register parameter handlers via entry point plugins. This functionality is accessible by initializing a ``ForceField`` with the ``load_plugins=True`` keyword argument.
- `PR #582 <https://github.com/openforcefield/openff-toolkit/pull/582>`_: Added fractional bond order interpolation
  Adds `return_topology` kwarg to
  :py:meth:`Forcefield.create_openmm_system <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system>`,
  which returns the processed topology along with the OpenMM ``System`` when ``True`` (default ``False``).

Tests added
"""""""""""
- `PR #558 <https://github.com/openforcefield/openff-toolkit/pull/558>`_: Adds tests ensuring
  that the new Topology particle indexing system are properly implemented, and that TopologyVirtualSites
  reference the correct TopologyAtoms.
- `PR #469 <https://github.com/openforcefield/openff-toolkit/pull/469>`_: Added round-trip SMILES test
  to add coverage for :py:meth:`Molecule.from_smiles <openff.toolkit.topology.Molecule.from_smiles>`.
- `PR #469 <https://github.com/openforcefield/openff-toolkit/pull/469>`_: Added tests for unique atom
  naming behavior in  :py:meth:`Topology.to_openmm <openff.toolkit.topology.Topology.to_openmm>`, as
  well as tests of the ``ensure_unique_atom_names=False`` kwarg disabling this behavior.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_: Added tests for
  :py:meth:`Molecule.hill_formula <openff.toolkit.topology.Molecule.hill_formula>` and
  :py:meth:`Molecule.to_hill_formula <openff.toolkit.topology.Molecule.to_hill_formula>` for the
  various supported input types.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_: Added round-trip test for
  :py:meth:`Molecule.from_qcschema <openff.toolkit.topology.Molecule.from_qcschema>` and
  :py:meth:`Molecule.to_qcschema <openff.toolkit.topology.Molecule.to_qcschema>`.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_: Added tests for
  :py:meth:`Molecule.is_isomorphic_with <openff.toolkit.topology.Molecule.is_isomorphic_with>` and
  :py:meth:`Molecule.are_isomorphic <openff.toolkit.topology.Molecule.are_isomorphic>`
  with various levels of isomorphic graph matching.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_: Added toolkit dependent tests
  for :py:meth:`Molecule.canonical_order_atoms <openff.toolkit.topology.Molecule.canonical_order_atoms>`
  due to differences in the algorithms used.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_: Added a test for
  :py:meth:`Molecule.from_mapped_smiles <openff.toolkit.topology.Molecule.from_mapped_smiles>` using
  the molecule from issue #412 to ensure it is now fixed.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_: Added a test for
  :py:meth:`Molecule.remap <openff.toolkit.topology.Molecule.remap>`, this also checks for expected
  error when the mapping is not complete.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_: Added tests for
  :py:meth:`Molecule.from_pdb_and_smiles <openff.toolkit.topology.Molecule.from_pdb_and_smiles>`
  to check for a correct combination of smiles and PDB and incorrect combinations.
- `PR #509 <https://github.com/openforcefield/openff-toolkit/pull/509>`_: Added test for
  :py:meth:`Molecule.chemical_environment_matches <openff.toolkit.topology.Molecule.chemical_environment_matches>`
  to check that the complete set of matches is returned.
- `PR #509 <https://github.com/openforcefield/openff-toolkit/pull/509>`_: Added test for
  :py:meth:`Forcefield.create_openmm_system <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system>`
  to check that a protein system can be created.
- `PR #506 <https://github.com/openforcefield/openff-toolkit/pull/506>`_: Added a test for the molecule
  identified in issue #513 as losing aromaticity when converted to rdkit.
- `PR #506 <https://github.com/openforcefield/openff-toolkit/pull/506>`_: Added a verity of toolkit dependent tests
  for identifying rotatable bonds while ignoring the user requested types.
- `PR #521 <https://github.com/openforcefield/openff-toolkit/pull/521>`_: Added toolkit independent round-trip InChI
  tests which add coverage for :py:meth:`Molecule.to_inchi <openff.toolkit.topology.Molecule.to_inchi>` and
  :py:meth:`Molecule.from_inchi <openff.toolkit.topology.Molecule.from_inchi>`. Also added coverage for bad inputs and
  :py:meth:`Molecule.to_inchikey <openff.toolkit.topology.Molecule.to_inchikey>`.
- `PR #529 <https://github.com/openforcefield/openff-toolkit/pull/529>`_: Added to XYZ file coverage tests.
- `PR #563 <https://github.com/openforcefield/openff-toolkit/pull/563>`_: Added `LibraryCharges` parameterization test
  for monatomic ions in ``test_forcefields/ion_charges.offxml``.
- `PR #543 <https://github.com/openforcefield/openff-toolkit/pull/543>`_: Added tests to assure that state enumeration can
  correctly find molecules tautomers, stereoisomers and protomers when possible.
- `PR #573 <https://github.com/openforcefield/openff-toolkit/pull/573>`_: Added test for ``quacpac`` error output
  for ``quacpac`` failure in ``Molecule.compute_partial_charges_am1bcc``.
- `PR #579 <https://github.com/openforcefield/openff-toolkit/pull/579>`_: Adds regression tests to ensure RDKit can be
  be used to write multi-model PDB files.
- `PR #582 <https://github.com/openforcefield/openff-toolkit/pull/582>`_: Added fractional bond order interpolation tests,
  tests for :py:class:`ValidatedDict <openff.toolkit.utils.collections.ValidatedDict>`.


Bugfixes
""""""""
- `PR #558 <https://github.com/openforcefield/openff-toolkit/pull/558>`_: Fixes a bug where
  :py:meth:`TopologyVirtualSite.atoms <openff.toolkit.topology.TopologyVirtualSite.atoms>` would
  not correctly apply ``TopologyMolecule`` atom ordering on top of the reference molecule ordering,
  in cases where the same molecule appears multiple times, but in a different order, in the same Topology.
- `Issue #460 <https://github.com/openforcefield/openff-toolkit/issues/460>`_: Creates unique atom
  names in :py:meth:`Topology.to_openmm <openff.toolkit.topology.Topology.to_openmm>` if the existing
  ones are not unique. The lack of unique atom names had been causing problems in workflows involving
  downstream tools that expect unique atom names.
- `Issue #448 <https://github.com/openforcefield/openff-toolkit/issues/448>`_: We can now make molecules
  from mapped smiles using :py:meth:`Molecule.from_mapped_smiles <openff.toolkit.topology.Molecule.from_mapped_smiles>`
  where the order will correspond to the indeing used in the smiles.
  Molecules can also be re-indexed at any time using the
  :py:meth:`Molecule.remap <openff.toolkit.topology.Molecule.remap>`.
- `Issue #462 <https://github.com/openforcefield/openff-toolkit/issues/462>`_: We can now instance the
  :py:class:`Molecule <openff.toolkit.topology.Molecule>` from a QCArchive entry record instance or dictionary
  representation.
- `Issue #412 <https://github.com/openforcefield/openff-toolkit/issues/412>`_: We can now instance the
  :py:class:`Molecule <openff.toolkit.topology.Molecule>` using
  :py:meth:`Molecule.from_mapped_smiles <openff.toolkit.topology.Molecule.from_mapped_smiles>`. This resolves
  an issue caused by RDKit considering atom map indices to be a distinguishing feature of an atom, which led
  to erroneous definition of chirality (as otherwise symmetric substituents would be seen as different).
  We anticipate that this will reduce the number of times you need to
  type ``allow_undefined_stereo=True`` when processing molecules that do not actually contain stereochemistrty.
- `Issue #513 <https://github.com/openforcefield/openff-toolkit/issues/513>`_: The
  :py:meth:`Molecule.to_rdkit <openff.toolkit.topology.Molecule.to_rdkit>` now re-sets the aromaticity model
  after sanitizing the molecule.
- `Issue #500 <https://github.com/openforcefield/openff-toolkit/issues/500>`_: The
  :py:meth:`Molecule.find_rotatable_bonds <openff.toolkit.topology.Molecule.find_rotatable_bonds>` has been added
  which returns a list of rotatable :py:class:`Bond <openff.toolkit.topology.Bond>` instances for the molecule.
- `Issue #491 <https://github.com/openforcefield/openff-toolkit/issues/491>`_: We can now parse large molecules without hitting a match limit cap.
- `Issue #474 <https://github.com/openforcefield/openff-toolkit/issues/474>`_: We can now  convert molecules to InChI and
  InChIKey and from InChI.
- `Issue #523 <https://github.com/openforcefield/openff-toolkit/issues/523>`_: The
  :py:meth:`Molecule.to_file <openff.toolkit.topology.Molecule.to_file>` method can now correctly write to ``MOL``
  files, in line with the supported file type list.
- `Issue #568 <https://github.com/openforcefield/openff-toolkit/issues/568>`_: The
  :py:meth:`Molecule.to_file <openff.toolkit.topology.Molecule.to_file>` can now correctly write multi-model PDB files
  when using the RDKit backend toolkit.


Examples added
""""""""""""""
- `PR #591 <https://github.com/openforcefield/openff-toolkit/pull/591>`_ and
  `PR #533 <https://github.com/openforcefield/openff-toolkit/pull/533>`_: Adds an
  `example notebook and utility to compute conformer energies <https://github.com/openforcefield/openff-toolkit/blob/master/examples/conformer_energies>`_.
  This example is made to be reverse-compatible with the 0.6.0 OpenFF Toolkit release.
- `PR #472 <https://github.com/openforcefield/openff-toolkit/pull/472>`_: Adds an example notebook
  `QCarchive_interface.ipynb <https://github.com/openforcefield/openff-toolkit/blob/master/examples/QCArchive_interface/QCarchive_interface.ipynb>`_
  which shows users how to instance the :py:class:`Molecule <openff.toolkit.topology.Molecule>` from
  a QCArchive entry level record and calculate the energy using RDKit through QCEngine.



0.6.0 - Library Charges
-----------------------

This release adds support for a new SMIRKS-based charge assignment method,
`Library Charges <https://openforcefield.github.io/standards/standards/smirnoff/#librarycharges-library-charges-for-polymeric-residues-and-special-solvent-models>`_.
The addition of more charge assignment methods opens the door for new types of
experimentation, but also introduces several complex behaviors and failure modes.
Accordingly, we have made changes
to the charge assignment infrastructure to check for cases when partial charges do
not sum to the formal charge of the molecule, or when no charge assignment method is able
to generate charges for a molecule. More detailed explanation of the new errors that may be raised and
keywords for overriding them are in the "Behavior Changed" section below.


With this release, we update ``test_forcefields/tip3p.offxml`` to be a working example of assigning LibraryCharges.
However, we do not provide any force field files to assign protein residue ``LibraryCharges``.
If you are interested in translating an existing protein FF to SMIRNOFF format or developing a new one, please
feel free to contact us on the `Issue tracker <https://github.com/openforcefield/openff-toolkit/issues>`_ or open a
`Pull Request <https://github.com/openforcefield/openff-toolkit/pulls>`_.


New features
""""""""""""
- `PR #433 <https://github.com/openforcefield/openff-toolkit/pull/433>`_: Closes
  `Issue #25 <https://github.com/openforcefield/openff-toolkit/issues/25>`_ by adding
  initial support for the
  `LibraryCharges tag in the SMIRNOFF specification <https://openforcefield.github.io/standards/standards/smirnoff/#librarycharges-library-charges-for-polymeric-residues-and-special-solvent-models>`_
  using
  :py:class:`LibraryChargeHandler <openff.toolkit.typing.engines.smirnoff.parameters.LibraryChargeHandler>`.
  For a molecule to have charges assigned using Library Charges, all of its atoms must be covered by
  at least one ``LibraryCharge``. If an atom is covered by multiple ``LibraryCharge`` s, then the last
  ``LibraryCharge`` matched will be applied (per the hierarchy rules in the SMIRNOFF format).

  This functionality is thus able to apply per-residue charges similar to those in traditional
  protein force fields. At this time, there is no concept of "residues" or "fragments" during
  parametrization, so it is not possible to assign charges to `some` atoms in a molecule using
  ``LibraryCharge`` s, but calculate charges for other atoms in the same molecule using a different
  method. To assign charges to a protein, LibraryCharges SMARTS must be provided for
  the residues and protonation states in the molecule, as well as for any capping groups
  and post-translational modifications that are present.

  It is valid for ``LibraryCharge`` SMARTS to `partially` overlap one another. For example, a molecule
  consisting of atoms ``A-B-C`` connected by single bonds could be matched by a SMIRNOFF
  ``LibraryCharges`` section containing two ``LibraryCharge`` SMARTS: ``A-B`` and ``B-C``. If
  listed in that order, the molecule would be assigned the ``A`` charge from the ``A-B`` ``LibraryCharge``
  element and the ``B`` and ``C`` charges from the ``B-C`` element. In testing, these types of
  partial overlaps were found to frequently be sources of undesired behavior, so it is recommended
  that users define whole-molecule ``LibraryCharge`` SMARTS whenever possible.

- `PR #455 <https://github.com/openforcefield/openff-toolkit/pull/455>`_: Addresses
  `Issue #393 <https://github.com/openforcefield/openff-toolkit/issues/393>`_ by adding
  :py:meth:`ParameterHandler.attribute_is_cosmetic <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler.attribute_is_cosmetic>`
  and
  :py:meth:`ParameterType.attribute_is_cosmetic <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType.attribute_is_cosmetic>`,
  which return True if the provided attribute name is defined for the queried object
  but does not correspond to an allowed value in the SMIRNOFF spec.

Behavior changed
""""""""""""""""
- `PR #433 <https://github.com/openforcefield/openff-toolkit/pull/433>`_: If a molecule
  can not be assigned charges by any charge-assignment method, an
  ``openforcefield.typing.engines.smirnoff.parameters.UnassignedMoleculeChargeException``
  will be raised. Previously, creating a system without either ``ToolkitAM1BCCHandler`` or
  the ``charge_from_molecules`` keyword argument to ``ForceField.create_openmm_system`` would
  produce an OpenMM ``System`` where the molecule has zero charge on all atoms. However, given that we
  will soon be adding more options for charge assignment, it is important that
  failures not be silent. Molecules with zero charge can still be produced by setting the
  ``Molecule.partial_charges`` array to be all zeroes, and including the molecule in the
  ``charge_from_molecules`` keyword argument to ``create_openmm_system``.
- `PR #433 <https://github.com/openforcefield/openff-toolkit/pull/433>`_: Due to risks
  introduced by permitting charge assignment using partially-overlapping ``LibraryCharge`` s,
  the toolkit will now raise a
  ``openforcefield.typing.engines.smirnoff.parameters.NonIntegralMoleculeChargeException``
  if the sum of partial charges on a molecule are found to be more than 0.01 elementary charge units
  different than the molecule's formal charge. This exception can be overridden by providing
  the ``allow_nonintegral_charges=True`` keyword argument to ``ForceField.create_openmm_system``.




Tests added
"""""""""""
- `PR #430 <https://github.com/openforcefield/openff-toolkit/pull/430>`_: Added test for
  Wiberg Bond Order implemented in OpenEye Toolkits. Molecules taken from
  DOI:10.5281/zenodo.3405489 . Added by Sukanya Sasmal.
- `PR #569 <https://github.com/openforcefield/openff-toolkit/pull/569>`_: Added round-trip tests for more serialization formats (dict, YAML, TOML, JSON, BSON, messagepack, pickle). Note that some are unsupported, but the tests raise the appropriate error.


Bugfixes
""""""""
- `PR #431 <https://github.com/openforcefield/openff-toolkit/pull/431>`_: Fixes an issue
  where ``ToolkitWrapper`` objects would improperly search for functionality in the
  ``GLOBAL_TOOLKIT_REGISTRY``, even though a specific ``ToolkitRegistry`` was requested for an
  operation.
- `PR #439 <https://github.com/openforcefield/openff-toolkit/pull/439>`_: Fixes
  `Issue #438 <https://github.com/openforcefield/openff-toolkit/issues/438>`_, by replacing
  call to NetworkX ``Graph.node`` with call to ``Graph.nodes``, per
  `2.4 migration guide <https://networkx.github.io/documentation/stable/release/release_2.4.html>`_.

Files modified
""""""""""""""
- `PR #433 <https://github.com/openforcefield/openff-toolkit/pull/433>`_: Updates
  the previously-nonfunctional ``test_forcefields/tip3p.offxml`` to a functional state
  by updating it to the SMIRNOFF
  0.3 specification, and specifying atomic charges using the ``LibraryCharges`` tag.


0.5.1 - Adding the parameter coverage example notebook
------------------------------------------------------

This release contains a new notebook example,
`check_parameter_coverage.ipynb <https://github.com/openforcefield/openff-toolkit/blob/master/examples/check_dataset_parameter_coverage/check_parameter_coverage.ipynb>`_,
which loads sets of molecules, checks whether they are parameterizable,
and generates reports of chemical motifs that are not.
It also fixes several simple issues, improves warnings and docstring text,
and removes unused files.

The parameter coverage example notebook goes hand-in-hand with the
release candidate of our initial force field,
`openff-1.0.0-RC1.offxml <https://github.com/openforcefield/openforcefields>`_
, which will be temporarily available until the official force
field release is made in October.
Our goal in publishing this notebook alongside our first major refitting is to allow interested
users to check whether there is parameter coverage for their molecules of interest.
If the force field is unable to parameterize a molecule, this notebook will generate
reports of the specific chemistry that is not covered. We understand that many organizations
in our field have restrictions about sharing specific molecules, and the outputs from this
notebook can easily be cropped to communicate unparameterizable chemistry without revealing
the full structure.

The force field release candidate is in our new refit force field package,
`openforcefields <https://github.com/openforcefield/openforcefields>`_.
This package is now a part of the Open Force Field Toolkit conda recipe, along with the original
`smirnoff99Frosst <https://github.com/openforcefield/smirnoff99Frosst>`_ line of force fields.

Once the ``openforcefields`` conda package is installed, you can load the release candidate using:

``ff = ForceField('openff-1.0.0-RC1.offxml')``

The release candidate will be removed when the official force field,
``openff-1.0.0.offxml``, is released in early October.

Complete details about this release are below.

Example added
"""""""""""""
- `PR #419 <https://github.com/openforcefield/openff-toolkit/pull/419>`_: Adds
  an example notebook
  `check_parameter_coverage.ipynb <https://github.com/openforcefield/openff-toolkit/blob/master/examples/check_dataset_parameter_coverage/check_parameter_coverage.ipynb>`_
  which shows how to use the toolkit to check a molecule
  dataset for missing parameter coverage, and provides functionality to output
  tagged SMILES and 2D drawings of the unparameterizable chemistry.


New features
""""""""""""
- `PR #419 <https://github.com/openforcefield/openff-toolkit/pull/419>`_: Unassigned
  valence parameter exceptions now include a list of tuples of
  :py:class:`TopologyAtom <openff.toolkit.topology.TopologyAtom>`
  which were unable to be parameterized (``exception.unassigned_topology_atom_tuples``)
  and the class of the
  :py:class:`ParameterHandler <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>`
  that raised the exception (``exception.handler_class``).
- `PR #425 <https://github.com/openforcefield/openff-toolkit/pull/425>`_: Implements
  Trevor Gokey's suggestion from
  `Issue #411 <https://github.com/openforcefield/openff-toolkit/issues/411>`_, which
  enables pickling of
  :py:class:`ForceFields <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>`
  and
  :py:class:`ParameterHandlers <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>`.
  Note that, while XML representations of ``ForceField``\ s are stable and conform to the SMIRNOFF
  specification, the pickled ``ForceField``\ s that this functionality enables are not guaranteed
  to be compatible with future toolkit versions.

Improved documentation and warnings
"""""""""""""""""""""""""""""""""""
- `PR #425 <https://github.com/openforcefield/openff-toolkit/pull/425>`_: Addresses
  `Issue #410 <https://github.com/openforcefield/openff-toolkit/issues/410>`_, by explicitly
  having toolkit warnings print ``Warning:`` at the beginning of each warning, and adding
  clearer language to the warning produced when the OpenEye Toolkits can not be loaded.
- `PR #425 <https://github.com/openforcefield/openff-toolkit/pull/425>`_: Addresses
  `Issue #421 <https://github.com/openforcefield/openff-toolkit/issues/421>`_ by
  adding type/shape information to all Molecule partial charge and conformer docstrings.
- `PR #425 <https://github.com/openforcefield/openff-toolkit/pull/425>`_: Addresses
  `Issue #407 <https://github.com/openforcefield/openff-toolkit/issues/421>`_ by
  providing a more extensive explanation of why we don't use RDKit's mol2 parser
  for molecule input.

Bugfixes
""""""""
- `PR #419 <https://github.com/openforcefield/openff-toolkit/pull/419>`_: Fixes
  `Issue #417 <https://github.com/openforcefield/openff-toolkit/issues/417>`_ and
  `Issue #418 <https://github.com/openforcefield/openff-toolkit/issues/418>`_, where
  :py:meth:`RDKitToolkitWrapper.from_file <openff.toolkit.utils.toolkits.RDKitToolkitWrapper.from_file>`
  would disregard the ``allow_undefined_stereo`` kwarg and skip the first molecule
  when reading a SMILES file.


Files removed
"""""""""""""
- `PR #425 <https://github.com/openforcefield/openff-toolkit/pull/425>`_: Addresses
  `Issue #424 <https://github.com/openforcefield/openff-toolkit/issues/424>`_ by
  deleting the unused files ``openforcefield/typing/engines/smirnoff/gbsaforces.py``
  and ``openforcefield/tests/test_smirnoff.py``. ``gbsaforces.py`` was only used internally
  and ``test_smirnoff.py`` tested unsupported functionality from before the 0.2.0 release.




0.5.0 - GBSA support and quality-of-life improvements
-----------------------------------------------------

This release adds support for the
`GBSA tag in the SMIRNOFF specification <https://openforcefield.github.io/standards/standards/smirnoff/#gbsa>`_.
Currently, the ``HCT``, ``OBC1``, and ``OBC2`` models (corresponding to AMBER keywords
``igb=1``, ``2``, and ``5``, respectively) are supported, with the ``OBC2`` implementation being
the most flexible. Unfortunately, systems produced
using these keywords are not yet transferable to other simulation packages via ParmEd, so users are restricted
to using OpenMM to simulate systems with GBSA.

OFFXML files containing GBSA parameter definitions are available,
and can be loaded in addition to existing parameter sets (for example, with the command
``ForceField('test_forcefields/smirnoff99Frosst.offxml', 'test_forcefields/GBSA_OBC1-1.0.offxml')``).
A manifest of new SMIRNOFF-format GBSA files is below.


Several other user-facing improvements have been added, including easier access to indexed attributes,
which are now accessible as ``torsion.k1``, ``torsion.k2``, etc. (the previous access method
``torsion.k`` still works as well). More details of the new features and several bugfixes are listed below.

New features
""""""""""""
- `PR #363 <https://github.com/openforcefield/openff-toolkit/pull/363>`_: Implements
  :py:class:`GBSAHandler <openff.toolkit.typing.engines.smirnoff.parameters.GBSAHandler>`,
  which supports the
  `GBSA tag in the SMIRNOFF specification <https://openforcefield.github.io/standards/standards/smirnoff/#gbsa>`_.
  Currently, only GBSAHandlers with ``gb_model="OBC2"`` support
  setting non-default values for the ``surface_area_penalty`` term (default ``5.4*calories/mole/angstroms**2``),
  though users can zero the SA term for ``OBC1`` and ``HCT`` models by setting ``sa_model="None"``.
  No model currently supports setting ``solvent_radius`` to any value other than ``1.4*angstroms``.
  Files containing experimental SMIRNOFF-format implementations of ``HCT``, ``OBC1``, and ``OBC2`` are
  included with this release (see below). Additional details of these models, including literature references,
  are available on the
  `SMIRNOFF specification page <https://openforcefield.github.io/standards/standards/smirnoff/#supported-generalized-born-gb-models>`_.

    .. warning :: The current release of ParmEd
      `can not transfer GBSA models produced by the Open Force Field Toolkit
      to other simulation packages
      <https://github.com/ParmEd/ParmEd/blob/3.2.0/parmed/openmm/topsystem.py#L148-L150>`_.
      These GBSA forces are currently only computable using OpenMM.

- `PR #363 <https://github.com/openforcefield/openff-toolkit/pull/363>`_: When using
  :py:meth:`Topology.to_openmm() <openff.toolkit.topology.Topology.to_openmm>`, periodic
  box vectors are now transferred from the Open Force Field Toolkit Topology
  into the newly-created OpenMM Topology.
- `PR #377 <https://github.com/openforcefield/openff-toolkit/pull/377>`_: Single indexed parameters in
  :py:class:`ParameterHandler <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>`
  and :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>`
  can now be get/set through normal attribute syntax in addition to the list syntax.
- `PR #394 <https://github.com/openforcefield/openff-toolkit/pull/394>`_: Include element and atom name
  in error output when there are missing valence parameters during molecule parameterization.

Bugfixes
""""""""
- `PR #385 <https://github.com/openforcefield/openff-toolkit/pull/385>`_: Fixes
  `Issue #346 <https://github.com/openforcefield/openff-toolkit/issues/346>`_ by
  having :py:meth:`OpenEyeToolkitWrapper.compute_partial_charges_am1bcc <openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.compute_partial_charges_am1bcc>`
  fall back to using standard AM1-BCC if AM1-BCC ELF10 charge generation raises
  an error about "trans COOH conformers"
- `PR #399 <https://github.com/openforcefield/openff-toolkit/pull/399>`_: Fixes
  issue where
  :py:class:`ForceField <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>`
  constructor would ignore ``parameter_handler_classes`` kwarg.
- `PR #400 <https://github.com/openforcefield/openff-toolkit/pull/400>`_: Makes
  link-checking tests retry three times before failing.



Files added
"""""""""""
- `PR #363 <https://github.com/openforcefield/openff-toolkit/pull/363>`_: Adds
  ``test_forcefields/GBSA_HCT-1.0.offxml``, ``test_forcefields/GBSA_OBC1-1.0.offxml``,
  and ``test_forcefields/GBSA_OBC2-1.0.offxml``, which are experimental implementations
  of GBSA models. These are primarily used in validation tests against OpenMM's models, and
  their version numbers will increment if bugfixes are necessary.

0.4.1 - Bugfix Release
----------------------

This update fixes several toolkit bugs that have been reported by the community.
Details of these bugfixes are provided below.

It also refactors how
:py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>`
and
:py:class:`ParameterHandler <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>`
store their attributes, by introducing
:py:class:`ParameterAttribute <openff.toolkit.typing.engines.smirnoff.parameters.ParameterAttribute>`
and
:py:class:`IndexedParameterAttribute <openff.toolkit.typing.engines.smirnoff.parameters.IndexedParameterAttribute>`.
These new attribute-handling classes provide a consistent backend which should simplify manipulation of parameters
and implementation of new handlers.

Bug fixes
"""""""""
- `PR #329 <https://github.com/openforcefield/openff-toolkit/pull/329>`_: Fixed a
  bug where the two
  :py:class:`BondType <openff.toolkit.typing.engines.smirnoff.parameters.BondHandler.BondType>`
  parameter attributes ``k`` and ``length`` were treated as indexed attributes. (``k`` and
  ``length`` values that correspond to specific bond orders will be indexed under
  ``k_bondorder1``, ``k_bondorder2``, etc when implemented in the future)
- `PR #329 <https://github.com/openforcefield/openff-toolkit/pull/329>`_: Fixed a
  bug that allowed setting indexed attributes to single values instead of strictly lists.
- `PR #370 <https://github.com/openforcefield/openff-toolkit/pull/370>`_: Fixed a
  bug in the API where
  :py:class:`BondHandler <openff.toolkit.typing.engines.smirnoff.parameters.BondHandler>`,
  :py:class:`ProperTorsionHandler <openff.toolkit.typing.engines.smirnoff.parameters.ProperTorsionHandler>`
  , and
  :py:class:`ImproperTorsionHandler <openff.toolkit.typing.engines.smirnoff.parameters.ImproperTorsionHandler>`
  exposed non-functional indexed parameters.
- `PR #351 <https://github.com/openforcefield/openff-toolkit/pull/351>`_: Fixes
  `Issue #344 <https://github.com/openforcefield/openff-toolkit/issues/344>`_,
  in which the main :py:class:`FrozenMolecule <openff.toolkit.topology.FrozenMolecule>`
  constructor and several other Molecule-construction functions ignored or did not
  expose the ``allow_undefined_stereo`` keyword argument.
- `PR #351 <https://github.com/openforcefield/openff-toolkit/pull/351>`_: Fixes
  a bug where a molecule which previously generated a SMILES using one cheminformatics toolkit
  returns the same SMILES, even though a different toolkit (which would generate
  a different SMILES for the molecule) is explicitly called.
- `PR #354 <https://github.com/openforcefield/openff-toolkit/pull/354>`_: Fixes
  the error message that is printed if an unexpected parameter attribute is found while loading
  data into a :py:class:`ForceField <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>`
  (now instructs users to specify ``allow_cosmetic_attributes`` instead of ``permit_cosmetic_attributes``)
- `PR #364 <https://github.com/openforcefield/openff-toolkit/pull/364>`_: Fixes
  `Issue #362 <https://github.com/openforcefield/openff-toolkit/issues/362>`_ by
  modifying
  :py:meth:`OpenEyeToolkitWrapper.from_smiles <openff.toolkit.utils.toolkits.OpenEyeToolkitWrapper.from_smiles>`
  and
  :py:meth:`RDKitToolkitWrapper.from_smiles <openff.toolkit.utils.toolkits.RDKitToolkitWrapper.from_smiles>`
  to make implicit hydrogens explicit before molecule creation. These functions also
  now raise an error if the optional keyword ``hydrogens_are_explicit=True`` but the
  SMILES are interpreted by the backend cheminformatic toolkit as having implicit
  hydrogens.
- `PR #371 <https://github.com/openforcefield/openff-toolkit/pull/371>`_: Fixes
  error when reading early SMIRNOFF 0.1 spec files enclosed by a top-level ``SMIRFF`` tag.

.. note ::
  The enclosing ``SMIRFF`` tag is present only in legacy files.
  Since developing a formal specification, the only acceptable top-level tag value in a SMIRNOFF data structure is
  ``SMIRNOFF``.

Code enhancements
"""""""""""""""""
- `PR #329 <https://github.com/openforcefield/openff-toolkit/pull/329>`_:
  :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>`
  was refactored to improve its extensibility. It is now possible to create new parameter
  types by using the new descriptors
  :py:class:`ParameterAttribute <openff.toolkit.typing.engines.smirnoff.parameters.ParameterAttribute>`
  and
  :py:class:`IndexedParameterAttribute <openff.toolkit.typing.engines.smirnoff.parameters.IndexedParameterAttribute>`.
- `PR #357 <https://github.com/openforcefield/openff-toolkit/pull/357>`_: Addresses
  `Issue #356 <https://github.com/openforcefield/openff-toolkit/issues/356>`_ by raising
  an informative error message if a user attempts to load an OpenMM topology which
  is probably missing connectivity information.



Force fields added
""""""""""""""""""
- `PR #368 <https://github.com/openforcefield/openff-toolkit/pull/368>`_: Temporarily adds
  ``test_forcefields/smirnoff99frosst_experimental.offxml`` to address hierarchy problems, redundancies, SMIRKS
  pattern typos etc., as documented in `issue #367 <https://github.com/openforcefield/openff-toolkit/issues/367>`_.
  Will ultimately be propagated to an updated force field in the ``openforcefield/smirnoff99frosst`` repo.
- `PR #371 <https://github.com/openforcefield/openff-toolkit/pull/371>`_: Adds
  ``test_forcefields/smirff99Frosst_reference_0_1_spec.offxml``, a SMIRNOFF 0.1 spec file enclosed by the legacy
  ``SMIRFF`` tag. This file is used in backwards-compatibility testing.



0.4.0 - Performance optimizations and support for SMIRNOFF 0.3 specification
----------------------------------------------------------------------------

This update contains performance enhancements that significantly reduce the time to create OpenMM systems for topologies containing many molecules via :py:meth:`ForceField.create_openmm_system <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system>`.

This update also introduces the `SMIRNOFF 0.3 specification <https://open-forcefield-toolkit.readthedocs.io/en/0.4.0/smirnoff.html>`_.
The spec update is the result of discussions about how to handle the evolution of data and parameter types as further functional forms are added to the SMIRNOFF spec.


We provide methods to convert SMIRNOFF 0.1 and 0.2 force fields written with the XML serialization (``.offxml``) to the SMIRNOFF 0.3 specification.
These methods are called automatically when loading a serialized SMIRNOFF data representation written in the 0.1 or 0.2 specification.
This functionality allows the toolkit to continue to read files containing SMIRNOFF 0.2 spec force fields, and also implements backwards-compatibility for SMIRNOFF 0.1 spec force fields.


.. warning :: The SMIRNOFF 0.1 spec did not contain fields for several energy-determining parameters that are exposed in later SMIRNOFF specs.
  Thus, when reading SMIRNOFF 0.1 spec data, the toolkit must make assumptions about the values that should be added for the newly-required fields.
  The values that are added include 1-2, 1-3 and 1-5 scaling factors, cutoffs, and long-range treatments for nonbonded interactions.
  Each assumption is printed as a warning during the conversion process.
  Please carefully review the warning messages to ensure that the conversion is providing your desired behavior.



`SMIRNOFF 0.3 specification updates <https://open-forcefield-toolkit.readthedocs.io/en/0.4.0/smirnoff.html>`_
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
* The SMIRNOFF 0.3 spec introduces versioning for each individual parameter section, allowing asynchronous updates to the features of each parameter class.
  The top-level ``SMIRNOFF`` tag, containing information like ``aromaticity_model``, ``Author``, and ``Date``, still has a version (currently 0.3).
  But, to allow for independent development of individual parameter types, each section (such as ``Bonds``, ``Angles``, etc) now has its own version as well (currently all 0.3).
* All units are now stored in expressions with their corresponding values. For example, distances are now stored as ``1.526*angstrom``, instead of storing the unit separately in the section header.
* The current allowed value of the ``potential`` field for ``ProperTorsions`` and ``ImproperTorsions`` tags is no longer ``charmm``, but is rather ``k*(1+cos(periodicity*theta-phase))``.
  It was pointed out to us that CHARMM-style torsions deviate from this formula when the periodicity of a torsion term is 0, and we do not intend to reproduce that behavior.
* SMIRNOFF spec documentation has been updated with tables of keywords and their defaults for each parameter section and parameter type.
  These tables will track the allowed keywords and default behavior as updated versions of individual parameter sections are released.

Performance improvements and bugfixes
"""""""""""""""""""""""""""""""""""""

* `PR #329 <https://github.com/openforcefield/openff-toolkit/pull/329>`_: Performance improvements when creating systems for topologies with many atoms.
* `PR #347 <https://github.com/openforcefield/openff-toolkit/pull/347>`_: Fixes bug in charge assignment that occurs when charges are read from file, and reference and charge molecules have different atom orderings.


New features
""""""""""""

* `PR #311 <https://github.com/openforcefield/openff-toolkit/pull/311>`_: Several new experimental functions.

  * Adds :py:meth:`convert_0_2_smirnoff_to_0_3 <openff.toolkit.utils.utils.convert_0_2_smirnoff_to_0_3>`, which takes a SMIRNOFF 0.2-spec data dict, and updates it to 0.3.
    This function is called automatically when creating a ``ForceField`` from a SMIRNOFF 0.2 spec OFFXML file.
  * Adds :py:meth:`convert_0_1_smirnoff_to_0_2 <openff.toolkit.utils.utils.convert_0_1_smirnoff_to_0_2>`, which takes a SMIRNOFF 0.1-spec data dict, and updates it to 0.2.
    This function is called automatically when creating a ``ForceField`` from a SMIRNOFF 0.1 spec OFFXML file.
  * NOTE: The format of the "SMIRNOFF data dict" above is likely to change significantly in the future.
    Users that require a stable serialized ForceField object should use the output of :py:meth:`ForceField.to_string('XML') <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.to_string>` instead.
  * Adds :py:class:`ParameterHandler <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>` and :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>` :py:meth:`add_cosmetic_attribute <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType.add_cosmetic_attribute>` and :py:meth:`delete_cosmetic_attribute <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType.delete_cosmetic_attribute>` functions.
    Once created, cosmetic attributes can be accessed and modified as attributes of the underlying object (eg. ``ParameterType.my_cosmetic_attrib = 'blue'``)
    These functions are experimental, and we are interested in feedback on how cosmetic attribute handling could be improved. (`See Issue #338 <https://github.com/openforcefield/openff-toolkit/issues/338>`_)
    Note that if a new cosmetic attribute is added to an object without using these functions, it will not be recognized by the toolkit and will not be written out during serialization.
  * Values for the top-level ``Author`` and ``Date`` tags are now kept during SMIRNOFF data I/O.
    If multiple data sources containing these fields are read, the values are concatenated using "AND" as a separator.


API-breaking changes
""""""""""""""""""""
* :py:meth:`ForceField.to_string <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.to_string>` and :py:meth:`ForceField.to_file <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField.to_file>` have had the default value of their ``discard_cosmetic_attributes`` kwarg set to False.
* :py:class:`ParameterHandler <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>` and :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>` constructors now expect the ``version`` kwarg (per the SMIRNOFF spec change above)
  This requirement can be skipped by providing the kwarg ``skip_version_check=True``
* :py:class:`ParameterHandler <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>` and :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>` functions no longer handle ``X_unit`` attributes in SMIRNOFF data (per the SMIRNOFF spec change above).
* The scripts in ``utilities/convert_frosst`` are now deprecated.
  This functionality is important for provenance and will be migrated to the ``openforcefield/smirnoff99Frosst`` repository in the coming weeks.
* :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>` ``._SMIRNOFF_ATTRIBS`` is now :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>` ``._REQUIRED_SPEC_ATTRIBS``, to better parallel the structure of the ``ParameterHandler`` class.
* :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>` ``._OPTIONAL_ATTRIBS`` is now :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>` ``._OPTIONAL_SPEC_ATTRIBS``, to better parallel the structure of the ``ParameterHandler`` class.
* Added class-level dictionaries :py:class:`ParameterHandler <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>` ``._DEFAULT_SPEC_ATTRIBS`` and :py:class:`ParameterType <openff.toolkit.typing.engines.smirnoff.parameters.ParameterType>` ``._DEFAULT_SPEC_ATTRIBS``.

0.3.0 - API Improvements
------------------------

Several improvements and changes to public API.

New features
""""""""""""

* `PR #292 <https://github.com/openforcefield/openff-toolkit/pull/292>`_: Implement ``Topology.to_openmm`` and remove ``ToolkitRegistry.toolkit_is_available``
* `PR #322 <https://github.com/openforcefield/openff-toolkit/pull/322>`_: Install directories for the lookup of OFFXML files through the entry point group ``openforcefield.smirnoff_forcefield_directory``. The ``ForceField`` class doesn't search in the ``data/forcefield/`` folder anymore (now renamed ``data/test_forcefields/``), but only in ``data/``.

API-breaking Changes
""""""""""""""""""""
* `PR #278 <https://github.com/openforcefield/openff-toolkit/pull/278>`_: Standardize variable/method names
* `PR #291 <https://github.com/openforcefield/openff-toolkit/pull/291>`_: Remove ``ForceField.load/to_smirnoff_data``, add ``ForceField.to_file/string`` and ``ParameterHandler.add_parameters``. Change behavior of ``ForceField.register_X_handler`` functions.

Bugfixes
""""""""
* `PR #327 <https://github.com/openforcefield/openff-toolkit/pull/327>`_: Fix units in tip3p.offxml (note that this file is still not loadable by current toolkit)
* `PR #325 <https://github.com/openforcefield/openff-toolkit/pull/325>`_: Fix solvent box for provided test system to resolve periodic clashes.
* `PR #325 <https://github.com/openforcefield/openff-toolkit/pull/325>`_: Add informative message containing Hill formula when a molecule can't be matched in ``Topology.from_openmm``.
* `PR #325 <https://github.com/openforcefield/openff-toolkit/pull/325>`_: Provide warning or error message as appropriate when a molecule is missing stereochemistry.
* `PR #316 <https://github.com/openforcefield/openff-toolkit/pull/316>`_: Fix formatting issues in GBSA section of SMIRNOFF spec
* `PR #308 <https://github.com/openforcefield/openff-toolkit/pull/308>`_: Cache molecule SMILES to improve system creation speed
* `PR #306 <https://github.com/openforcefield/openff-toolkit/pull/306>`_: Allow single-atom molecules with all zero coordinates to be converted to OE/RDK mols
* `PR #313 <https://github.com/openforcefield/openff-toolkit/pull/313>`_: Fix issue where constraints are applied twice to constrained bonds

0.2.2 - Bugfix release
----------------------

This release modifies an example to show how to parameterize a solvated system, cleans up backend code, and makes several improvements to the README.

Bugfixes
""""""""
* `PR #279 <https://github.com/openforcefield/openff-toolkit/pull/279>`_: Cleanup of unused code/warnings in main package ``__init__``
* `PR #259 <https://github.com/openforcefield/openff-toolkit/pull/259>`_: Update T4 Lysozyme + toluene example to show how to set up solvated systems
* `PR #256 <https://github.com/openforcefield/openff-toolkit/pull/256>`_ and `PR #274 <https://github.com/openforcefield/openff-toolkit/pull/274>`_: Add functionality to ensure that links in READMEs resolve successfully


0.2.1 - Bugfix release
----------------------

This release features various documentation fixes, minor bugfixes, and code cleanup.

Bugfixes
""""""""
* `PR #267 <https://github.com/openforcefield/openff-toolkit/pull/267>`_: Add neglected ``<ToolkitAM1BCC>`` documentation to the SMIRNOFF 0.2 spec
* `PR #258 <https://github.com/openforcefield/openff-toolkit/pull/258>`_: General cleanup and removal of unused/inaccessible code.
* `PR #244 <https://github.com/openforcefield/openff-toolkit/pull/244>`_: Improvements and typo fixes for BRD4:inhibitor benchmark

0.2.0 - Initial RDKit support
-----------------------------

This version of the toolkit introduces many new features on the way to a 1.0.0 release.

New features
""""""""""""

* Major overhaul, resulting in the creation of the `SMIRNOFF 0.2 specification <https://open-forcefield-toolkit.readthedocs.io/en/0.2.0/smirnoff.html>`_ and its XML representation
* Updated API and infrastructure for reference SMIRNOFF :class:`ForceField <openff.toolkit.typing.engines.smirnoff.forcefield.ForceField>` implementation
* Implementation of modular :class:`ParameterHandler <openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler>` classes which process the topology to add all necessary forces to the system.
* Implementation of modular :class:`ParameterIOHandler <openff.toolkit.typing.engines.smirnoff.io.ParameterIOHandler>` classes for reading/writing different serialized SMIRNOFF force field representations
* Introduction of :class:`Molecule <openff.toolkit.topology.Molecule>` and :class:`Topology <openff.toolkit.topology.Topology>` classes for representing molecules and biomolecular systems
* New :class:`ToolkitWrapper <openff.toolkit.utils.toolkits.ToolkitWrapper>` interface to RDKit, OpenEye, and AmberTools toolkits, managed by :class:`ToolkitRegistry <openff.toolkit.utils.toolkits.ToolkitRegistry>`
* API improvements to more closely follow `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ guidelines
* Improved documentation and examples

0.1.0
-----

This is an early preview release of the toolkit that matches the functionality described in the preprint describing the SMIRNOFF v0.1 force field format: `[DOI] <https://doi.org/10.1101/286542>`_.

New features
""""""""""""

This release features additional documentation, code comments, and support for automated testing.

Bugfixes
""""""""

Treatment of improper torsions
''''''''''''''''''''''''''''''

A significant (though currently unused) problem in handling of improper torsions was corrected.
Previously, non-planar impropers did not behave correctly, as six-fold impropers have two potential chiralities.
To remedy this, SMIRNOFF impropers are now implemented as three-fold impropers with consistent chirality.
However, current force fields in the SMIRNOFF format had no non-planar impropers, so this change is mainly aimed at future work.

:::
