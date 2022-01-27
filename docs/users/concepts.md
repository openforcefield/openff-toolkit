# Core concepts

OpenFF [`Molecule`](openff.toolkit.topology.Molecule)
: A graph representation of a molecule containing enough information to unambiguously parametrize it.
  Required data fields for a `Molecule` are:

  - `atoms`: element (integer), formal_charge (integer), is_aromatic (boolean), stereochemistry (`'R'`/`'S'`/`None`)
  - `bonds`: order (integer), is_aromatic (boolean), stereochemistry (`'E'`/`'Z'`/`None`)

  There are several other optional attributes such as `conformers` and `partial_charges` that may be populated in the `Molecule` data structure.
  These are considered "optional" because they are not required for system creation, however if those fields are populated, the user MAY use them to override values that would otherwise be generated during system creation.

  A dictionary, `Molecule.properties` is exposed, which is a Python dict that can be populated with arbitrary data.
  This data should be considered cosmetic and should not affect system creation.
  Whenever possible, molecule serialization or format conversion should preserve this data.

OpenFF `System`
: An object that contains everything needed to calculate a molecular system's energy, except the atomic coordinates.
  Note that this does not exist yet, and that OpenMM `System` objects are being used for this purpose right now. Development is underway on [GitHub](https://github.com/openforcefield/openff-system).

OpenFF [`Topology`](openff.toolkit.topology.Topology)
: An object that efficiently holds many OpenFF `Molecule` objects.
  The atom indexing in a `Topology` may differ from those of the underlying `Molecule`s

OpenFF [`TopologyMolecule`](openff.toolkit.topology.TopologyMolecule)
: The efficient data structures that make up an OpenFF `Topology`.
  There is one `TopologyMolecule` for each instance of a chemical species in a `Topology`.
  However, each unique chemical species has a single OpenFF `Molecule` representing it, which may be shared by multiple `TopologyMolecules`.
  `TopologyMolecule`s contain an atom index map, as several copies of the same chemical species in a `Topology` may be present with different atom orderings.
  This data structure allows the OpenFF toolkit to only parametrize each unique `Molecule` once, and then write a copy of the assigned parameters out for each of the `Molecule` in the `Topology` (accounting for atom indexing differences in the process).

OpenFF [`ForceField`](openff.toolkit.typing.engines.smirnoff.forcefield.ForceField)
: An object generated from an OFFXML file (or other source of SMIRNOFF data).
  Most information from the SMIRNOFF data source is stored in this object's several `ParameterHandler`s, however some top-level SMIRNOFF data is stored in the `ForceField` object itself.
