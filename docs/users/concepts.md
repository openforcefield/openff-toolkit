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

OpenFF [`Topology`](openff.toolkit.topology.Topology)
: A collection of molecules, as well as optional box vector, positions, and other information. A `Topology` describes the chemistry of a system, but has no force field parameters.

OpenFF [`ForceField`](openff.toolkit.typing.engines.smirnoff.forcefield.ForceField)
: An object generated from an OFFXML file (or other source of SMIRNOFF data).
  Most information from the SMIRNOFF data source is stored in this object's several `ParameterHandler`s, however some top-level SMIRNOFF data is stored in the `ForceField` object itself.

OpenFF [`Interchange`](openff.interchange:index)
: A `Topology` that has been parametrized by a `ForceField`. An `Interchange` contains all the information needed to calculate an energy or start a simulation. `Interchange` also provides methods to export this information
to MD simulation engines.
