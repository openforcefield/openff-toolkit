# Built-in visualizations in `openforcefield`

PR [#560](https://github.com/openforcefield/openforcefield/pull/560) introduced rich representation for the `openforcefield.topology.Molecule` objects. This means you can visualize them in your Jupyter Notebooks:

We have implemented three backends:
- `rdkit`
- `openeye`
- `nglview` (requires conformers)

There are two ways to invoke the visualization:
- implicitly, by _evaluating_ the object in a cell
- explicitly, by calling `Molecule.visualize()`

This notebook demonstrates all the possible uses.