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

**Note**: The content in this notebook demonstrates usage of `nglview` with `openforcefield`.
This can be tricky to get working.
Install with:

    conda install -c conda-forge nglview
    
And configure for use with Jupyter with:

    jupyter-nbextension enable nglview --py --sys-prefix
    
To use with Jupyterlab, configure with:

    jupyter labextension install  nglview-js-widgets
    jupyter-labextension install @jupyter-widgets/jupyterlab-manager
