[![Build Status](https://travis-ci.org/open-forcefield-group/smirnoff.svg?branch=master)](https://travis-ci.org/open-forcefield-group/smirnoff?branch=master)

# Open Forcefield Group toolkit

This repository contains a number of tools from the [Open Force Field Group](http://github.com/open-forcefield-group) for the development and use of modern molecular mechanics forcefields based on direct chemical perception and parameterized with rigorous statistical methods.

This repository hosts tools that we have committed to stably maintain throughout their lifetimes:
* The [SMIRks Native Open Force Field (SMIRNOFF)](https://github.com/open-forcefield-group/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md) direct chemical perception forcefield specification language
* Tools for direct chemical environment perception and manipulation

## Installation

We recommend the [miniconda](http://conda.pydata.org/miniconda.html) Python distribution.
To install `miniconda` on `osx` with `bash`, this is:
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda3
export PATH="$HOME/miniconda3/bin:${PATH}"
```
These tools currently require the OpenEye toolkit (which requires a [license](https://www.eyesopen.com/licensing-philosophy) that is free for academics indenting to rapidly release results into the public domain):
```bash
pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits
```
Install `openforcefield` tools via conda:
```bash
conda install --yes -c omnia openforcefield
```

# Tools

## `SMIRNOFF`: SMIRks Native Open Force Field

This repository houses the SMIRNOFF SMIRKS-based force field format, along with classes to parameterize OpenMM systems given [SMIRNOFF `.ffxml` format files](https://github.com/open-forcefield-group/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md).

The SMIRNOFF force field format is documented [here](https://github.com/open-forcefield-group/smirnoff/blob/master/The-SMIRNOFF-force-field-format.md).

The SMIRNOFF forcefield format is available in sample form under `data/forcefield`, and is handled by `forcefield.py`.
 An example comparing SMIRNOFF versus AMBER energies for the parm@frosst forcefield is provided under
examples/SMIRNOFF_comparison, where two scripts can compare energies for a single molecule or for the entire AlkEthOH set.
Note that two forcefields are currently available in this format, `Frosst_AlkEtOH.ffxml`,
the parm@frosst forcefield as it should have been for this set, and `Frosst_AlkEtOH_parmAtFrosst.ffxml`,
the forcefield as it was actually implemented (containing several bugs as noted in the file itself).

It can also be of interest to know what SMIRNOFF parameters would be applied to particular molecules. Utility functionality for this is provided under `forcefield_labeler.py`, which has generally similar structure to `forcefield.py` but instead of providing OpenMM systems with parameters, it can be applied to specific molecules and returns information about what parameters would be applied.

## `ChemicalEnvironment`: Tools for chemical environment perception and manipulation

Documentation forthcoming

# Manifest

* `examples/` - some examples - look here to get started
* `openforcefield/` - openforcefield tools
* `devtools/` - continuous integration and packaging scripts and utilities
* `utilities/` - utilities; scripts to convert parm@frosst modified `frcmod` files to SMIRNOFF XML
* `oe_license.txt.enc` - encrypted OpenEye license for continuous integration testing
* `.travis.yml` - travis-ci continuous integration file

# Contributors

* [David L. Mobley (UCI)](https://github.com/davidlmobley)
* [John D. Chodera (MSKCC)](https://github.com/jchodera)
* [Caitlin Bannan (UCI)](https://github.com/bannanc)
* [Camila Zanette (UCI)](https://github.com/camizanette)
* [Christopher I. Bayly (OpenEye)](https://github.com/cbayly13)
* [Nathan M. Lim (UCI)](https://github.com/nathanmlim)
