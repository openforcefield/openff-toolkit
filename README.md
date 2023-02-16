| **Status** | [![GH Actions Status](https://github.com/openforcefield/openff-toolkit/workflows/CI/badge.svg)](https://github.com/openforcefield/openff-toolkit/actions?query=branch%3Amain+workflow%3ACI)  [![Codecov coverage](https://img.shields.io/codecov/c/github/openforcefield/openff-toolkit.svg?logo=Codecov&logoColor=white)](https://codecov.io/gh/openforcefield/openff-toolkit) |
| :------ |:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Latest Release** | [![Last release tag](https://img.shields.io/github/release-pre/openforcefield/openff-toolkit.svg)](https://github.com/openforcefield/openff-toolkit/releases)  [![Commits since release](https://img.shields.io/github/commits-since/openforcefield/openff-toolkit/0.11.4.svg)](https://github.com/openforcefield/openff-toolkit/releases/tag/0.11.4)  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7506404.svg)](https://doi.org/10.5281/zenodo.7506404)
| **Communication** | [![docs stable](https://img.shields.io/badge/docs-stable-5077AB.svg?logo=read%20the%20docs)](https://open-forcefield-toolkit.readthedocs.io/) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/openforcefield/openff-toolkit/stable?filepath=%2Fexamples%2F) [![user & dev discussions](https://img.shields.io/badge/user%20%26%20dev%20discussions-GitHub-red?logo=github)](https://github.com/openforcefield/discussions/discussions) |
| **Foundation** | [![license](https://img.shields.io/github/license/openforcefield/openff-toolkit.svg)](https://opensource.org/licenses/MIT) [![platforms](https://img.shields.io/badge/Platforms-Linux%2C%20MacOS-orange.svg)](https://open-forcefield-toolkit.readthedocs.io/en/stable/installation.html) [![python](https://img.shields.io/badge/python-3.8%2C%203.9-blue.svg)](https://open-forcefield-toolkit.readthedocs.io/en/stable/installation.html) [![Funding](https://img.shields.io/badge/Funding-Open%20Force%20Field%20Consortium-brightgreen.svg)](http://openforcefield.org)        |
| **Installation** | [![Releases](https://img.shields.io/badge/obtain-latest-green.svg)](https://github.com/openforcefield/openff-toolkit/releases) [![Conda](https://img.shields.io/conda/v/conda-forge/openff-toolkit.svg)](https://anaconda.org/conda-forge/openff-toolkit) [![Last updated](https://anaconda.org/conda-forge/openff-toolkit/badges/latest_release_relative_date.svg)](https://anaconda.org/conda-forge/openff-toolkit) [![Anaconda Cloud downloads](https://anaconda.org/conda-forge/openff-toolkit/badges/downloads.svg)](https://anaconda.org/conda-forge/openff-toolkit)          |

# The Open Force Field toolkit

The Open Force Field Toolkit, built by the [Open Force Field Initiative](http://openforcefield.org), is a Python toolkit for the development and application of modern molecular mechanics force fields based on direct chemical perception and rigorous statistical parameterization methods.

The toolkit currently covers two main areas we have committed to stably maintain throughout their lifetimes:
* Tools for using the [SMIRKS Native Open Force Field (SMIRNOFF) specification](https://openforcefield.github.io/standards/standards/smirnoff/)
* Tools for [direct chemical environment perception](https://dx.doi.org/10.1021/acs.jctc.8b00640) and manipulation

**Note**: Prior to version 0.9.0, this toolkit and its associated repository were named `openforcefield` and used different import paths. For details on this change and migration instructions, see the [release notes](https://open-forcefield-toolkit.readthedocs.io/en/stable/releasehistory.html#namespace-migration) of version 0.9.0.

## Documentation

[Documentation](https://open-forcefield-toolkit.readthedocs.io/en/stable/) for the Open Force Field Toolkit is hosted at [readthedocs](https://open-forcefield-toolkit.readthedocs.io/en/stable). Example notebooks are hosted on [binder](https://mybinder.org/v2/gh/openforcefield/openff-toolkit/stable?filepath=%2Fexamples%2F)

## How to cite

Please cite the OpenFF Toolkit using the [Zenodo record](https://zenodo.org/record/7506404) of the [latest release](https://zenodo.org/record/7506404) or the version that was used. The BibTeX reference of the latest release can be found [here](https://zenodo.org/record/7506404/export/hx#.Y7cYbOzMKrM).

## Installation

The Open Force Field Toolkit (`openff-toolkit`) is a Python toolkit, and supports Python 3.8 through 3.10.

### Installing via conda

Detailed installation instructions can be found [here](https://open-forcefield-toolkit.readthedocs.io/en/stable/installation.html).

### Installing from source

If you need to install via source, see the build and run package requirements listed in the [development conda recipe](https://github.com/openforcefield/openff-toolkit/blob/0.8.3/devtools/conda-recipe/meta.yaml).

## Force Fields

Two major force field development efforts have been undertaken by the Initiative, with results hosted in separate repositories.

* The [Open Force Fields repository](https://github.com/openforcefield/openff-forcefields/), which features the [Parsley](https://openforcefield.org/community/news/general/introducing-openforcefield-1.0/) and [Sage](https://openforcefield.org/community/news/general/sage2.0.0-release/) force field lines. These are the Open Force Field Initiative's efforts toward building _new_ force fields. The initial parameters are taken from smirnoff99Frosst, but software and data produced by the Initiative's efforts have been used to refit parameter values and add new SMIRKS-based parameters.
* The [smirnoff99Frosst repository](https://github.com/openforcefield/smirnoff99Frosst/), which is descended from AMBER's parm99 force field as well as Merck-Frosst's parm@frosst. This line of force fields does not aim to alter parameter values, but is instead a test of accurately converting an atom type-based force field to the SMIRNOFF format.

Force fields from both of these packages are available in their respective GitHub repositories and also as conda packages. Tables detailing the individual file names/versions within these force field lines are in the README of each repository. By default, installing the Open Force Field toolkit using `conda` or the single-file toolkit installers will also install these conda packages. A [plugin architecture](https://github.com/openforcefield/openff-toolkit/blob/main/FAQ.md#how-can-i-distribute-my-own-force-fields-in-smirnoff-format) is provided for other force field developers to produce python/conda packages that can be imported by the Open Force Field Toolkit as well.

# Toolkit features

## The SMIRKS Native Open Force Field (SMIRNOFF) format

This repository provides tools for using the [SMIRKS Native Open Force Field (SMIRNOFF) specification](https://openforcefield.github.io/standards/standards/smirnoff/), which currently supports an XML representation for force field definition files.

By convention, files containing [XML representations](https://en.wikipedia.org/wiki/XML) of SMIRNOFF force fields carry `.offxml` extensions.

Example SMIRNOFF `.offxml` force field definitions can be found in [`openff/toolkit/data/test_forcefields/`](https://github.com/openforcefield/openff-toolkit/tree/main/openff/toolkit/data/test_forcefields). These force fields are for testing only, and we neither record versions of these files, nor do we guarantee their correctness or completeness.

### Working with SMIRNOFF parameter sets

The `ForceField` class in the OpenFF Toolkit is essentially a drop-in replacement for the [OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField).

```python
# Load a molecule into the OpenFF Molecule object
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import get_data_file_path
sdf_file_path = get_data_file_path('molecules/ethanol.sdf')
molecule = Molecule.from_file(sdf_file_path)

# Create an OpenFF Topology object from the molecule
from openff.toolkit.topology import Topology
topology = Topology.from_molecules(molecule)

# Load the latest OpenFF force field release: version 2.0.0, codename "Sage"
from openff.toolkit.typing.engines.smirnoff import ForceField
forcefield = ForceField('openff-2.0.0.offxml')

# Create an OpenMM system representing the molecule with SMIRNOFF-applied parameters
openmm_system = forcefield.create_openmm_system(topology)

```
Detailed examples of using SMIRNOFF with the toolkit can be found [in the documentation](https://open-forcefield-toolkit.readthedocs.io/en/stable/examples.html).

# Frequently asked questions (FAQ)

See [`FAQ.md`](FAQ.md) for answers to a variety of common problems, such as:
* Why do I need to provide molecules corresponding to the components of my system, or a `Topology` with bond orders?
* Can I use an AMBER, CHARMM, or gromacs topology/coordinate file as a starting point for applying a SMIRNOFF force field?
* What if I am starting from a PDB file?

# Contributors

For a full list of contributors, see the [GitHub Contributors page](https://github.com/openforcefield/openff-toolkit/graphs/contributors).
