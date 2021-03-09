| **Status** | [![GH Actions Status](https://github.com/openforcefield/openforcefield/workflows/CI/badge.svg)](https://github.com/openforcefield/openforcefield/actions?query=branch%3Amaster+workflow%3ACI)  [![Codecov coverage](https://img.shields.io/codecov/c/github/openforcefield/openforcefield.svg?logo=Codecov&logoColor=white)](https://codecov.io/gh/openforcefield/openforcefield) [![LGTM analysis](https://img.shields.io/lgtm/grade/python/g/openforcefield/openforcefield.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/openforcefield/openforcefield/context:python) |
| :------ | :------- |
| **Latest Release** | [![Last release tag](https://img.shields.io/github/release-pre/openforcefield/openforcefield.svg)](https://github.com/openforcefield/openforcefield/releases)  [![Commits since release](https://img.shields.io/github/commits-since/openforcefield/openforcefield/0.8.2.svg)](https://github.com/openforcefield/openforcefield/releases/tag/0.8.2) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4323602.svg)](https://doi.org/10.5281/zenodo.4323602) |
| **Communication** | [![docs latest](https://img.shields.io/badge/docs-latest-5077AB.svg?logo=read%20the%20docs)](https://open-forcefield-toolkit.readthedocs.io/en/latest/) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/openforcefield/openforcefield/master?filepath=%2Fexamples%2F) [![dev chat on slack](https://img.shields.io/badge/dev_chat-on_slack-808493.svg?logo=slack)](https://join.slack.com/t/openforcefieldgroup/shared_invite/enQtNjA4MTMxMDg0MDAxLWY3Y2Q5NDY4MmU1OTIzMDhiYzFjOWFkZGFjN2Y4N2Q4OTRkOWNjODVhMDMxMzkwMDcxNDA5MjYyNjJjYjE2NTM) |
| **Foundation** | [![license](https://img.shields.io/github/license/openforcefield/openforcefield.svg)](https://opensource.org/licenses/MIT) [![platforms](https://img.shields.io/badge/Platforms-Linux%2C%20MacOS-orange.svg)](https://open-forcefield-toolkit.readthedocs.io/en/latest/installation.html) [![python](https://img.shields.io/badge/python-3.6%2C%203.7-blue.svg)](https://open-forcefield-toolkit.readthedocs.io/en/latest/installation.html) [![Funding](https://img.shields.io/badge/Funding-Open%20Force%20Field%20Consortium-brightgreen.svg)](http://openforcefield.org) |
| **Installation** | [![Releases](https://img.shields.io/badge/obtain-latest-green.svg)](https://github.com/openforcefield/openforcefield/releases) [![Conda](https://img.shields.io/conda/v/omnia/openforcefield.svg)](https://anaconda.org/omnia/openforcefield) [![Last updated](https://anaconda.org/omnia/openforcefield/badges/latest_release_relative_date.svg)](https://anaconda.org/omnia/openforcefield) [![Anaconda Cloud downloads](https://anaconda.org/omnia/openforcefield/badges/downloads.svg)](https://anaconda.org/omnia/openforcefield) |

# The Open Force Field toolkit

The Open Force Field Toolkit, built by the [Open Force Field Initiative](http://openforcefield.org), is a Python toolkit for the development and application of modern molecular mechanics force fields based on direct chemical perception and rigorous statistical parameterization methods.

The toolkit currently covers two main areas we have committed to stably maintain throughout their lifetimes:
* Tools for using the [SMIRKS Native Open Force Field (SMIRNOFF) specification](https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html)
* Tools for [direct chemical environment perception](https://dx.doi.org/10.1021/acs.jctc.8b00640) and manipulation

## Documentation

[Documentation](https://open-forcefield-toolkit.readthedocs.io/en/latest/) for the `openforcefield` toolkit is hosted at [readthedocs](https://open-forcefield-toolkit.readthedocs.io/en/latest). Example notebooks are hosted on [binder](https://mybinder.org/v2/gh/openforcefield/openforcefield/master?filepath=%2Fexamples%2F)

## Installation

`openforcefield` is a Python toolkit, and supports Python 3.6 and 3.7.

### Installing via conda

Detailed installation instructions can be found [here](https://open-forcefield-toolkit.readthedocs.io/en/latest/installation.html).

### Installing from source

If you need to install via source, see the build and run package requirements listed in the [development conda recipe](https://github.com/openforcefield/openforcefield/blob/084-update/devtools/conda-recipe/meta.yaml).

## Force Fields

Two major force field development efforts have been undertaken by the Initiative, with results hosted in separate repositories.

* The [Open Force Fields repository](https://github.com/openforcefield/openforcefields/), which features the [Parsley force field line](https://openforcefield.org/community/news/general/introducing-openforcefield-1.0/). This is the Open Force Field Initiative's first effort toward building a _new_ force field. The initial parameters are taken from smirnoff99Frosst, but software and data produced by the Initiative's efforts are being used to refit parameter values and add new SMIRKS-based parameters.
* The [smirnoff99Frosst repository](https://github.com/openforcefield/smirnoff99Frosst/), which is descended from AMBER's parm99 forcefield as well as Merck-Frosst's parm@frosst. This line of force fields does not aim to alter parameter values, but is instead a test of accurately converting an atom type-based force field to the SMIRNOFF format.

Force fields from both of these packages are available in their respective GitHub repositories and also as conda packages. Tables detailing the individual file names/versions within these force field lines are in the README of each repository. By default, installing the Open Force Field toolkit using `conda` or the single-file toolkit installers will also install these conda packages. A plugin architecture is provided for other force field developers to produce python/conda packages that can be imported by the Open Force Field Toolkit as well.

# Toolkit features

## The SMIRKS Native Open Force Field (SMIRNOFF) format

This repository provides tools for using the [SMIRKS Native Open Force Field (SMIRNOFF) specification](https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html), which currently supports an XML representation for force field definition files.

By convention, files containing [XML representations](https://en.wikipedia.org/wiki/XML) of SMIRNOFF force fields carry `.offxml` extensions.

Example SMIRNOFF `.offxml` force field definitions can be found in [`openforcefield/data/test_forcefields/`](https://github.com/openforcefield/openforcefield/tree/084-update/openforcefield/data/test_forcefields). These force fields are for testing only, and we neither record versions of these files, nor do we guarantee their correctness or completeness.

### Working with SMIRNOFF parameter sets

The SMIRNOFF `ForceField` class is essentially a drop-in replacement for the [OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField).

```python
# Load a molecule into the openforcefield Molecule object
from openforcefield.topology import Molecule
from openforcefield.utils import get_data_file_path
sdf_file_path = get_data_file_path('molecules/ethanol.sdf')
molecule = Molecule.from_file(sdf_file_path)

# Create an openforcefield Topology object from the molecule
from openforcefield.topology import Topology
topology = Topology.from_molecules(molecule)

# Load the latest OpenFF force field definition
from openforcefield.typing.engines.smirnoff import ForceField
forcefield = ForceField('openff-1.2.0.offxml')

# Create an OpenMM system representing the molecule with SMIRNOFF-applied parameters
openmm_system = forcefield.create_openmm_system(topology)

```
Detailed examples of using SMIRNOFF with the toolkit can be found [in the documentation](https://open-forcefield-toolkit.readthedocs.io/en/latest/examples.html).

# Frequently asked questions (FAQ)

See [`FAQ.md`](FAQ.md) for answers to a variety of common problems, such as:
* Why do I need to provide molecules corresponding to the components of my system, or a `Topology` with bond orders?
* Can I use an AMBER, CHARMM, or gromacs topology/coordinate file as a starting point for applying a SMIRNOFF force field?
* What if I am starting from a PDB file?

# Contributors

For a full list of contributors, see the [GitHub Contributors page](https://github.com/openforcefield/openforcefield/graphs/contributors).
