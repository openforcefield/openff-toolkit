[![Build Status](https://travis-ci.org/openforcefield/openforcefield.svg?branch=master)](https://travis-ci.org/openforcefield/openforcefield?branch=master)
[![Documentation](https://readthedocs.org/projects/open-forcefield-toolkit/badge/?version=latest)](http://open-forcefield-toolkit.readthedocs.io/en/latest/?badge=latest)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/openforcefield/openforcefield.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/openforcefield/openforcefield/context:python)
[![codecov](https://codecov.io/gh/openforcefield/openforcefield/branch/master/graph/badge.svg)](https://codecov.io/gh/openforcefield/openforcefield)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Install with Conda](https://anaconda.org/omnia/openforcefield/badges/installer/conda.svg)](https://anaconda.org/omnia/openforcefield)
[![Anaconda Cloud platforms](https://anaconda.org/omnia/openforcefield/badges/platforms.svg)](https://anaconda.org/omnia/openforcefield)
[![Anaconda Cloud downloads](https://anaconda.org/omnia/openforcefield/badges/downloads.svg)](https://anaconda.org/omnia/openforcefield)
[![Funding](https://img.shields.io/badge/Funding-Open%20Force%20Field%20Consortium-brightgreen.svg)](http://openforcefield.org)

# The Open Force Field toolkit

The Open Force Field Toolkit, built by the [Open Force Field Initiative](http://openforcefield.org), is a Python toolkit for the development and application of modern molecular mechanics force fields based on direct chemical perception and rigorous statistical parameterization methods.

The toolkit currently covers two main areas we have committed to stably maintain throughout their lifetimes:
* Tools for using the [SMIRKS Native Open Force Field (SMIRNOFF) specification](https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html)
* Tools for [direct chemical environment perception](https://dx.doi.org/10.1021/acs.jctc.8b00640) and manipulation

## Documentation

[Documentation](https://open-forcefield-toolkit.readthedocs.io/en/latest/) for the `openforcefield` toolkit is hosted at [readthedocs](https://open-forcefield-toolkit.readthedocs.io/en/latest).

## Installation

`openforcefield` is a Python toolkit, and supports Python 3.6 and 3.7.

### Installing via conda

Detailed installation instructions can be found [here](https://open-forcefield-toolkit.readthedocs.io/en/latest/installation.html).

### Installing from source

If you need to install via source, see the build and run package requirements listed in the [development conda recipe](https://github.com/openforcefield/openforcefield/blob/master/devtools/conda-recipe/meta.yaml).

# Toolkit features

## The SMIRKS Native Open Force Field (SMIRNOFF) format

This repository provides tools for using the [SMIRKS Native Open Force Field (SMIRNOFF) specification](https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html), which currently supports an XML representation for force field definition files.

By convention, files containing [XML representations](https://en.wikipedia.org/wiki/XML) of SMIRNOFF force fields carry `.offxml` extensions.

Example SMIRNOFF `.offxml` force field definitions can be found in [`openforcefield/data/forcefield`](https://github.com/openforcefield/openforcefield/tree/master/openforcefield/data/forcefield).

### Working with SMIRNOFF parameter sets

The SMIRNOFF `ForceField` class is essentially a drop-in replacement for the [OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField).

```python
# Load a molecule into the openforcefield Molecule object
from openforcefield.topology import Molecule
from openforcefield.utils import get_data_filename
sdf_file_path = get_data_filename('molecules/ethanol.sdf')
molecule = Molecule.from_file(sdf_file_path)

# Create an openforcefield Topology object from the molecule
from openforcefield.topology import Topology
topology = Topology.from_molecules(molecule)

# Load the smirnoff99Frosst SMIRNOFF force field definition
from openforcefield.typing.engines.smirnoff import ForceField
forcefield = ForceField('smirnoff99Frosst.offxml')

# Create an OpenMM system representing the molecule with SMIRNOFF-applied parameters
openmm_system = forcefield.create_openmm_system(topology)

# Load a SMIRNOFF small molecule forcefield for alkanes, ethers, and alcohols
forcefield = ForceField('Frosst_AlkEthOH_parmAtFrosst.offxml')
```
Detailed examples of using SMIRNOFF with the toolkit can be found [in the documentation](https://open-forcefield-toolkit.readthedocs.io/en/latest/examples.html).

## Chemical environments: Chemical environment perception and manipulation

The `ChemicalEnvironments` class can be used to parse and manipulate [tagged SMARTS strings](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html) or single-fragment [SMIRKS strings](http://www.daylight.com/dayhtml/doc/theory/theory.smirks.html) representing chemical environments with tagged atoms.

### Working with chemical environments

```python
from openforcefield.typing.chemistry import environment

smirks = "[#6X3,#7:1]~;@[#8;r:2]~;@[#6X3,#7:3]"
angle = environment.AngleChemicalEnvironment(smirks=smirks)
print(angle.asSMIRKS())
# "[#6X3,#7:1]~;@[#8;r:2]~;@[#6X3,#7:3]"

# add a new atom
atom3 = angle.selectAtom(3)
alpha_ORtypes = [('#8', ['X2'])]
alpha_bondANDtypes = ['!@']
alpha = angle.addAtom(atom3, bondANDtypes = alpha_bondANDtypes, newORtypes = alpha_ORtypes)
print(alpha.asSMIRKS()) # smirks for atom only
# "[#8X2H1;R0]"
print(angle.asSMIRKS())
# "[#6X3,#7:1]~;@[#8;r:2]~;@[#6X3,#7:3]~;!@[#8X2]"
```
Daylight provides detailed specifications of the [SMILES](http://www.daylight.com/dayhtml_tutorials/languages/smiles/index.html), [SMARTS](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html), and [SMIRKS](http://www.daylight.com/dayhtml_tutorials/languages/smirks/index.html) languages.

# Frequently asked questions (FAQ)

See [`FAQ.md`](FAQ.md) for answers to a variety of common problems, such as:
* Why do I need to provide molecules corresponding to the components of my system, or a `Topology` with bond orders?
* Can I use an AMBER, CHARMM, or gromacs topology/coordinate file as a starting point for applying a SMIRNOFF force field?
* What if I am starting from a PDB file?

# Contributors

For a full list of contributors, see the [GitHub Contributors page](https://github.com/openforcefield/openforcefield/graphs/contributors).
