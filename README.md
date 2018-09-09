[![Build Status](https://travis-ci.org/openforcefield/openforcefield.svg?branch=master)](https://travis-ci.org/openforcefield/openforcefield?branch=master)
[![Documentation](https://readthedocs.org/projects/open-forcefield-toolkit/badge/?version=topology)](http://open-forcefield-toolkit.readthedocs.io/en/topology/?badge=topology)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Install with Conda](https://anaconda.org/omnia/openforcefield/badges/installer/conda.svg)](https://anaconda.org/omnia/openforcefield)
[![Anaconda Cloud platforms](https://anaconda.org/omnia/openforcefield/badges/platforms.svg)](https://anaconda.org/omnia/openforcefield)
[![Anaconda Cloud downloads](https://anaconda.org/omnia/openforcefield/badges/downloads.svg)](https://anaconda.org/openforcefield/openforcefield)
[![Funding](https://img.shields.io/badge/Funding-Open%20Force%20Field%20Consortium-brightgreen.svg)](http://openforcefield.org)

# The Open Force Field toolkit

This repository contains a number of tools from the [Open Force Field Initiative](http://openforcefield.org) for the development and use of modern molecular mechanics forcefields based on direct chemical perception and parameterized with rigorous statistical methods.

This repository hosts tools that we have committed to stably maintain throughout their lifetimes:
* Tools for using the [SMIRKS Native Open Force Field (SMIRNOFF) format](https://open-forcefield-toolkit.readthedocs.io/en/topology/smirnoff.html)
* Tools for direct chemical environment perception and manipulation

## Documentation

[Documentation](https://open-forcefield-toolkit.readthedocs.io/en/topology) for the `openforcefield` toolkit is hosted at [readthedocs](https://open-forcefield-toolkit.readthedocs.io/en/topology).

## Installation

`openforcefield` is a Python toolkit, and supports Python 2.7, 3.5, and 3.6.

### Installing via conda

We recommend the this project be used with the [miniconda](http://conda.pydata.org/miniconda.html) Python distribution for automatic installation of dependencies.

To install `miniconda` on `linux`:
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:${PATH}"
```
or, on `osx`:
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:${PATH}"
```

Next, configure your installation to add the `omnia` and `conda-forge` channels:
```bash
conda config --add channels omnia --add channels conda-forge
conda update --yes --all
```

Finally, install `openforcefield` toolkit:
```bash
conda install --yes openforcefield
```

### Optional dependencies

This toolkit can optionally make use of the OpenEye toolkit (which requires a [license](https://www.eyesopen.com/licensing-philosophy) that is free for academics intending to release results into the public domain):
```bash
conda install --yes -c openeye openeye-toolkits
```
Currently, the OpenEye toolkit provides features for generating AM1-BCC charges.

### Installing via source

If you need to install via source, see the build and run package requirements listed in the [development conda recipe](https://github.com/openforcefield/openforcefield/blob/topology/devtools/conda-recipe/meta.yaml).

# Features of the `openforcefield` toolkit

## SMIRNOFF: The SMIRKS Native Open Force Field format

This repository provides tools for using the [SMIRKS Native Open Force Field (SMIRNOFF) format](https://open-forcefield-toolkit.readthedocs.io/en/topology/smirnoff.html), which currently supports an XML representation for force field definition files.

By convention, files containing XML representations of SMIRNOFF force fields carry `.offxml` extensions.

Example SMIRNOFF `.offxml` force field definitions can be found in `openforcefield/data/forcefield`.

Detailed examples of using SMIRNOFF with the toolkit can be found [in the documentation](https://open-forcefield-toolkit.readthedocs.io/en/topology/examples.html).

### Using SMIRNOFF

The SMIRNOFF `ForceField` class is essentially a drop-in replacement for the [OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField).

```python
# Load a molecule into the openforcefield Molecule object
from openforcefield.topology import Molecule
molecule = Molecule.from_file(mol_filename)

# Create an openforcefield Topology object from the molecule
from openforcefield.topology import Topology
topology = Topology.from_molecules(molecule)

# Load the smirnoff99Frosst SMIRNOFF force field definition
from openforcefield.typing.engines.smirnoff import ForceField
forcefield = ForceField('smirnoff99Frosst.offxml')

# Create an OpenMM system representing the molecule with SMIRNOFF-applied parameters
openmm_system = forcefield.create_openmm_system(topology)

# Load a SMIRNOFF small molecule forcefield for alkanes, ethers, and alcohols
forcefield = ForceField(offxml_filename)
```

## `ChemicalEnvironment`: Tools for chemical environment perception and manipulation

The `ChemicalEnvironments` class can be used to parse and manipulate [tagged SMARTS strings](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html) or single-fragment [SMIRKS strings](http://www.daylight.com/dayhtml/doc/theory/theory.smirks.html) representing chemical environments with tagged atoms.

### Example usage

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
If you are not familiar with the SMIRKS language, take a look at these Daylight resources:
* [SMILES](http://www.daylight.com/dayhtml_tutorials/languages/smiles/index.html)
* [SMARTS](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
* [SMIRKS](http://www.daylight.com/dayhtml_tutorials/languages/smirks/index.html)

For more detailed examples, see the README and the [`examples/chemicalEnvironments/using_environment.ipynb`](https://github.com/openforcefield/openforcefield/blob/master/examples/chemicalEnvironments/using_environments.ipynb) example notebook.

# Frequently asked questions (FAQs)

See [`Chemical-starting-points.md`](Chemical-starting-points.md) for answers to a variety of common problems, such as:
* Why do I need to provide molecules corresponding to the components of my system, or a `Topology` with bond orders?
* Can I use an AMBER, CHARMM, or gromacs topology/coordinate file as a starting point for applying a SMIRNOFF force field?
* What if I am starting from a PDB file?

# Manifest

* `examples/` - some examples - look here to get started; see especially `host_guest_simulation` for a detailed worked example of SMIRNOFF simulation of host-guest binding.
* `openforcefield/` - openforcefield tools
* `devtools/` - continuous integration and packaging scripts and utilities
* `utilities/` - utilities; scripts to convert parm@frosst modified `frcmod` files to SMIRNOFF XML
* `oe_license.txt.enc` - encrypted OpenEye license for continuous integration testing
* `.travis.yml` - travis-ci continuous integration file
* `The-SMIRNOFF-force-field-format.md` - specifications for the SMIRNOFF force field format
* `Chemical-starting-points.md` - discussion of appropriate starting points for applying a SMIRNOFF force field to a system

# Contributors

* [David L. Mobley (UCI)](https://github.com/davidlmobley)
* [John D. Chodera (MSKCC)](https://github.com/jchodera)
* [Caitlin Bannan (UCI)](https://github.com/bannanc)
* [Camila Zanette (UCI)](https://github.com/camizanette)
* [Christopher I. Bayly (OpenEye)](https://github.com/cbayly13)
* [Nathan M. Lim (UCI)](https://github.com/nathanmlim)
