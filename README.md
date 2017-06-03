[![Build Status](https://travis-ci.org/open-forcefield-group/openforcefield.svg?branch=master)](https://travis-ci.org/open-forcefield-group/openforcefield?branch=master)

# Open Forcefield Group toolkit

This repository contains a number of tools from the [Open Force Field Group](http://github.com/open-forcefield-group) for the development and use of modern molecular mechanics forcefields based on direct chemical perception and parameterized with rigorous statistical methods.

This repository hosts tools that we have committed to stably maintain throughout their lifetimes:
* The [SMIRks Native Open Force Field (SMIRNOFF)](https://github.com/open-forcefield-group/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md) direct chemical perception forcefield specification language
* Tools for direct chemical environment perception and manipulation

## Installation
We currently support Python 2.7, 3.5 and 3.6.

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
conda install --yes -c conda-forge -c omnia openforcefield
```

# Tools

## `SMIRNOFF`: SMIRks Native Open Force Field

This repository houses the SMIRNOFF SMIRKS-based force field format, along with classes to parameterize OpenMM systems given [SMIRNOFF `.ffxml` format files](https://github.com/open-forcefield-group/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md).

The SMIRNOFF force field format is documented [here](https://github.com/open-forcefield-group/smirnoff/blob/master/The-SMIRNOFF-force-field-format.md).

The SMIRNOFF forcefield format is available in sample form under `data/forcefield`, and is handled by `forcefield.py`.
 An example comparing SMIRNOFF versus AMBER energies for the parm@frosst forcefield is provided under
examples/SMIRNOFF_comparison, where two scripts can compare energies for a single molecule or for the entire AlkEthOH set.
Note that two forcefields are currently available in this format, `Frosst_AlkEthOH.ffxml`,
the parm@frosst forcefield as it should have been for this set, and `Frosst_AlkEthOH_parmAtFrosst.ffxml`,
the forcefield as it was actually implemented (containing several bugs as noted in the file itself).

It can also be of interest to know what SMIRNOFF parameters would be applied to particular molecules. Utility functionality for this is provided under `forcefield_labeler.py`, which has generally similar structure to `forcefield.py` but instead of providing OpenMM systems with parameters, it can be applied to specific molecules and returns information about what parameters would be applied.

### Example usage

The SMIRNOFF `ForceField` class is essentially a drop-in replacement for the [OpenMM `ForceField` class](http://docs.openmm.org/7.1.0/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField), with the additional requirement that an OpenEye `OEMol`-compatible object must also be provided to allow for chemical environment perception (and optionally charges).
For example, if we have an `OEMol` named `mol`, we can create an OpenMM `System` object with the following code:
```python
# Import the SMIRNOFF forcefield engine and some useful tools
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils import read_molecules, get_data_filename, generateTopologyFromOEMol 

# read in molecule from file in openforcefield/data/molecules/
mols = read_molecules('benzene.mol2')

# Get positions and topology in OpenMM-compatible format
topology = generateTopologyFromOEMol(mols[0])

# Load a SMIRNOFF small molecule forcefield for alkanes, ethers, and alcohols
FF_filename = get_data_filename('forcefield/Frosst_AlkEthOH_parmAtFrosst.ffxml')
forcefield = ForceField(FF_filename)

# Create the OpenMM system, additionally specifying a list of OEMol objects for the unique molecules in the system
system = forcefield.createSystem(topology, mols)
```
See `examples/SMIRNOFF_simulation/` for a complete example of how SMIRNOFF can be used for small molecule vacuum simulations, and `examples/mixedFF_structure` for how to set up a system which uses an AMBER forcefield (in this case, AMBER99SB-ILDN) for a protein in combination with SMIRNOFF for a small molecules. Via ParmEd, this can be translated into GROMACS, AMBER, or CHARMM formats for use elsewhere (and additional formats via InterMol).

## `ChemicalEnvironment`: Tools for chemical environment perception and manipulation

ChemicalEnvironments are a python class used to parse and manipulate SMIRKS strings. 
They were created with the goal of being able to automatically sample over chemical perceptions space. 
Someday they will be used to generate SMIRKS patterns for SMIRKS Native-Open Force Fields parameters. 
These are initiated with SMIRKS strings for single molecules fragements`*` and then the information is stored for each atom and bond in the initial fragment. 

`*` NOTE SMIRKS can be used to show how a reaction would happen between fragments in different molecules. This is done with `'.'` between molecules and `'>>'` to indicate a reaction. Chemical Environments can only parse SMIRKS strings for fragments of a single molecule.  
 
```python
from openforcefield.typing.chemistry import environment

smirks = "[#6X3,#7:1]~;@[#8;r:2]~;@[#6X3,#7:3]"
angle = environment.AngleChemicalEnvironment(smirks = smirks)
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

For more detailed examples see README and `using_environment.ipynb` in  `examples/chemicalEnvironments/` 

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
