(installation)=

# Installation

(installation/mamba)=

## Installing via `mamba`

The simplest way to install the Open Force Field Toolkit is via [Mamba](https://github.com/mamba-org/mamba), a drop-in replacement for the [Conda](https://docs.conda.io/en/latest/) package manager.
We publish [packages](https://github.com/conda-forge/openff-toolkit-feedstock) via [`conda-forge`](https://conda-forge.org/).
With Mamba installed, use it to install the OpenFF Toolkit into a new environment:

```shell-session
$ mamba create -n openff-toolkit -c conda-forge openff-toolkit
```

To use the new environment in a shell session, you must first activate it:

```shell-session
$ mamba activate openff-toolkit
```

If you do not have Mamba or Conda installed, see the [ecosystem installation documentation].

:::{note}
Installation via the Mamba package manager is the recommended method since all dependencies are automatically fetched and installed for you.
:::

[ecosystem installation documentation]: inv:openff.docs#install

(installation/platforms)=

### OS support

The OpenFF Toolkit is pure Python, and we expect it to work on any platform that supports its dependencies.
Our automated testing takes place on both (x86) MacOS and Ubuntu Linux.

(installation/source)=

## Installing from source

The OpenFF Toolkit has a lot of dependencies, so we strongly encourage installation with a package manager. The [developer's guide](install_dev) describes setting up a development environment. If you're sure you want to install from source, check the [`conda-forge` recipe](https://github.com/conda-forge/openff-toolkit-feedstock/blob/main/recipe/meta.yaml) for current dependencies, install them, download and extract the source distribution from [GitHub](https://github.com/openforcefield/openff-toolkit/releases), and then run `pip`:

```shell-session
$ cd openff-toolkit
$ python -m pip install .
```

## Optional dependencies (toolkits)

The OpenFF Toolkit outsources many common computational chemistry algorithms to other toolkits.
Only one such toolkit is needed to gain access to all of the OpenFF Toolkit's features.
If more than one is available, the Toolkit allows the user to specify their preference with the `toolkit_registry` argument to most functions and methods.

The `openff-toolkit` package installs everything needed to run the toolkit, including the optional dependencies RDKit and AmberTools.
To install only the hard dependencies and provide your own optional dependencies, install the `openff-toolkit-base` package.

The OpenFF Toolkit requires an external toolkit for most functions.
Though a "built-in" toolkit is provided, it implements only a small number of functions and is intended primarily for testing.

There are certain differences in toolkit behavior between RDKit/AmberTools and OpenEye when reading a small fraction of molecules, and we encourage you to report any unexpected behavior that may be caused by toolkit differences to our [issue tracker](https://github.com/openforcefield/openff-toolkit/issues).

### RDKit

RDKit is a free and open source chemistry toolkit installed by default with the `openff-toolkit` package.
It provides most of the functionality that the OpenFF Toolkit relies on.

### AmberTools

AmberTools is a collection of free tools provided with the Amber MD software and installed by default with the `openff-toolkit` package.
It provides a free implementation of functionality required by OpenFF Toolkit and not provided by RDKit.

(installation/openeye)=

### OpenEye

The OpenFF Toolkit can optionally make use of the [OpenEye toolkit](https://www.eyesopen.com/toolkit-development) if the user has a license key installed.
Academic laboratories intending to release results into the public domain can [obtain a free license key](https://www.eyesopen.com/licensing-philosophy), while other users (including academics intending to use the software for purposes of generating protected intellectual property) must [pay to obtain a license](https://www.eyesopen.com/pricing).

To install the OpenEye toolkits:

```shell-session
$ mamba install -c openeye openeye-toolkits
```

Though OpenEye can be installed for free, using it requires a license file.
No essential `openff-toolkit` release capabilities *require* the OpenEye toolkit, but the Open Force Field developers make use of it in parameterizing new open source force fields.

### Check installed toolkits

All available toolkits are automatically registered in the `GLOBAL_TOOLKIT_REGISTRY`. The available toolkits and their versions can be inspected through the `registered_toolkit_versions` dictionary attribute:

```python
from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY
print(GLOBAL_TOOLKIT_REGISTRY.registered_toolkit_versions)
# {'The RDKit': '2022.03.5', 'AmberTools': '22.0', 'Built-in Toolkit': None}
```
