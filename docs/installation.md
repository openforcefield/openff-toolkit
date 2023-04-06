(installation)=

# Installation

(installation/conda)=

## Installing via `conda`

The simplest way to install the Open Force Field Toolkit is via the [conda](https://docs.conda.io/en/latest/) package manager.
We publish [packages](https://github.com/conda-forge/openff-toolkit-feedstock) via [`conda-forge`](https://conda-forge.org/). 
With Conda installed, use it to install the OpenFF Toolkit into a new environment:

```shell-session
$ conda create -n openff-toolkit -c conda-forge openff-toolkit
```

If you have Mamba installed, it is often faster:

```shell-session
$ mamba create -n openff-toolkit -c conda-forge openff-toolkit
```

If you do not have Mamba or Conda installed, see the [ecosystem installation documentation]. To use the new environment in a shell session, activate it:

```shell-session
$ conda activate openff-toolkit
```

:::{note}
Installation via the Conda package manager is the preferred method since all dependencies are automatically fetched and installed for you.
:::

[ecosystem installation documentation]: openff.docs:installation

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

## Single-file installer

As of release 0.4.1, single-file installers are available for each Open Force Field Toolkit release.
The single-file installer packages an entire Conda distribution complete will all dependencies needed to run the Toolkit.
These are provided primarily for users who do not have access to the Anaconda cloud for installing packages.
These installers have few requirements beyond a Linux or MacOS operating system and will, in one command, produce a functional Python executable containing the Open Force Field Toolkit, as well as all required dependencies.
The installers are very similar to the widely-used Miniconda `*.sh` files.
Accordingly, installation using the "single-file installer" does not require root access.

The installers are between 200 and 300 MB each, and can be downloaded from the "Assets" section of the Toolkit's [GitHub Releases page](https://github.com/openforcefield/openff-toolkit/releases/).
They are generated using a [workflow](https://github.com/openforcefield/toolkit-installer-constructor) leveraging the Conda Constructor utility.

Please report any installer difficulties to the [OFF Toolkit issue tracker](https://github.com/openforcefield/openff-toolkit/issues), as we hope to make this a major distribution channel for the toolkit moving forward.

### Installation

Download the appropriate installer (`openff-toolkit-<X.Y.Z>-<py3x>-<your platform>-x86_64.sh.zip`) from the "Assets" section at the bottom of the desired [release on GitHub](https://github.com/openforcefield/openff-toolkit/releases/).
Then, install the toolkit with the following command:

```shell-session
$ bash openff-toolkit-<X.Y.Z>-py37-<your platform>-x86_64.sh
```

and follow the prompts.

:::{note}
You must have write access to the installation directory.
This is generally somewhere in the user's home directory.
When prompted, we recommend NOT initializing the single-file installer.
:::

:::{warning}
We recommend that you do not install this package as root.
Conda is intended to support on-the-fly creation of several independent environments, and managing a multi-user Conda installation is [complicated.](https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/admin-multi-user-install.html)
:::

### Usage

Any time you want to use the single-file installer Conda environment in a terminal, run

```shell-session
$ source <install_directory>/etc/profile.d/conda.sh
$ conda activate base
```

Once the single-file installer `base` environment is activated, your shell session will use the Python installation (and any other provided executables) from it. Note that the single-file installer is a rare case in which we recommend altering the `base` environment, as the entire Conda distribution is centered on it. For more information about Conda environments, see [](openff.docs:managing_environments).

## Optional dependencies (toolkits)

The OpenFF Toolkit outsources many common computational chemistry algorithms to other toolkits. 
Only one such toolkit is needed to gain access to all of the OpenFF Toolkit's features. 
If more than one is available, the Toolkit allows the user to specify their preference with the `toolkit_registry` argument to most functions and methods.

The `openff-toolkit` package installs everything needed to run the toolkit, including the optional dependencies RDKit and AmberTools. 
To install only the hard dependencies and provide your own optional dependencies, install the `openff-toolkit-base` package.

The OpenFF Toolkit requires an external toolkit for most functions. 
Though a `builtin` toolkit is provided, it implements only a small number of functions and is intended primarily for testing.

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
$ conda install -c openeye -c conda-forge openeye-toolkits
```

Though OpenEye can be installed for free, using it requires a license file. 
No essential `openff-toolkit` release capabilities *require* the OpenEye toolkit, but the Open Force Field developers make use of it in parameterizing new open source force fields.


### Check installed toolkits

All available toolkits are automatically registered in the `GLOBAL_TOOLKIT_REGISTRY`. The available toolkits and their versions can be inspected through the `registered_toolkit_versions` dictionary:

```python
from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY
print(GLOBAL_TOOLKIT_REGISTRY.registered_toolkit_versions)
# {'The RDKit': '2022.03.5', 'AmberTools': '22.0', 'Built-in Toolkit': None}
```
