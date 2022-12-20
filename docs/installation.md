(installation)=

# Installation

(installation/conda)=

## Installing via `conda`

The simplest way to install the Open Force Field Toolkit is via the [conda](https://docs.conda.io/en/latest/) package manager.
We publish [packages](https://github.com/conda-forge/openff-toolkit-feedstock) via [`conda-forge`](https://conda-forge.org/).

If you are using the [Anaconda](https://www.anaconda.com/products/individual#Downloads) scientific Python distribution, you already have the `conda` package manager installed.
If not, the quickest way to get started is to install the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) distribution, a lightweight, minimal installation of Python and the Conda package manager.
See the [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) documentation for detailed installation instructions.
We recommend [Miniforge](https://github.com/conda-forge/miniforge#readme), a drop-in replacement for Miniconda that uses the community-run `conda-forge` channel by default.

Once Conda is installed, use it to install the OpenFF Toolkit:

```shell-session
$ conda install -c conda-forge openff-toolkit
```

:::{note}
Installation via the Conda package manager is the preferred method since all dependencies are automatically fetched and installed for you.
:::

(installation/platforms)=

### OS support

The OpenFF Toolkit is pure Python, and we expect it to work on any platform that supports its dependencies.
Our automated testing takes place on both (x86) MacOS and Ubuntu Linux.


(installation/windows)=

#### Windows

For Windows support, we recommend using the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) (WSL) to run a Linux system integrated into Windows.
We strongly suggest using WSL2, if your hardware supports it, for a smoother experience.
WSL2 requires virtualization support in hardware.
This is available on most modern CPUs, but may require activation in the BIOS.

Once WSL is configured, installing and using the Toolkit is done exactly as it would be for Linux.
Note that by default, Jupyter Notebook will not be able to open a browser window and will log an error on startup; just ignore the error and open the link it provides in your ordinary Windows web browser.

:::{note}
WSL2 [does support](https://docs.microsoft.com/en-us/windows/wsl/tutorials/gpu-compute) GPU compute, at least with nvidia cards, but setting it up [takes some work](https://developer.nvidia.com/cuda/wsl).
:::


(installation/m1)=

#### macOS with M1 chips

The Toolkit supports Apple Silicon (M1, M2, etc.), but only through [Rosetta] because (as of August 2022) some upstream dependencies are not yet built for the `osx-arm64` architecture. 
Conda can be configured to [use Rosetta] with the `CONDA_SUBDIR=osx-64` environment variable or the `subdir` Conda config variable. 
We recommend using this on a per-environment basis so that it persists across updates and new installs, but does not affect existing setups:

```shell-session
$ CONDA_SUBDIR=osx-64 conda create --name openff -c conda-forge openff-toolkit
$ conda activate
$ conda config --env --set subdir osx-64
```

Alternatively, make this setting the global default by updating the system Conda config:

```
$ conda config --system --set subdir osx-64
```

Note that this will affect how Conda behaves with other environments.

[Rosetta]: https://support.apple.com/en-au/HT211861
[use Rosetta]: https://conda-forge.org/docs/user/tipsandtricks.html#installing-apple-intel-packages-on-apple-silicon

(conda_envs)=

### Conda environments

Conda environments that mix packages from the default channels and `conda-forge` can become inconsistent; to prevent this mixing, we recommend using `conda-forge` for all packages. The easiest way to do this is to install Conda with [Miniforge](https://github.com/conda-forge/miniforge#readme).

If you already have a complex environment, or you wish to install a version of the Toolkit that is incompatible with other software you have installed, you can install the Toolkit into a new environment:

```shell-session
$ conda create -c conda-forge --name offtk openff-toolkit
```

An environment must be activated in any new shell session to use the software installed in it:

```shell-session
$ conda activate offtk
```

### Upgrading your installation

To update an earlier `conda` installation of `openff-toolkit` to the latest release version, you can use `conda update`:

```shell-session
$ conda update -c conda-forge openff-toolkit
```

Note that this may update other packages or install new packages if the most recent release of the Toolkit requires it.

(installation/source)=

## Installing from source

The OpenFF Toolkit has a lot of dependencies, so we strongly encourage installation with a package manager. The [developer's guide](install_dev) describes setting up a development environment. If you're sure you want to install from source, check the [`conda-forge` recipe](https://github.com/conda-forge/openff-toolkit-feedstock/blob/main/recipe/meta.yaml) for current dependencies, install them, download and extract the source distribution from [GitHub](https://github.com/openforcefield/openff-toolkit/releases), and then run `pip`:

```shell-session
$ cd openff-toolkit
$ python -m pip install .
```

## Single-file installer

As of release 0.4.1, single-file installers are available for each Open Force Field Toolkit release.
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
When prompted, we recommend NOT making modifications to your `bash_profile`.
:::

:::{warning}
We recommend that you do not install this package as root.
Conda is intended to support on-the-fly creation of several independent environments, and managing a multi-user Conda installation is [complicated.](https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/admin-multi-user-install.html)
:::

### Usage

Any time you want to use this Conda environment in a terminal, run

```shell-session
$ source <install_directory>/etc/profile.d/conda.sh
$ conda activate base
```

Once the `base` environment is activated, your system will default to use Python (and other executables) from the newly installed Conda environment. For more information about Conda environments, see {ref}`conda_envs`

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
