(installation)=

# Installation

## Installing via `conda`

The simplest way to install the Open Force Field Toolkit is via the [conda](https://docs.conda.io/en/latest/) package manager.
Packages are provided on the [`conda-forge` channel](https://conda-forge.org/) for Linux, OS X, and Win platforms.
The [`openff-toolkit` Anaconda Cloud page](https://anaconda.org/conda-forge/openff-toolkit) has useful instructions and [download statistics](https://anaconda.org/conda-forge/openff-toolkit/files).

If you are using the [Anaconda](https://www.anaconda.com/products/individual#Downloads) scientific Python distribution, you already have the `conda`package manager installed.
If not, the quickest way to get started is to install the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) distribution, a lightweight minimal installation of Anaconda Python.

On Linux, you can install the Python 3 version into `$HOME/miniconda3` with (on `bash` systems):

```shell-session
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash ./Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
$ source ~/miniconda3/etc/profile.d/conda.sh
$ conda activate base
```

On MacOS, use the `osx` binary

```shell-session
$ curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O
$ bash ./Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda3
$ source ~/miniconda3/etc/profile.d/conda.sh
$ conda activate base
```

You may want to add the new `source ~/miniconda3/etc/profile.d/conda.sh` line to your `~/.bashrc` file to ensure Anaconda Python can enabled in subsequent terminal sessions.
`conda activate base` will need to be run in each subsequent terminal session to return to the environment where the toolkit will be installed.

Note that `openff-toolkit` will be installed into this local Python installation, so that you will not need to worry about disrupting existing Python installations.

:::{note}
Installation via the Conda package manager is the preferred method since all dependencies are automatically fetched and installed for you.
:::

### Required dependencies

The `openff-toolkit` makes use of the [Conda Forge ](https://conda-forge.org/) free and open source community package repository:

```shell-session
$ conda config --add channels conda-forge
$ conda update --all
```

This only needs to be done once.

:::{note}
If automation is required, provide the `--yes` argument to `conda update` and `conda install` comamnds.
More information on the `conda` command-line API can be found in the [Conda online documentation ](https://conda.io/docs/commands.html).
:::

### Release build

You can install the latest stable release build of `openff-toolkit` via the `conda` package manager with

```shell-session
$ conda config --add channels conda-forge
$ conda install openff-toolkit
```

This version is recommended for all users not actively developing new force field parameterization algorithms.

:::{note}
The Conda package manager will install dependencies from binary packages automatically, including difficult-to-install packages such as OpenMM, numpy, and scipy. This is really the easiest way to get started.
:::

### Upgrading your installation

To update an earlier `conda` installation of `openff-toolkit` to the latest release version, you can use `conda update`:

```shell-session
$ conda update openff-toolkit
```

### Optional dependencies

This toolkit can optionally make use of the [OpenEye toolkit ](https://www.eyesopen.com/toolkit-development) if the user has a license key installed.
Academic laboratories intending to release results into the public domain can [obtain a free license key ](https://www.eyesopen.com/licensing-philosophy), while other users (including academics intending to use the software for purposes of generating protected intellectual property) must [pay to obtain a license ](https://www.eyesopen.com/pricing).

To install the OpenEye toolkits (provided you have a valid license file):

```shell-session
$ conda install --yes -c openeye openeye-toolkits
```

No essential `openff-toolkit` release capabilities *require* the OpenEye toolkit, but the Open Force Field developers make use of it in parameterizing new open source force fields.
It is known that there are certain differences in toolkit behavior between RDKit and OpenEye when reading a small fraction of molecules, and we encourage you to report any unexpected behavior that may be caused by toolkit differences to our [issue tracker ](https://github.com/openforcefield/openff-toolkit/issues).

## Alternative method: Single-file installer

As of release 0.4.1, single-file installers are available for each Open Force Field Toolkit release.
These are provided primarily for users who do not have access to the Anaconda cloud for installing packages.
These installers have few requirements beyond a Linux or OSX operating system and will, in one command, produce a functional Python executable containing the Open Force Field Toolkit, as well as all required dependencies.
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
Conda is intended to support on-the-fly creation of several independent environments, and managing a multi-user Conda installation is [somewhat involved.](https://docs.conda.io/projects/conda/en/latest/user-guide/configuration/admin-multi-user-install.html)
:::

### Usage

Any time you want to use this Conda environment in a terminal, run

```shell-session
$ source <install_directory>/etc/profile.d/conda.sh
$ conda activate base
```

Once the `base` environment is activated, your system will default to use Python (and other executables) from the newly installed Conda environment.

### Installing optional OpenEye toolkits

We're waiting on permission to redistribute the OpenEye toolkits inside the single-file installer, so for now the installers only ship with the open-source backend (RDKit+AmberTools).
With this in mind, the Conda environment created by the installer *contains the Conda package manager itself*, which can be used to install the OpenEye toolkits if you have access to the Anaconda cloud.

```shell-session
$ conda install -c openeye openeye-toolkits
```

:::{note}
The OpenEye Toolkits Conda package still requires a valid OpenEye license file to run.
:::

