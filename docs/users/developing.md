(developing)=

# Developing for the toolkit

:::{contents}
  :local:
:::

This guide is written with the understanding that our contributors are NOT professional software developers, but are instead computational chemistry trainees and professionals.
With this in mind, we aim to use a minimum of bleeding-edge technology and alphabet soup, and we will define any potentially unfamiliar processes or technologies the first time they are mentioned.
We enforce use of certain practices (tests, formatting, coverage analysis, documentation) primarily because they are worthwhile upfront investments in the long-term sustainability of this project.
The resources allocated to this project will come and go, but we hope that following these practices will ensure that minimal developer time will maintain this software far into the future.

The process of contributing to the OpenFF Toolkit is more than just writing code.
Before contributing, it is a very good idea to start a discussion on the [Issue tracker](https://github.com/openforcefield/openff-toolkit/issues) about the functionality you'd like to add.
This Issue discussion will help us decide with you where in the codebase it should go, any overlapping efforts with other developers, and what the user experience should be.
Please note that the OpenFF Toolkit is intended to be used primarily as one piece of larger workflows, and that simplicity and reliability are two of our primary goals.
Often, the cost/benefit of new features must be discussed, as a complex codebase is harder to maintain.
When new functionality is added to the OpenFF Toolkit, it becomes our responsibility to maintain it, so it's important that we understand contributed code and are in a position to keep it up to date.

## Overview

### Philosophy

- The *core functionality* of the OpenFF Toolkit is to combine an Open Force Field [`ForceField`](openff.toolkit.typing.engines.smirnoff.forcefield.ForceField) and [`Topology`](openff.toolkit.topology.Topology) to create an OpenMM [`System`](http://docs.openmm.org/latest/api-python/generated/openmm.openmm.System.html).
- An OpenMM [`System`](http://docs.openmm.org/latest/api-python/generated/openmm.openmm.System.html) contains *everything* needed to compute the potential energy of a system, except the coordinates and (optionally) box vectors.
- The OpenFF toolkit employs a modular "plugin" architecture wherever possible, providing a standard interface for contributed features.


### Terminology

For high-level toolkit concepts and terminology important for both development and use of the Toolkit, see the [core concepts page](concepts).

#### SMIRNOFF and the OpenFF Toolkit

SMIRNOFF data
: A hierarchical data structure that complies with the [SMIRNOFF specification](smirnoff).
  This can be serialized in many formats, including XML (OFFXML).
  The subsections in a SMIRNOFF data source generally correspond to one energy term in the functional form of a force field.

Cosmetic attribute
: Data in a SMIRNOFF data source that does not correspond to a known attribute.
  These have no functional effect, but several programs use the extensibility of the OFFXML format to define additional attributes for their own use, and their workflows require the OFF toolkit to process the files while retaining these keywords.

#### Development Infrastructure

Continuous Integration (CI)
: Tests that run frequently while the code is undergoing changes, ensuring that the codebase still installs and has the intended behavior.
  Currently, we use a service called [GitHub Actions](https://github.com/features/actions) for this.
  CI jobs run every time a commit is made to the `main` branch of the `openff-toolkit` GitHub repository or in a PR opened against it.
  These runs start by booting virtual machines that mimic brand new Linux and macOS computers.
  They then follow build instructions (see the `.github/workflows/CI.yml` file) to install the toolkit.
  After installing the OpenFF Toolkit and its dependencies, these virtual machines run our test suite.
  If the tests all pass, the build "passes" (returns a green check mark on GitHub).

  If all the tests for a specific change to the `main` branch return green, then we know that the change has not broken the toolkit's existing functionality.
  When proposing code changes, we ask that contributors open a Pull Request (PR) on GitHub to merge their changes into the `main` branch.
  When a pull request is open, CI will run on the latest set of proposed changes and indicate whether they are safe to merge through status checks, summarized as a green check mark or red cross.

CodeCov
: An extension to our testing framework that reports the fraction of our source code lines that were run during the tests (our "code coverage").
  This functionality is actually the combination of several components -- GitHub Actions runners run the tests using `pytest` with the `pytest-cov` plugin, and then coverage reports are uploaded to [CodeCov's website](https://codecov.io).
  This analysis is re-run alongside the rest of our CI, and a badge showing our coverage percentage is in the project README.

"Looks Good To Me" (LGTM)
: A service that analyzes the code in our repository for simple style and formatting issues.
  This service assigns a letter grade to the codebase, and a badge showing our LGTM report is in the project README.

ReadTheDocs (RTD)
: A service that compiles and renders the package's documentation (from the `docs/` folder).
  The documentation itself can be accessed from the ReadTheDocs badge in the README.
  It is compiled by RTD alongside the other CI checks, and the compiled documentation for a pull request can be viewed by clicking the "details" link after the status.

### User Experience

One important aspect of how we make design decisions is by asking "who do we envision using this software, and what would they want it to do here?".
There is a wide range of possible users, from non-chemists, to students/trainees, to expert computational medicinal chemists.
We have decided to build functionality intended for use by *expert medicinal chemists*, and whenever possible, add fatal errors if the toolkit risks doing the wrong thing.
So, for example, if a molecule is loaded with an odd ionization state, we assume that the user has input it this way intentionally.

This design philosophy inevitably has trade-offs --- For example, the OpenFF Toolkit will give the user a hard time if they try to load a "dirty" molecule dataset, where some molecules have errors or are not described in enough detail for the toolkit to unambiguously parametrize them.
If there is risk of misinterpreting the molecule (for example, bond orders being undefined or chiral centers without defined stereochemistry), the toolkit should raise an error that the user can override.
In this regard we differ from RDKit, which is more permissive in the level of detail it requires when creating molecules.
This makes sense for RDKit's use cases, as several of its analyses can operate with a lower level of detail about the molecules.
Often, the same design decision is the best for all types of users, but when we do need to make trade-offs, we assume the user is an expert.

At the same time, we aim for "automagic" behavior whenever a decision will clearly go one way over another.
System parameterization is an inherently complex topic, and the OFF toolkit would be nearly unusable if we required the user to explicitly approve every aspect of the process.
For example, if a `Topology` has its `box_vectors` attribute defined, we assume that the resulting OpenMM `System` should be periodic.



## Modular design features

There are a few areas where we've designed the toolkit with extensibility in mind.
Adding functionality at these interfaces should be considerably easier than in other parts of the toolkit, and we encourage experimentation and contribution on these fronts.

These features have occasionally confusing names.
"Parameter" here refers to a single value in a force field, as it is generally used in biophysics; it does not refer to an argument to a function.
"Attribute" is used to refer to an XML attribute, which allows data to be defined for a particular tag; it does not refer to a member of a Python class or object.
For example, in the following XML excerpt the `<SMIRNOFF>` tag has the attributes `version` and `aromaticity_model`:

```xml
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
  ...
</SMIRNOFF>
```

"Member" is used here to describe Python attributes.
This terminology is borrowed for the sake of clarity in this section from languages like C++ and Java.


### `ParameterAttribute`
A [`ParameterAttribute`](openff.toolkit.typing.engines.smirnoff.parameters.ParameterAttribute) is a single value that can be validated at runtime.

A `ParameterAttribute` can be instantiated as Python class or instance members to define the kinds of value that a particular parameter can take.
They are used in the definitions of both `ParameterHandler` and `ParameterType`.
The sorts of values a `ParameterAttribute` can take on are restricted by runtime validation.
This validation is highly customizable, and may do things like allowing only certain values for a string or enforcing the correct units or array dimensions on the value; in fact, the validation can be defined using arbitrary code.
The name of a `ParameterAttribute` should correspond exactly to the corresponding attribute in an OFFXML file.

#### `IndexedParameterAttribute`
An [`IndexedParameterAttribute`](openff.toolkit.typing.engines.smirnoff.parameters.IndexedParameterAttribute) is a  `ParameterAttribute` with a sequence of values, rather than just one. Each value in the sequence is indexed by an integer.

The exact name of an `IndexedParameterAttribute` is NOT expected to appear verbatim in a OFFXML file, but instead should appear with a numerical integer suffix.
For example the `IndexedParameterAttribute` `k` should only appear as `k1`, `k2`, `k3`, and so on in an OFFXML.
The current implementation requires this indexing to start at 1 and subsequent values be contiguous (no skipping numbers), but does not enforce an upper limit on the integer.

For example, dihedral torsions are often parameterized as the sum of several sine wave terms.
Each of the parameters of the sine wave `k`, `periodicity`, and `phase` is implemented as an `IndexedParameterAttribute`.

#### `MappedParameterAttribute`
A [`MappedParameterAttribute`](openff.toolkit.typing.engines.smirnoff.parameters.MappedParameterAttribute) is a  `ParameterAttribute` with several values, with some arbitrary mapping to access values.

#### `IndexedMappedParameterAttribute`
An [`IndexedMappedParameterAttribute`](openff.toolkit.typing.engines.smirnoff.parameters.IndexedMappedParameterAttribute) is a `ParameterAttribute` with a sequence of maps of values.

### `ParameterHandler`
[`ParameterHandler`](openff.toolkit.typing.engines.smirnoff.parameters.ParameterHandler) is a generic base class for objects that perform parameterization for one section in a SMIRNOFF data source.
A `ParameterHandler` has the ability to produce one component of an OpenMM `System`.
Extend this class to add a support for a new force or energy term to the toolkit.

Each `ParameterHandler`-derived class MUST implement the following methods and define the following attributes:

- Class members `ParameterAttributes`: These correspond to the header-level attributes in a SMIRNOFF data source.
For example, the `Bonds` tag in the SMIRNOFF spec has an optional `fractional_bondorder_method` field, which corresponds to the line  `fractional_bondorder_method = ParameterAttribute(default=None)` in the `BondHandler` class definition.
The `ParameterAttribute` and `IndexedParameterAttribute` classes offer considerable flexibility for validating inputs.
Defining these attributes at the class level implements the corresponding behavior in the default `__init__` function.
- Class members `_MIN_SUPPORTED_SECTION_VERSION` and `_MAX_SUPPORTED_SECTION_VERSION`.
`ParameterHandler` versions allow us to evolve `ParameterHandler` behavior in a controlled, recorded way.
Force field development is experimental by nature, and it is unlikely that the initial choice of header attributes is suitable for all use cases.
Recording the "versions" of a SMIRNOFF spec tag allows us to encode the default behavior and API of a specific generation of a `ParameterHandler`, while allowing the safe addition of new attributes and behaviors.
If these attributes are not defined, defaults in the base class will apply and updates introducing new versions may break the existing code.

Each `ParameterHandler`-derived class MAY implement:
  - `_KWARGS`: Keyword arguments passed to `ForceField.create_openmm_system` are validated against the `_KWARGS` lists of each ParameterHandler that the ForceField owns.
    If present, these keyword arguments and their values will be passed on to the `ParameterHandler`.
  - `_TAGNAME`: The name of the SMIRNOFF OFFXML tag used to parameterize the class.
    This tag should appear in the top level within the [`<SMIRNOFF>`](https://openforcefield.github.io/standards/standards/smirnoff/#the-enclosing-smirnoff-tag) tag; see the [Parameter generators](https://openforcefield.github.io/standards/standards/smirnoff/#parameter-generators) section of the SMIRNOFF specification.
  - `_INFOTYPE`: The `ParameterType` subclass used to parse the elements in the `ParameterHandler`'s parameter list.
  - `_DEPENDENCIES`: A list of `ParameterHandler` subclasses that, when present, must run before this one.
    Note that this is *not* a list of `ParameterHandler` subclasses that are required by this one.
    Ideally, child classes of `ParameterHandler` are entirely independent, and energy components of a force field form distinct terms; when this is impossible, `_DEPENDENCIES` may be used to guarantee execution order.
  - `to_dict`: converts the `ParameterHandler` to a hierarchical dict compliant with the SMIRNOFF specification.
    The default implementation of this function should suffice for most developers.
  - `check_handler_compatibility`: Checks whether this `ParameterHandler` is "compatible" with another.
    This function is used when a `ForceField` is attempted to be constructed from *multiple* SMIRNOFF data sources, and it is necessary to check that two sections with the same tag name can be combined in a sane way.
    For example, if the user instructed two `vdW` sections to be read, but the sections defined different vdW potentials, then this function should raise an Exception indicating that there is no safe way to combine the parameters.
    The default implementation of this function should suffice for most developers.

### `ParameterType`
[`ParameterType`](openff.toolkit.typing.engines.smirnoff.parameters.ParameterType) is a base class for the SMIRKS-based parameters of a `ParameterHandler`.
Extend this alongside `ParameterHandler` to define and validate the force field parameters of a new force.
This is analogous to ParmEd's XType classes, like [BondType](https://parmed.github.io/ParmEd/html/api/parmed/parmed.html?highlight=bondtype#parmed.BondType). A `ParameterType` should correspond to a single SMARTS-based parameter.

For example, the Lennard-Jones potential can be parameterized through either the size `ParameterAttribute` `sigma` or `r_min`, alongside the energy `ParameterAttribute` `epsilon`. Both options are handled through the [`vdWType`](openff.toolkit.typing.engines.smirnoff.parameters.vdWHandler.vdWType) class, a subclass of `ParameterType`.


% TODO : fill in the modular components below
%    Molecule.to_X
%    Molecule.from_X
%    Force field directories


%  TODO : fill in the sections below
%  Molecule definition
%  '''''''''''''''''''
%
%  Required stereochemistry
%  ''''''''''''''''''''''''
%
%  Conformation dependence
%  '''''''''''''''''''''''
%
%
%
%  Reliance on external dependencies
%  '''''''''''''''''''''''''''''''''
%
%
%
%  ForceField file paths
%  '''''''''''''''''''''

%  TODO : expand this section
%  Documentation
%  '''''''''''''
%  If you define a new class, add new files to autodoc

### Non-bonded methods as implemented in OpenMM

The SMIRNOFF specification describes the contents of a force field, which can be implemented in a number of different ways in different molecular simulation engines.
The OpenMM implementation provided by the OpenFF Toolkit either produces an `openmm.System` containing a `openmm.NonbondedForce` object or raises an exception depending on how the non-bonded parameters are specified.
Exceptions are raised when parameters are incompatible with OpenMM (`IncompatibleParameterError`) or otherwise against spec (`SMIRNOFFSpecError`), and also when they are appropriate for the spec but not yet implemented in the toolkit (`SMIRNOFFSpecUnimplementedError`).
This table describes which [`NonbondedMethod`](http://docs.openmm.org/latest/userguide/application.html?highlight=ljpme#nonbonded-interactions) is used in the produced `NonbondedForce`, or else which exception is raised.

| vdw_method | electrostatics_method | periodic | OpenMM Nonbonded method or exception    | Common case |
|------------|-----------------------|----------|-----------------------------------------|-------------|
| cutoff     | Coulomb               | True     | raises `IncompatibleParameterError`     |             |
| cutoff     | Coulomb               | False    | `openmm.NonbondedForce.NoCutoff`        |             |
| cutoff     | reaction-field        | True     | raises `SMIRNOFFSpecUnimplementedError` |             |
| cutoff     | reaction-field        | False    | raises `SMIRNOFFSpecError`              |             |
| cutoff     | PME                   | True     | `openmm.NonbondedForce.PME`             | *           |
| cutoff     | PME                   | False    | `openmm.NonbondedForce.NoCutoff`        |             | 
| LJPME      | Coulomb               | True     | raises `IncompatibleParameterError`     |             |
| LJPME      | Coulomb               | False    | `openmm.NonbondedForce.NoCutoff`        |             | 
| LJPME      | reaction-field        | True     | raises `IncompatibleParameterError`     |             |
| LJPME      | reaction-field        | False    | raises `SMIRNOFFSpecError`              |             |
| LJPME      | PME                   | True     | `openmm.NonbondedForce.LJPME`           |             |
| LJPME      | PME                   | False    | `openmm.NonbondedForce.NoCutoff`        |             | 

Notes:
* The most commonly-used case (including the Parsley line) is in the fifth row (cutoff vdW, PME electrostatics, periodic topology) and marked with an asterisk.
* For all cases included a non-periodic topology, `openmm.NonbondedForce.NoCutoff` is currently used.
* Electrostatics method `reaction-field` can only apply to periodic systems, however it is not currently implemented.
* LJPME (particle mesh ewald for LJ/vdW interactions) is not yet fully described in the SMIRNOFF specification.
* In the future, the OpenFF Toolkit may create multiple `CustomNonbondedForce` objects in order to better de-couple vdW and electrostatic interactions.

## Contributing

We always welcome [GitHub pull requests](https://github.com/openforcefield/openff-toolkit/pulls).
For bug fixes, major feature additions, or refactoring, please raise an issue on the [GitHub issue tracker](http://github.com/openforcefield/openff-toolkit/issues) first to ensure the design will be amenable to current developer plans.
Development of new toolkit features generally proceeds in the following stages:

* Begin a discussion on the [GitHub issue tracker](http://github.com/openforcefield/openff-toolkit/issues) to determine big-picture "what should this feature do?" and "does it fit in the scope of the OpenFF Toolkit?"
    * ["... typically, for existing water models, we want to assign library charges"](https://github.com/openforcefield/openff-toolkit/issues/25)
* Start identifying details of the implementation that will be clear from the outset
    * ["Create a new "special section" in the SMIRNOFF format (kind of analogous to the BondChargeCorrections section) which allows SMIRKS patterns to specify use of library charges for specific groups](https://github.com/openforcefield/openff-toolkit/issues/25#issue-225173968)
    * ["Following #86, here's how library charges might work: ..."](https://github.com/openforcefield/openff-toolkit/issues/25#issuecomment-354636391)
* Create a branch or fork for development
    * The OpenFF Toolkit has one unusual aspect of its CI build process, which is that certain functionality requires the OpenEye toolkits, so the builds must contain a valid OpenEye license file.
    An OpenEye license is stored as an encrypted token within the `openforcefield` organization on GitHub.
    For security reasons, builds run from forks cannot access this key.
    Therefore, tests that depend on the OpenEye Toolkits will be skipped on forks.
    Contributions run on forks are still welcome, especially as features that do not interact directly with the OpenEye Toolktis are not likely affected by this limitation.


(install_dev)=
### Setting up a development environment

1. Install the `conda` package manager as part of the Anaconda Distribution from [here](https://www.anaconda.com/distribution/)

2. Set up conda environment:

    ```shell
    git clone https://github.com/openforcefield/openff-toolkit
    cd openff-toolkit/
    # Create a conda environment with dependencies from env/YAML file
    conda env create -n openff-dev -f devtools/conda-envs/test_env.yaml
    conda activate openff-dev
    # Perform editable/dirty dev install
    pip install -e .
    ```

3. Obtain and store Open Eye license somewhere like `~/.oe_license.txt`.
   Optionally store the path in environmental variable `OE_LICENSE`, i.e. using a command like `echo
   "export OE_LICENSE=/Users/yournamehere/.oe_license.txt" >> ~/.bashrc`

#### Building the Docs

The documentation is composed of two parts, a hand-written user guide and an auto-generated API reference.
Both are compiled by Sphinx, and can be automatically served and regenerated on changes with `sphinx-autobuild`.
Documentation for released versions is available at [ReadTheDocs](https://open-forcefield-toolkit.readthedocs.io/en/latest/).
ReadTheDocs also builds the documentation for each Pull Request opened on GitHub and keeps the output for 90 days.

To add the documentation dependencies to your existing `openff-dev` Conda environment:

```shell
# Add the documentation requirements to your Conda environment
conda env update --name openff-dev --file docs/environment.yml
conda install --name openff-dev -c conda-forge sphinx-autobuild
```

To build the documentation from scratch:

```shell
# Build the documentation
# From the openff-toolkit root directory
conda activate openff-dev
cd docs
make html
# Documentation can be found in docs/_build/html/index.html
```

To watch the source directory for changes and automatically rebuild the documentation and refresh your browser:

```shell
# Host the docs on a local HTTP server and rebuild when a source file is changed
# Works best when the docs have already been built
# From the openff-toolkit root directory
conda activate openff-dev
sphinx-autobuild docs docs/_build/html --watch openff
# Then navigate your web browser to http://localhost:8000
```

### Style guide

Development for the `openff-toolkit` conforms to the recommendations given by the [Software Development Best Practices for Computational Chemistry](https://github.com/choderalab/software-development) guide.

The naming conventions of classes, functions, and variables follows [PEP8](https://www.python.org/dev/peps/pep-0008/), consistently with the best practices guide.
The naming conventions used in this library not covered by PEP8 are:
- Use `file_path`, `file_name`, and `file_stem` to indicate `path/to/stem.extension`, `stem.extension`, and `stem` respectively, consistent with the variables in the [`pathlib` module](https://docs.python.org/3/library/pathlib.html) of the standard library.
- Use `n_x` to abbreviate "number of $x$" (e.g. `n_atoms`, `n_molecules`).

We place a high priority on code cleanliness and readability, even if code could be written more compactly.
For example, 15-character variable names are fine. Triply nested list comprehensions are not.

The `openff-toolkit` has adopted code formatting tools ("linters") to maintain consistent style and remove the burden of adhering to these standards by hand.
Currently, two are employed:
1. [Black](https://black.readthedocs.io/), the uncompromising code formatter, automatically formats code with a consistent style.
1. [isort](https://timothycrosley.github.io/isort/), sorts imports

There is a step in CI that uses these tools to check for a consistent style (see the file [`.github/workflows/lint.yml`](https://github.com/openforcefield/openff-toolkit/blob/main/.github/workflows/lint.yml)).
These checks will use the most recent versions of each linter.
To ensure that changes follow these standards, you can install and run these tools locally:

```shell
conda install black isort -c conda-forge
black openff
isort openff
```

Anything not covered above is currently up to personal preference, but this may change as new linters are added.

### Pre-commit

The [`pre-commit`](https://pre-commit.com/) tool can _optionally_ be used to automate some or all of the above style checks. It automatically runs other programs ("hooks") when you run `git commit`. It aborts the commit with an exit code if any of the hooks fail, i.e. if `black` reformats code. This project uses `pre-commit ci`, a free service that enforces style on GitHub using the `pre-commit` framework in CI.

To use `pre-commit` locally, first install it:

```shell
conda conda install pre-commit -c conda-forge  # also available via pip
```

Then, install the pre-commit hooks (note that it installs the linters into an isolated virtual environment, not the current conda environment):

```
pre-commit install
```

Hooks will now run automatically before commits. Once installed, they should run in a total of a few seconds.

If `pre-commit` is not used by the developer and style issues are found, a `pre-commit.ci` bot may commit directly to a PR to make these fixes. This bot should only ever alter styl and never make functional changes to code.

Note that tests (too slow) and type-checking (weird reasons) are not run by `pre-commit`. You should still manually run tests before commiting code.

## Checking code coverage locally

Run something like

```shell
$ python -m pytest --cov=openff --cov-config=setup.cfg --cov-report html openff
$ open htmlcov/index.html
```

to see a code coverage report. This uses [Coverage.py](https://coverage.readthedocs.io/) and also requires a pytest plugin (`pytest-cov`), both of which should be installed if you are using the provided conda environments. CodeCov provides this as an automated service with their own web frontend using a similar file generated by our CI. However, it can be useful to run this locally to check coverage changes (i.e. "is this function I added sufficiently coverged by tests?") without waiting for CI to complete.

## Supported Python versions

The OpenFF Toolkit roughly follows [NEP 29](https://numpy.org/neps/nep-0029-deprecation_policy.html).
As of April 2023 this means Python 3.9-3.10 are officially supported (3.11 is missing some upstream dependencies).
We develop, test, and distribute on macOS and Linux-based operating systems.
We do not currently support Windows.
Some CI builds run using only RDKit as a backend, some run using only OpenEye Toolkits, and some run using both installed at once.

The CI matrix is currently as follows:

:::{eval-rst}
+-----------------------+------------+-----------+-------------+------------+-----------+-------------+
|                       | Linux                                | macOS                                |
+-----------------------+------------+-----------+-------------+------------+-----------+-------------+
|                       | RDKit      | OpenEye   | RDKit + OE  | RDKit      | OpenEye   | RDKit + OE  |
+=======================+============+===========+=============+============+===========+=============+
| Python 3.8 and older  | No support after April 2023                                                 |
+-----------------------+------------+-----------+-------------+------------+-----------+-------------+
| Python 3.9            | **Test**   | **Test**  | **Test**    | **Test**   | **Test**  | **Test**    |
+-----------------------+------------+-----------+-------------+------------+-----------+-------------+
| Python 3.10           | **Test**   | Skip      | Skip        | **Test**   | Skip      | Skip        |
+-----------------------+------------+-----------+-------------+------------+-----------+-------------+
| Python 3.11           | **Test**   | Skip      | Skip        | **Test**   | Skip      | Skip        |
+-----------------------+------------+-----------+-------------+------------+-----------+-------------+
| Python 3.12 and newer | Pending official releases and upstream support                              | 
+-----------------------+------------+-----------+-------------+------------+-----------+-------------+
:::
