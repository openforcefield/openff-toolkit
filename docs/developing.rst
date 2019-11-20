.. _developing:

Developing for the toolkit
**************************

Overview
""""""""

Introduction
''''''''''''

This guide is written with the understanding that our contributors are NOT professional software developers, but are instead computational chemistry trainees and professionals.
With this in mind, we aim to use a minimum of bleeding-edge technology and alphabet soup, and we will define any potentially unfamiliar processes or technologies the first time they are mentioned.
We enforce use of certain practices (tests, formatting, coverage analysis, documentation) primarily because they are worthwhile upfront investments in the long-term sustainability of this project.
The resources allocated to this project will come and go, but we hope that following these practices will ensure that minimal developer time will maintain this software far into the future.

The process of contributing to the OFF toolkit is more than just writing code.
Before contributing, it is a very good idea to start a discussion on the Issue tracker about the functionality you'd like to add.
This Issue will help us identify where in the codebase it should go, any overlapping efforts with other developers, and what the user experience should be.
Please note that the OFF toolkit is intended to be used primarily as one piece of larger workflows, and that simplicity and reliability are two of our primary goals.
Often, the cost/benefit of new features must be discussed, as a complex codebase is harder to maintain.
When new functionality is added to the OFF Toolkit, it becomes our responsibility to maintain it, so it's important that we understand contributed code and are in a position to keep it up to date.

Philosophy
''''''''''

- The `core functionality` of the OFF Toolkit is to combine a ForceField and a Topology to create a System
- A System contains `everything` needed to compute the potential energy of a system, except the coordinates.
- The OFF toolkit employs a modular "plugin" architecture wherever possible, providing a standard interface for contributed features.


Terminology
'''''''''''

Open Force Field Toolkit Concepts

- ``OFF Molecule``: A graph representation of a molecule containing enough information to unambiguously parametrize it.
   Required data fields for an ``OFF Molecule`` are:
    - Atoms: element (integer), formal_charge (integer), is_aromatic (boolean), stereochemistry (R/S/None)
    - Bonds: order (integer), is_aromatic (boolean), stereochemistry (E/Z/None),
- ``OFF System``: An object that contains everything needed to calculate a molecular system's energy, except the atomic coordinates.
  Note that this does not exist yet, and that OpenMM System objects are being used for this purpose right now.
- ``OFF Topology``: An object that efficiently holds many OFF Molecules.
- ``OFF ForceField``: An object generated from an OFFXML file (or other source of SMIRNOFF data).
  Most information from the SMIRNOFF data source is stored in an OFF ForceField's several ParameterHandlers, however some top-level SMIRNOFF data is stored in the ForceField itself.
- ``SMIRNOFF data``: A hierarchical data structure that complies with the SMIRNOFF specification.
  This can be serialized in many formats, including XML (OFFXML).
  The subsections in a SMIRNOFF data source generally correspond to one energy term in the functional form of a force field.
- ``ParameterHandler``: An object that has the ability to produce one component of a System, corresponding to one subsection in a SMIRNOFF data source.
  Most ParameterHandlers contain a list of ``ParameterType`` objects.
- ``ParameterType``: An object corresponding to a single SMARTS-based parameter.

Development Infrastructure

- ``CI``: "Continuous integration" testing.
  The process of testing that the codebase still installs and has the intended behavior.
  Currently, we use a service called "Travis" for this.
  Every time we change the ``master`` branch of the openforcefield Github repository, a set of virtual machines that mimic brand new Linux and Mac OSX computers are created, and follow build instructions specified in the repo's ``.travis.yml`` file to install the toolkit.
  After installing the OFF toolkit and its dependencies, these virtual machines run our test suite.
  If the tests all pass, the build returns a green chekc mark.
  If all the tests for a specific change to the ``master`` branch return green, then we know that the change has not broken the toolkit's existing functionality.
  When proposing code changes, we ask that contributors open a Pull Request (PR) on GitHub to merge their changes into the ``master`` branch.
  When a pull request is open, Travis will test the proposed changes and indicate whether they are safe to merge.
- ``CodeCov``: Code coverage.
  An extension to our testing framework that reports the fraction of our source code lines that were run during the tests.
  This functionality is actually the combination of several components -- Travis CI runs the tests using the ``pytest-cov`` package, and then uploads the results to the website codecov.io.
  This analysis is re-run with each change to the ``master`` branch, and a badge showing our coverage percentage is in the project README.
- ``LGTM``: "Looks Good To Me".
  A service that analyzes the code in our repository for simple style and formatting issues.
  This service assigns a letter grade to codebases, and a badge showing our LGTM report is in the project README.
- ``RTD``: ReadTheDocs.
  A service that compiles and renders the packages documentation (from the ``docs/`` folder).
  The documentation itself can be accessed from the ReadTheDocs badge in the README.

Modular design features
'''''''''''''''''''''''

- ParameterHandler
- ParameterType
- ToolkitRegistry
- Molecule.to/from_object
-


Conformation dependence
'''''''''''''''''''''''



Reliance on external dependencies
'''''''''''''''''''''''''''''''''



ForceField file paths
'''''''''''''''''''''

User Experience
'''''''''''''''

One important aspect of how we make design decisions is by asking "who do we envision using this software, and what would they want it to do here?".
There is a wide range of possible users, from non-chemists, to students/trainees, to expert computational medicinal chemists.
We have decided to build functionality intended for use by `expert medicinal chemists`, and whenever possible, add fatal errors if the toolkit risks doing the wrong thing.
So, for example, if a molecule is loaded with an odd ionization state, we assume that the user has input it this way intentionally.
This design philosophy invariably has tradeoffs -- For example, the OFF Toolkit will give the user a hard time if they try to load a "dirty" molecule dataset, where some molecules have errors or are not described in enough detail for the toolkit to unambiguously parametrize them.
If there is risk of misinterpreting the molecule (for example, bond orders being undefined or chiral centers without defined stereochemistry), the toolkit should raise an error that the user can override.
In this regard we differ from RDKit, which is more permissive in the level of detail it requires when creating molecules.
This makes sense for RDKit's use cases, as several of its analyses can operate with a lower level of detail about the molecules.
Often, the same design decision is the best for all types of users, and there is no need for discussion.
But when we do need to make tradeoffs, "assume the user is an expert" is our guiding principle.

At the same time, we aim for "automagic" behavior whenever a decision will clearly go one way over another.
System parametrization is an inherently complex topic, and the OFF toolkit would be nearly unusable if we required the user to explicitly approve every aspect of the process.
For example, if a ``Topology`` has its ``box_vectors`` attribute defined, we assume that the resulting ``System`` should be periodic.



Setting up a development environment
""""""""""""""""""""""""""""""""""""



Steps
"""""

Development of new toolkit features generally proceeds in the following stages:

* Begin a discussion on the `GitHub issue tracker <http://github.com/openforcefield/openforcefield/issues>`_ to determine big-picture "what should this feature do?" and "does it fit in the scope of the OFF Toolkit?"
    * `"... typically, for existing water models, we want to assign library charges" <https://github.com/openforcefield/openforcefield/issues/25>`_
* Start identifying details of the implementation that will be clear from the outset
    * `"Create a new "special section" in the SMIRNOFF format (kind of analogous to the BondChargeCorrections section) which allows SMIRKS patterns to specify use of library charges for specific groups <https://github.com/openforcefield/openforcefield/issues/25#issue-225173968>`_
    * `"Following #86, here's how library charges might work: ..." <https://github.com/openforcefield/openforcefield/issues/25#issuecomment-354636391>`_
* Create a branch or fork for development
    * The OFF Toolkit has one unusual aspect of its CI build process, which is that certain functionality requires the OpenEye toolkits, so the builds must contain a valid OpenEye license file.
      An encrypted OpenEye license is present in the OFF Toolkit GitHub repository, as ``oe_license.txt.enc``.
      Only Travis has the decryption key for this file.
      However, this setup poses the risk that anyone who can run Travis builds could simply print the contents of the license after decryption, which would put us in violation of our academic contract with OpenEye.
      For this reason, the OpenEye-dependent tests will be skipped on forks.
    * Note that creating a fork will prevent the OpenEye license from being decrypted on Travis, so a few options are possible:
        * If you aren't working on OpenEye-dependent functionality


Developing in forks -- License stuff

Contributing
""""""""""""

We always welcome `GitHub pull requests <https://github.com/openforcefield/openforcefield/pulls>`_.
For bug fixes, major feature additions, or refactoring, please raise an issue on the `GitHub issue tracker <http://github.com/openforcefield/openforcefield/issues>`_ first to ensure the design will be amenable to current developer plans.

How can I become a developer?
"""""""""""""""""""""""""""""

If you would like to contribute, please post an issue on the `GitHub issue tracker <http://github.com/openforcefield/openforcefield/issues>`_ describing the contribution you would like to make to start a discussion.

Style guide
"""""""""""

Development for the ``openforcefield`` toolkit conforms to the recommendations given by the `Software Development Best Practices for Computational Chemistry <https://github.com/choderalab/software-development>`_ guide.

The naming conventions of classes, functions, and variables follows `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_, consistently with the best practices guide. The naming conventions used in this library not covered by PEP8 are:
- Use ``file_path``, ``file_name``, and ``file_stem`` to indicate ``path/to/stem.extension``, ``stem.extension``, and ``stem`` respectively, consistently with the variables in the standard ``pathlib`` library.
- Use ``n_x`` to abbreviate "number of X` (e.g. `n_atoms`, `n_molecules`).

We place a high priority on code cleanliness and readability.
So, 15-character variable names are fine.
Triply nested list comprehensions are not.


Anything not covered above is up to personal preference.
To remove the human friction from code formatting, we will likely adopt a standard formatter like `black <https://github.com/psf/black>`_ in the near future.

