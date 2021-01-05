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

- The *core functionality* of the OFF Toolkit is to combine an Open Force Field ``ForceField`` and ``Topology`` to create an OpenMM ``System``.
- An OpenMM ``System`` contains *everything* needed to compute the potential energy of a system, except the coordinates.
- The OFF toolkit employs a modular "plugin" architecture wherever possible, providing a standard interface for contributed features.


Terminology
'''''''''''

Open Force Field Toolkit Concepts


``OFF Molecule``
  A graph representation of a molecule containing enough information to unambiguously parametrize it.
  Required data fields for an ``OFF Molecule`` are:

  - ``atoms``: element (integer), formal_charge (integer), is_aromatic (boolean), stereochemistry (R/S/None)
  - ``bonds``: order (integer), is_aromatic (boolean), stereochemistry (E/Z/None)

  There are several other optional attributes such as ``conformers`` and ``partial_charges`` that may be populated in the ``Molecule`` data structure.
  These are considered "optional" because they are not required for system creation, however if those fields are populated, the user MAY use them to override values that would otherwise be generated during system creation.

  A dictionary, ``Molecule.properties`` is exposed, which is a Python dict that can be populated with arbitrary data.
  This data should be considered cosmetic and should not affect system creation.
  Whenever possible, molecule serialization or format conversion should preserve this data.

``OFF System``
  An object that contains everything needed to calculate a molecular system's energy, except the atomic coordinates.
  Note that this does not exist yet, and that OpenMM ``System`` objects are being used for this purpose right now.

``OFF Topology``
  An object that efficiently holds many OFF ``Molecule`` objects.
  The atom indexing in a ``Topology`` may differ from those of the underlying ``Molecule``s

``OFF TopologyMolecule``
  The efficient data structures that make up an OFF Topology.
  There is one TopologyMolecule for each instance of a chemical species in a Topology.
  However, each unique chemical species has a single OFF Molecule representing it, which may be shared by multiple TopologyMolecules.
  TopologyMolecules contain an atom index map, as several copies of the same chemical species in a Topology may be present with different atom orderings.
  This data structure allows the OFF toolkit to only parametrize each unique Molecule once, and then write a copy of the assigned parameters out for each of the Molecule in the Topology (accounting for atom indexing differences in the process).


``OFF ForceField``
  An object generated from an OFFXML file (or other source of SMIRNOFF data).
  Most information from the SMIRNOFF data source is stored in this object's several ``ParameterHandler``s, however some top-level SMIRNOFF data is stored in the ``ForceField`` object itself.

``SMIRNOFF data``
  A hierarchical data structure that complies with the SMIRNOFF specification.
  This can be serialized in many formats, including XML (OFFXML).
  The subsections in a SMIRNOFF data source generally correspond to one energy term in the functional form of a force field.

``ParameterHandler``
  An object that has the ability to produce one component of an OpenMM ``System``, corresponding to one subsection in a SMIRNOFF data source.
  Most ``ParameterHandler`` objects contain a list of ``ParameterType`` objects.

``ParameterType``
  An object corresponding to a single SMARTS-based parameter.

``Cosmetic attribute``
  Data in a SMIRNOFF data source that does not correspond to a known attribute.
  These have no functional effect, but several programs use the extensibility of the OFFXML format to define additional attributes for their own use, and their workflows require the OFF toolkit to process the files while retaining these keywords.

Development Infrastructure

``CI``
    "Continuous integration" testing.

    Services that run frequently while the code is undergoing changes, ensuring that the codebase still installs and has the intended behavior.
    Currently, we use a service called `Travis CI <https://travis-ci.org>`_ for this.
    Every time we make commits to the ``master`` branch of the openff-toolkit Github repository, a set of virtual machines that mimic brand new Linux and Mac OSX computers are created, and follow build instructions specified in the repo's ``.travis.yml`` file to install the toolkit.
    After installing the OFF toolkit and its dependencies, these virtual machines run our test suite.
    If the tests all pass, the build "passes" (returns a green check mark on GitHub).
    If all the tests for a specific change to the ``master`` branch return green, then we know that the change has not broken the toolkit's existing functionality.
    When proposing code changes, we ask that contributors open a Pull Request (PR) on GitHub to merge their changes into the ``master`` branch.
    When a pull request is open, CI will run on the latest set of proposed changes and indicate whether they are safe to merge through status checks, summarized as a green check mark or red X.

``CodeCov``
  Code coverage.

  An extension to our testing framework that reports the fraction of our source code lines that were run during the tests.
  This functionality is actually the combination of several components -- Travis CI runs the tests using the ``pytest-cov`` package, and then uploads the results to the website codecov.io.
  This analysis is re-run with each change to the ``master`` branch, and a badge showing our coverage percentage is in the project README.

``LGTM``
  "Looks Good To Me".

  A service that analyzes the code in our repository for simple style and formatting issues.
  This service assigns a letter grade to codebases, and a badge showing our LGTM report is in the project README.

``RTD``
  ReadTheDocs.

  A service that compiles and renders the packages documentation (from the ``docs/`` folder).
  The documentation itself can be accessed from the ReadTheDocs badge in the README.

Modular design features
'''''''''''''''''''''''

There are a few areas where we've designed the toolkit with extensibility in mind.
Adding functionality at these interfaces should be considerably easier than in other parts of the toolkit, and we encourage experimentation and contribution on these fronts.

ParameterHandler
    A generic base class for objects that perform parametrization for one section in a SMIRNOFF data source.

    Each ParameterHandler-derived class MUST implement:
        - ``create_force(self, system, topology, **kwargs)``: takes an OpenMM ``System`` and a OpenFF ``Topology`` as input, as well as optional keyword arguments, and modifies the ``System`` to contain the appropriate parameters.
        - Class-level ``ParameterAttributes`` and ``IndexedParameterAttributes``: These correspond to the header-level attributes in a SMIRNOFF data source.
          For example,, the ``Bonds`` tag in the SMIRNOFF spec has an optional ``fractional_bondorder_method`` field, which corresponds to the line  ``fractional_bondorder_method = ParameterAttribute(default=None)`` in the ``BondHandler`` class definition.
          The ``ParameterAttribute`` and ``IndexedParameterAttribute`` classes offer considerable flexibility for validating inputs.
          Defining these attributes at the class level implements the corresponding behavior in the default ``__init__`` function.
        - Class-level definitions ``_MAX_SUPPORTED_SECTION_VERSION`` and ``_MAX_SUPPORTED_SECTION_VERSION``.
          ParameterHandler versions allow us to evolve ParameterHandler behavior in a controlled, recorded way.
          Force field development is experimental by nature, and it is unlikely that the initial choice of header attributes is suitable for all use cases.
          Recording the "versions" of a SMIRNOFF spec tag allows us to encode the default behavior and API of a specific generation of ParameterHandlers, while allowing the safe addition of new attributes and behaviors.
    - Each ParameterHandler-derived class MAY implement:
        - ``known_kwargs``: Keyword arguments passed to ``ForceField.create_openmm_system`` are validated against the ``known_kwargs`` lists of each ParameterHandler that the ForceField owns.
          If present, these kwargs and their values will be passed on to the ParameterHandler.
        - ``to_dict``: converts the ParameterHandler to a hierarchical dict compliant with the SMIRNOFF specification.
          The default implementation of this function should suffice for most developers.
        - ``check_handler_compatibility``: Checks whether this ParameterHandler is "compatible" with another.
          This function is used when a ForceField is attempted to be constructed from *multiple* SMIRNOFF data sources, and it is necessary to check that two sections with the same tagname can be combined in a sane way.
          For example, if the user instructed two ``vdW`` sections to be read, but the sections defined different vdW potentials, then this function should raise an Exception indicating that there is no safe way to combine the parameters.
          The default implementation of this function should suffice for most developers.
        - ``postprocess_system``: operates identically to ``create_force``, but is run after each ParameterHandlers' ``create_force`` has already been called.
          The default implementation of this method simply does nothing, and should suffice for most developers.

.. TODO : fill in the modular components below

ParameterType

   ToolkitRegistry

       ``ToolkitRegistry.from_object``  / ``ToolkitRegistry.from_smiles`` / ``OpenEyeToolkitWrapper.from_openeye`` / ``RDKitToolkitWrapper.from_rdkit``
        - These methods are a bit strange because they are effectively classmethods for ``FrozenMolecule`` and ``Molecule`` subclasses.
          In `PR #583 <https://github.com/openforcefield/openff-toolkit/pull/583>`_, jaimergp raised a concern that effectively boils down to "if I subclass ``Molecule`` into a new class, ``MyMol``, then I expect ``MyMol.from_rdkit`` to return an instance of ``MyMol``, not ``Molecule``.
          However, before this PR, methods like ``ToolkitRegistry.from_smiles`` didn't have any way to know what type of object they should return, and instead always returned ``Molecule`` objects.
          So as of  `PR #583 <https://github.com/openforcefield/openff-toolkit/pull/583>`_, ToolkitRegistry methods that produce a Molecule must take a private parameter, ``_cls``, indicating the type of object to return.
          This parameter should be of type ``type`` and should subclass ``FrozenMolecule``, or otherwise expose ``Molecule._add_atom``, ``._add_bond``, ``.add_conformer``, and ``.partial_charges``.


   Molecule.to_X

   Molecule.from_X

   Force field directories


.. TODO : fill in the sections below
.. Molecule definition
   '''''''''''''''''''
   
   Required stereochemistry
   ''''''''''''''''''''''''
   
   Conformation dependence
   '''''''''''''''''''''''
   
   
   
   Reliance on external dependencies
   '''''''''''''''''''''''''''''''''
   
   
   
   ForceField file paths
   '''''''''''''''''''''

.. TODO : expand this section
.. Documentation
   '''''''''''''
   If you define a new class, add new files to autodoc

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

1. Install the ``conda`` package manager as part of the Anaconda Distribution from `here <https://www.anaconda.com/distribution/>`_

2. Set up conda environment

.. code-block:: shell

    $ # Create a conda environment with the Open Force Field toolkit and its dependencies
    $ conda create --name openff-dev -c conda-forge -c openeye openff-toolkit openeye-toolkits
    $ conda activate openff-dev
    $ # Remove (only) the toolkit and replace it with a local install
    $ conda remove --force openff-toolkit
    $ git clone https://github.com/openforcefield/openff-toolkit
    $ cd openff-toolkit
    $ pip install -e .

3. Obtain and store Open Eye license somewhere like ``~/.oe_license.txt``.
   Optionally store the path in environmental variable ``OE_LICENSE``, i.e. using a command like ``echo
   "export OE_LICENSE=/Users/yournamehere/.oe_license.txt" >> ~/.bashrc``


Development Process
"""""""""""""""""""

Development of new toolkit features generally proceeds in the following stages:

* Begin a discussion on the `GitHub issue tracker <http://github.com/openforcefield/openff-toolkit/issues>`_ to determine big-picture "what should this feature do?" and "does it fit in the scope of the OFF Toolkit?"
    * `"... typically, for existing water models, we want to assign library charges" <https://github.com/openforcefield/openff-toolkit/issues/25>`_
* Start identifying details of the implementation that will be clear from the outset
    * `"Create a new "special section" in the SMIRNOFF format (kind of analogous to the BondChargeCorrections section) which allows SMIRKS patterns to specify use of library charges for specific groups <https://github.com/openforcefield/openff-toolkit/issues/25#issue-225173968>`_
    * `"Following #86, here's how library charges might work: ..." <https://github.com/openforcefield/openff-toolkit/issues/25#issuecomment-354636391>`_
* Create a branch or fork for development
    * The OFF Toolkit has one unusual aspect of its CI build process, which is that certain functionality requires the OpenEye toolkits, so the builds must contain a valid OpenEye license file.
      An encrypted OpenEye license is present in the OFF Toolkit GitHub repository, as ``oe_license.txt.enc``.
      Only Travis has the decryption key for this file.
      However, this setup poses the risk that anyone who can run Travis builds could simply print the contents of the license after decryption, which would put us in violation of our academic contract with OpenEye.
      For this reason, the OpenEye-dependent tests will be skipped on forks.
    * Note that creating a fork will prevent the OpenEye license from being decrypted on Travis


Contributing
""""""""""""

We always welcome `GitHub pull requests <https://github.com/openforcefield/openff-toolkit/pulls>`_.
For bug fixes, major feature additions, or refactoring, please raise an issue on the `GitHub issue tracker <http://github.com/openforcefield/openff-toolkit/issues>`_ first to ensure the design will be amenable to current developer plans.

How can I become a developer?
"""""""""""""""""""""""""""""

If you would like to contribute, please post an issue on the `GitHub issue tracker <http://github.com/openforcefield/openff-toolkit/issues>`_ describing the contribution you would like to make to start a discussion.

Style guide
"""""""""""

Development for the ``openff-toolkit`` conforms to the recommendations given by the `Software Development Best Practices for Computational Chemistry <https://github.com/choderalab/software-development>`_ guide.

The naming conventions of classes, functions, and variables follows `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_, consistently with the best practices guide. The naming conventions used in this library not covered by PEP8 are:
- Use ``file_path``, ``file_name``, and ``file_stem`` to indicate ``path/to/stem.extension``, ``stem.extension``, and ``stem`` respectively, consistently with the variables in the standard ``pathlib`` library.
- Use ``n_x`` to abbreviate "number of X` (e.g. `n_atoms`, `n_molecules`).

We place a high priority on code cleanliness and readability, even if code could be written more compactly. For example, 15-character variable names are fine. Triply nested list comprehensions are not.

The ``openff-toolkit`` is in the process of adopting code formatting tools ("linters") to maintain consistent style and remove the burden of adhering to these standards by hand. Currently, two are employed:
1. `Black <https://black.readthedocs.io/>`_, the uncompromising code formatter, automatically formats code with a consistent style.
1. `isort <https://timothycrosley.github.io/isort/>`_, sorts imports

There is a step in CI that uses these tools to check for a consistent style. These checks will use the most recent versions of each linter. To ensure that changes follow these standards, you can install and run these tools locally:

.. code-block:: shell

    $ conda install black isort -c conda-forge
    $ black openff
    $ isort openff

Anything not covered above is currently up to personal preference, but may change as new linters are added.
