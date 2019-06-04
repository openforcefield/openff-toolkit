Release History
===============

Releases follow the ``major.minor.micro`` scheme recommended by `PEP440 <https://www.python.org/dev/peps/pep-0440/#final-releases>`_, where

* ``major`` increments denote a change that may break API compatibility with previous ``major`` releases
* ``minor`` increments add features but do not break API compatibility
* ``micro`` increments represent bugfix releases or improvements in documentation



0.4.0 - SMIRNOFF Spec Update and Performance Optimization
---------------------------------------------------------

This update introduces the 0.3 SMIRNOFF spec.
This spec upgrade is the result of discussions about how to handle the evolution of data and parameter types as further functional forms are added to the SMIRNOFF spec.

The toolkit now contains functionality to "upgrade" the spec of a data source.
These methods are called automatically when loading an 0.1 or 0.2-spec data source.
This functionality allows the toolkit to continue to read 0.2-spec files, and also implements backwards-compatibility for 0.1-spec files.


WARNING: The 0.1 SMIRNOFF spec did not contain fields for several energy-determining parameters that are exposed in the 0.2 spec.
Thus, the 0.1-to-0.2 spec conversion process makes assumptions about the values that should be added for the newly-required fields.
These assumptions are intended to reproduce the spec update process for the smirnoff99Frosst family of force fields.
If you use the spec conversion functionality for other force field families, please carefully review the warning messages printed during spec conversion to ensure they are providing your desired behavior.


This update also greatly improves performance when running ``create_openmm_system`` on large topologies.

SMIRNOFF Spec Changes
"""""""""""""""""""""
* All individual sections are now versioned.
  The top-level ``SMIRNOFF`` tag, containing information like ``aromaticity_model``, ``Author``, and ``Date``, still has a version (currently 0.3).
  But, to allow for independent development of individual parameter types, each section (`such as ``Bonds``, ``Angles``, etc) has its own version as well (currently all 0.3).
* All units are now stored in expressions with their corresponding values. For example, distances are now stored as ``1.526*angstrom``, instead of storing the unit separately.
* The ``potential`` field for ``ProperTorsions`` and ``ImproperTorsions`` tags is no longer ``charmm``, but is rather ``k*(1+cos(periodicity*theta-phase))``.
  It was pointed out to us that CHARMM-style torsions deviate from this formula when the periodicity of a torsion term is 0.

New features
""""""""""""
Adds the following new functions:

* `PR #311 <https://github.com/openforcefield/openforcefield/pull/311>`_:
    * ``utils/utils.py``: Adds ``convert_0_2_smirnoff_to_0_3``, which takes a 0.2-spec SMIRNOFF data dict, and upgrades it to 0.3.
      This function is called automatically when creating a ``ForceField`` from a 0.2 spec OFFXML file.
    * ``utils/utils.py``: Adds ``convert_0_1_smirnoff_to_0_2``, which takes a 0.1-spec SMIRNOFF data dict, and upgrades it to 0.2.
      This function is called automatically when creating a ``ForceField`` from a 0.1 spec OFFXML file.
    * ``typing/engines/smirnoff/parameters.py``: Adds ``ParameterHandler`` and ``ParameterType`` ``add_cosmetic_attribute`` and ``delete_cosmetic_attribute`` functions.
      Once created, these can be accessed and modified as attributes of the underlying object (eg. ``ParameterType.my_cosmetic_attrib = 'blue'``)
      These functions are experimental, and we are interested in feedback on how  cosmetic attribute handling could be improved. (`See Issue #338 <https://github.com/openforcefield/openforcefield/issues/338>`_)
* `PR #329 <https://github.com/openforcefield/openforcefield/pull/329>`_: Improved runtime for system creation in large topologies.



API-breaking changes
""""""""""""""""""""
* ``ParameterType`` and ``ParameterHandler`` constructors now expect the ``version`` kwarg (per the SMIRNOFF spec change above)
  This requirement can be skipped by providing the kwarg ``skip_version_check=True``
* ``ParameterType`` and ``ParameterHandler`` functions no longer handle ``X_unit`` attributes in SMIRNOFF data (per the SMIRNOFF spec change above).
* The scripts in ``utilities/convert_frosst`` are now deprecated.
  This functionality is importance for provenance and will be migrated to the ``openforcefield/smirnoff99Frosst`` repository in the coming weeks.

0.3.0 - API Improvements
------------------------

Several improvements and changes to public API.

New features
""""""""""""

* `PR #292 <https://github.com/openforcefield/openforcefield/pull/292>`_: Implement ``Topology.to_openmm`` and remove ``ToolkitRegistry.toolkit_is_available``
* `PR #322 <https://github.com/openforcefield/openforcefield/pull/322>`_: Install directories for the lookup of OFFXML files through the entry point group ``openforcefield.smirnoff_forcefield_directory``. The ``ForceField`` class doesn't search in the ``data/forcefield/`` folder anymore (now renamed ``data/test_forcefields/``), but only in ``data/``.

API-breaking Changes
""""""""""""""""""""
* `PR #278 <https://github.com/openforcefield/openforcefield/pull/278>`_: Standardize variable/method names
* `PR #291 <https://github.com/openforcefield/openforcefield/pull/291>`_: Remove ``ForceField.load/to_smirnoff_data``, add ``ForceField.to_file/string`` and ``ParameterHandler.add_parameters``. Change behavior of ``ForceField.register_X_handler`` functions.

Bugfixes
"""""""" 
* `PR #327 <https://github.com/openforcefield/openforcefield/pull/327>`_: Fix units in tip3p.offxml (note that this file is still not loadable by current toolkit)
* `PR #325 <https://github.com/openforcefield/openforcefield/pull/325>`_: Fix solvent box for provided test system to resolve periodic clashes.
* `PR #325 <https://github.com/openforcefield/openforcefield/pull/325>`_: Add informative message containing Hill formula when a molecule can't be matched in ``Topology.from_openmm``.
* `PR #325 <https://github.com/openforcefield/openforcefield/pull/325>`_: Provide warning or error message as appropriate when a molecule is missing stereochemistry.
* `PR #316 <https://github.com/openforcefield/openforcefield/pull/316>`_: Fix formatting issues in GBSA section of SMIRNOFF spec
* `PR #308 <https://github.com/openforcefield/openforcefield/pull/308>`_: Cache molecule SMILES to improve system creation speed
* `PR #306 <https://github.com/openforcefield/openforcefield/pull/306>`_: Allow single-atom molecules with all zero coordinates to be converted to OE/RDK mols
* `PR #313 <https://github.com/openforcefield/openforcefield/pull/313>`_: Fix issue where constraints are applied twice to constrained bonds

0.2.2 - Bugfix release
----------------------

This release modifies an example to show how to parameterize a solvated system, cleans up backend code, and makes several improvements to the README.

Bugfixes
""""""""
* `PR #279 <https://github.com/openforcefield/openforcefield/pull/279>`_: Cleanup of unused code/warnings in main package ``__init__``
* `PR #259 <https://github.com/openforcefield/openforcefield/pull/259>`_: Update T4 Lysozyme + toluene example to show how to set up solvated systems
* `PR #256 <https://github.com/openforcefield/openforcefield/pull/256>`_ and `PR #274 <https://github.com/openforcefield/openforcefield/pull/274>`_: Add functionality to ensure that links in READMEs resolve successfully


0.2.1 - Bugfix release
----------------------

This release features various documentation fixes, minor bugfixes, and code cleanup.

Bugfixes
""""""""
* `PR #267 <https://github.com/openforcefield/openforcefield/pull/267>`_: Add neglected ``<ToolkitAM1BCC>`` documentation to the SMIRNOFF 0.2 spec
* `PR #258 <https://github.com/openforcefield/openforcefield/pull/258>`_: General cleanup and removal of unused/inaccessible code.
* `PR #244 <https://github.com/openforcefield/openforcefield/pull/244>`_: Improvements and typo fixes for BRD4:inhibitor benchmark

0.2.0 - Initial RDKit support
-----------------------------

This version of the toolkit introduces many new features on the way to a 1.0.0 release.

New features
""""""""""""

* Major overhaul, resulting in the creation of the `SMIRNOFF 0.2 specification <https://open-forcefield-toolkit.readthedocs.io/en/master/smirnoff.html>`_ and its XML representation
* Updated API and infrastructure for reference SMIRNOFF :class:`ForceField` implementation
* Implementation of modular :class:`ParameterHandler` classes which process the topology to add all necessary forces to the system.
* Implementation of modular :class:`ParameterIOHandler` classes for reading/writing different serialized SMIRNOFF forcefield representations
* Introduction of :class:`Molecule` and :class:`Topology` classes for representing molecules and biomolecular systems
* New :class:`ToolkitWrapper` interface to RDKit, OpenEye, and AmberTools toolkits, managed by :class:`ToolkitRegistry`
* API improvements to more closely follow `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ guidelines
* Improved documentation and examples

0.1.0
-----

This is an early preview release of the toolkit that matches the functionality described in the preprint describing the SMIRNOFF v0.1 force field format: `[DOI] <https://doi.org/10.1101/286542>`_.

New features
""""""""""""

This release features additional documentation, code comments, and support for automated testing.

Bugfixes
""""""""

Treatment of improper torsions
''''''''''''''''''''''''''''''

A significant (though currently unused) problem in handling of improper torsions was corrected.
Previously, non-planar impropers did not behave correctly, as six-fold impropers have two potential chiralities.
To remedy this, SMIRNOFF impropers are now implemented as three-fold impropers with consistent chirality.
However, current force fields in the SMIRNOFF format had no non-planar impropers, so this change is mainly aimed at future work.
