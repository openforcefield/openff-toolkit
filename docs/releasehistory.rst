Release History
===============

Releases follow the ``major.minor.micro`` scheme recommended by `PEP440 <https://www.python.org/dev/peps/pep-0440/#final-releases>`_, where

* ``major`` increments denote a change that may break API compatibility with previous ``major`` releases
* ``minor`` increments add features but do not break API compatibility
* ``micro`` increments represent bugfix releases or improvements in documentation

0.2.0
-----

This version of the toolkit introduces many new features on the way to a 1.0.0 release.

New features
""""""""""""

* Major overhaul of the `SMIRNOFF specification <https://open-forcefield-toolkit.readthedocs.io/en/master/smirnoff.html>`_ and its XML representation
* Updated API and infrastructure for reference SMIRNOFF :class:`ForceField` implementation
* Implementation of modular :class:`ParameterHandler` classes which process the topology to add all necessary forces to the system..
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
