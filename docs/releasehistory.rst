Release History
===============

Releases follow the ``major.minor.micro`` scheme recommended by `PEP440 <https://www.python.org/dev/peps/pep-0440/#final-releases>`_, where

* ``major`` increments denote a change that may break API compatibility with previous ``major`` releases
* ``minor`` increments add features but do not break API compatibility
* ``micro`` increments represent bugfix releases or improvements in documentation

1.0.0
-----

This is the first major release of the openforcefield toolkit intended for production use.

New features
""""""""""""

* Complete overhaul of SMIRNOFF format
* Updated API and infrastructure for :class:`ForceField` with modular reader/writer support for different representations
* Introduction of new :class:`Topology` format

0.1.0
-----

This is an early preview release of the openforcefield toolkit that matches the functionality described in the preprint describing the SMIRNOFF v0.1 force field format: `[DOI] <https://doi.org/10.1101/286542>`_.

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
