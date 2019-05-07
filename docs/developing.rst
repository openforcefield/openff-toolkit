.. _developing:

Developing for the toolkit
**************************

Style guide
"""""""""""

Development for the ``openforcefield`` toolkit conforms to the recommendations given by the `Software Development Best Practices for Computational Chemistry <https://github.com/choderalab/software-development>`_ guide.

The naming conventions of classes, functions, and variables follows `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_, consistently with the best practices guide. The naming conventions used in this library not covered by PEP8 are:
- Use ``file_path``, ``file_name``, and ``file_stem`` to indicate ``path/to/stem.extension``, ``stem.extension``, and ``stem`` respectively, consistently with the variables in the standard ``pathlib`` library.
- Use ``n_x`` to abbreviate "number of X` (e.g. `n_atoms`, `n_molecules`).

Contributing
""""""""""""

We always welcome `GitHub pull requests <https://github.com/openforcefield/openforcefield/pulls>`_.
For bug fixes, major feature additions, or refactoring, please raise an issue on the `GitHub issue tracker <http://github.com/openforcefield/openforcefield/issues>`_ first to ensure the design will be amenable to current developer plans.

How can I become a developer?
"""""""""""""""""""""""""""""

If you would like to contribute, please post an issue on the `GitHub issue tracker <http://github.com/openforcefield/openforcefield/issues>`_ describing the contribution you would like to make to start a discussion.
