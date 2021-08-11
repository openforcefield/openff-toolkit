.. openforcefield documentation master file, created by
   sphinx-quickstart on Sun Dec  3 23:12:54 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Open Force Field Toolkit
========================

A modern, extensible library for molecular mechanics force field science from the `Open Force Field Initiative <http://openforcefield.org>`_

Copy-pasting from the `moldoc README <https://github.com/lukasturcani/moldoc>`_:

.. moldoc::

    # The content of a moldoc directive is just a Python script
    # which needs to define a moldoc_display_molecule variable.

    import moldoc.molecule as molecule

    moldoc_display_molecule = molecule.Molecule(
        atoms=(
            # molecule.Atom(atomic_number, position)
            molecule.Atom(6, (-0.06, -0.17, 0.)),
            molecule.Atom(17, (-1.35, 1.04, -0.04)),
            molecule.Atom(35, (1.65, 0.73, -0.06)),
            molecule.Atom(1, (-0.15, -0.88, -0.87)),
            molecule.Atom(1, (-0.09, -0.72, 0.97)),
        ),
        bonds=(
            # molecule.Bond(atom1_id, atom2_id, order)
            molecule.Bond(0, 1, 1),
            molecule.Bond(0, 2, 1),
            molecule.Bond(0, 3, 1),
            molecule.Bond(0, 4, 1),
        ),
    )


Using ``Molecule.to_moldoc()``:  # Make this link to the method

.. moldoc::

    from openff.toolkit.topology import Molecule

    molecule = Molecule.from_smiles("CCO")
    molecule.generate_conformers(n_conformers=1)

    moldoc_display_molecule = molecule.to_moldoc()


Getting started
---------------

.. toctree::
   :maxdepth: 1

   installation
   examples
   releasehistory
   faq

Using the toolkit
-----------------

.. toctree::
   :maxdepth: 1

   users/concepts
   users/molecule_cookbook
   users/smirnoff
   users/virtualsites
   users/developing

API documentation
-----------------

.. toctree::
  :maxdepth: 1

  topology
  typing
  utils
