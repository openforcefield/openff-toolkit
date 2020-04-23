Release History
===============

Releases follow the ``major.minor.micro`` scheme recommended by `PEP440 <https://www.python.org/dev/peps/pep-0440/#final-releases>`_, where

* ``major`` increments denote a change that may break API compatibility with previous ``major`` releases
* ``minor`` increments add features but do not break API compatibility
* ``micro`` increments represent bugfix releases or improvements in documentation

0.7.0 - Current development
---------------------------

Behavior changed
""""""""""""""""
- `PR #508 <https://github.com/openforcefield/openforcefield/pull/508>`_:
  In order to provide the same results for the same chemical species, regardless of input
  conformation, fractional bond order calculation methods now default to ignore input conformers
  and generate a new conformer of the molecule before running semiempirical calculations.
  Users can override this behavior by specifying the keyword argument
  ``use_conformers=molecule.conformers``
- `PR #281 <https://github.com/openforcefield/openforcefield/pull/281>`_: Closes
  `Issue #250 <https://github.com/openforcefield/openforcefield/issues/250>`_
  by adding support for partial charge I/O in SDF. The partial charges are stored as a property in the
  SDF molecule block under the tag ``<atom.dprop.PartialCharge>``.
- `PR #281 <https://github.com/openforcefield/openforcefield/pull/281>`_: If an OFFMol's
  ``partial_charges`` attribute is set to ``None`` (the default value), calling ``to_openeye`` will
  now produce a OE molecule with partial charges set to ``nan``. This would previously produce an OE
  molecule with partial charges of 0.0, which was a loss of information, since it wouldn't be clear
  whether the original OFFMol's partial charges were REALLY all-zero as opposed to ``None``. OpenEye toolkit
  wrapper methods such as ``from_smiles`` and ``from_file`` now produce OFFMols with
  ``partial_charges = None`` when appropriate (previously these would produce OFFMols with
  all-zero charges, for the same reasoning as above).
- `PR #281 <https://github.com/openforcefield/openforcefield/pull/281>`_: ``Molecule.to_rdkit``
  now sets partial charges on the RDAtom's ``PartialCharges`` property (this was previously set
  on the ``partial_charges`` property). If the OFFMol's partial_charges attribute is None, this property
  will not be defined on the RDAtoms.
- `PR #281 <https://github.com/openforcefield/openforcefield/pull/281>`_:
  Enforce the behavior during SDF I/O that a SDF may contain multiple MOLECULES, but that the OFF Toolkit
  DOES NOT assume that it contains multiple CONFORMERS of the same molecule. This is an
  important distinction, since otherwise there is ambiguity around whether properties of one
  entry in a SDF are shared among several molecule blocks or not, or how to resolve conflicts if properties
  are defined differently for several "conformers" of chemically-identical species (More info
  `here <https://docs.eyesopen.com/toolkits/python/oechemtk/oemol.html#dude-where-s-my-sd-data>`_.
  If the user requests the OFF
  Toolkit to write a multi-conformer ``Molecule`` to SDF, only the first conformer will be written.
  For more fine-grained control of writing properties, conformers, and partial charges, consider
  using ``Molecule.to_rdkit`` or ``Molecule.to_openeye`` and using the functionality offered by
  those packages.
- `PR #281 <https://github.com/openforcefield/openforcefield/pull/281>`_: Due to different
  constraints placed on the data types allowed by external toolkits, we make our best effort to
  preserve ``offmol.properties`` when converting molecules to other packages, but users should be aware that
  no guarantee of data integrity is made. The only data format for keys and values in the property dict that
  we will try to support through a roundtrip to another toolkit's Molecule object is ``string``.
- `PR #574 <https://github.com/openforcefield/openforcefield/pull/574>`_: Removed check that all
  partial charges are zero after assignment by `quacpac` when AM1BCC used for charge assignment.
  This check fails erroneously for cases in which the partial charge assignments are correctly all zero,
  such as for ``N#N``. It is also an unnecessary check given that `quacpac` will reliably indicate when
  it has failed to assign charges.

API-breaking changes
""""""""""""""""""""
- `PR #558 <https://github.com/openforcefield/openforcefield/pull/558>`_: Removes
  ``TopologyMolecule.topology_particle_start_index``, since the :py:class:`Topology <openforcefield.topology.Topology>`
  particle indexing system now orders :py:class:`TopologyVirtualSites <openforcefield.topology.TopologyVirtualSite>`
  after all atoms.
  :py:meth:`TopologyMolecule.topology_atom_start_index <openforcefield.topology.TopologyMolecule.topology_atom_start_index>`
  and
  :py:meth:`TopologyMolecule.topology_virtual_site_start_index <openforcefield.topology.TopologyMolecule.topology_virtual_site_start_index>`
  are still available to access the appropriate values in the respective topology indexing systems.
- `PR #508 <https://github.com/openforcefield/openforcefield/pull/508>`_:
  ``OpenEyeToolkitWrapper.compute_wiberg_bond_orders`` is now
  :py:meth:`OpenEyeToolkitWrapper.assign_fractional_bond_orders <openforcefield.utils.toolkits.OpenEyeToolkitWrapper.assign_fractional_bond_orders>`.
  The ``charge_model`` keyword is now ``bond_order_model``. The allowed values of this keyword have
  changed from ``am1`` and ``pm3`` to ``am1-wiberg`` and ``pm3-wiberg``, respectively.
- `PR #508 <https://github.com/openforcefield/openforcefield/pull/508>`_:
  ``Molecule.compute_wiberg_bond_orders`` is now
  :py:meth:`Molecule.assign_fractional_bond_orders <openforcefield.topology.Molecule.assign_fractional_bond_orders>`.


New features
""""""""""""
- `PR #508 <https://github.com/openforcefield/openforcefield/pull/508>`_:
  Adds :py:meth:`AmberToolsToolkitWrapper.assign_fractional_bond_orders <openforcefield.utils.toolkits.AmberToolsToolkitWrapper.assign_wiberg_bond_orders>`.
- `PR #469 <https://github.com/openforcefield/openforcefield/pull/469>`_:
  The :py:class:`Molecule <openforcefield.topology.Molecule>` adds
  :py:meth:`Molecule.has_unique_atom_names <openforcefield.topology.Molecule.has_unique_atom_names>`
  and :py:meth:`Molecule.has_unique_atom_names <openforcefield.topology.Molecule.generate_unique_atom_names>`.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_:
  The :py:class:`Molecule <openforcefield.topology.Molecule>` adds
  :py:meth:`Molecule.are_isomorphic <openforcefield.topology.Molecule.are_isomorphic>`
  and :py:meth:`Molecule.is_isomorphic_with <openforcefield.topology.Molecule.is_isomorphic_with>`
  and :py:meth:`Molecule.hill_formula <openforcefield.topology.Molecule.hill_formula>`
  and :py:meth:`Molecule.to_hill_formula <openforcefield.topology.Molecule.to_hill_formula>`
  and :py:meth:`Molecule.to_qcschema <openforcefield.topology.Molecule.to_qcschema>`
  and :py:meth:`Molecule.from_qcschema <openforcefield.topology.Molecule.from_qcschema>`
  and :py:meth:`Molecule.from_mapped_smiles <openforcefield.topology.Molecule.from_mapped_smiles>`
  and :py:meth:`Molecule.from_pdb_and_smiles <openforcefield.topology.Molecule.from_pdb_and_smiles>`
  and :py:meth:`Molecule.canonical_order_atoms <openforcefield.topology.Molecule.canonical_order_atoms>`
  and :py:meth:`Molecule.remap <openforcefield.topology.Molecule.remap>`
- `PR #506 <https://github.com/openforcefield/openforcefield/pull/506>`_:
   The :py:class:`Molecule <openforcefield.topology.Molecule>` adds
  :py:meth:`Molecule.find_rotatable_bonds <openforcefield.topology.Molecule.find_rotatable_bonds>`
- `PR #521 <https://github.com/openforcefield/openforcefield/pull/521>`_:
  Adds :py:meth:`Molecule.to_inchi <openforcefield.topology.Molecule.to_inchi>`
  and :py:meth:`Molecule.to_inchikey <openforcefield.topology.Molecule.to_inchikey>`
  and :py:meth:`Molecule.from_inchi <openforcefield.topology.Molecule.from_inchi>`
      .. warning::
         InChI was not designed as an molecule interchange format and using it as one is not recommended. Many round trip
         tests will fail when using this format due to a loss of information. We have also added support for fixed
         hydrogen layer nonstandard InChI which can help in the case of tautomers, but overall creating molecules from InChI should be
         avoided.
- `PR #529 <https://github.com/openforcefield/openforcefield/pull/529>`_: Adds the ability to write out to XYZ files via
  :py:meth:`Molecule.to_file <openforcefield.topology.Molecule.to_file>` Both single frame and multiframe XYZ files are supported.
  Note reading from XYZ files will not be supported due to the lack of connectivity information.
- `PR #535 <https://github.com/openforcefield/openforcefield/pull/535>`_: Extends the the API for the
  :py:meth:`Molecule.to_smiles <openforcefield.topology.Molecule.to_smiles>` to allow for the creation of cmiles
  identifiers through combinations of isomeric, explicit hydrogen and mapped smiles, the default settings will return
  isomeric explicit hydrogen smiles as expected.
        .. warning::
           Atom maps can be supplied to the properties dictionary to modify which atoms have their map index included,
           if no map is supplied all atoms will be mapped in the order they appear in the
           :py:class:`Molecule <openforcefield.topology.Molecule>`.
- `PR #563 <https://github.com/openforcefield/openforcefield/pull/563>`_:
  Adds ``test_forcefields/ion_charges.offxml``, giving `LibraryCharges` for monatomic ions.
- `PR #543 <https://github.com/openforcefield/openforcefield/pull/543>`_:
  Adds 3 new methods to the :py:class:`Molecule <openforcefield.topology.Molecule>` which allow the enumeration of molecule
  states these are :py:meth:`Molecule.enumerate_tautomers <openforcefield.topology.Molecule.enumerate_tautomers>`,
  py:meth:`Molecule.enumerate_stereoisomers <openforcefield.topology.Molecule.enumerate_stereoisomers>`,
  py:meth:`Molecule.enumerate_protomers <openforcefield.topology.Molecule.enumerate_protomers>`
      .. warning::
         Enumerate protomoers is currently only available through the OpenEye toolkit.
- `PR #573 <https://github.com/openforcefield/openforcefield/pull/573>`_:
  Adds ``quacpac`` error output to ``quacpac`` failure in ``Molecule.compute_partial_charges_am1bcc``.
- `PR #560 <https://github.com/openforcefield/openforcefield/issues/560>`_: Added visualization method to the the Molecule class.


Behavior changed
""""""""""""""""
- `PR #558 <https://github.com/openforcefield/openforcefield/pull/558>`_: The
  :py:class`Topology <openforcefield.topology.Topology>`
  particle indexing system now orders :py:class`TopologyVirtualSites <openforcefield.topology.TopologyVirtualSite>`
  after all atoms.
- `PR #469 <https://github.com/openforcefield/openforcefield/pull/469>`_:
  When running :py:meth:`Topology.to_openmm <openforcefield.topology.Topology.to_openmm>`, unique atom names
  are generated if the provided atom names are not unique (overriding any existing atom names). This
  uniqueness extends only to atoms in the same molecule. To disable this behavior, set the kwarg
  ``ensure_unique_atom_names=False``.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_:
  The :py:meth:`Molecule.__eq__ <openforcefield.topology.Molecule.__eq__>` now uses the new
  :py:meth:`Molecule.are_isomorphic <openforcefield.topology.Molecule.are_isomorphic>` to perform the
  similarity checking.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_:
  The :py:meth:`Topology.from_openmm <openforcefield.topology.Topology.from_openmm>` and
  :py:meth:`Topology.add_molecule <openforcefield.topology.Topology.add_molecule>` now use the
  :py:meth:`Molecule.are_isomorphic <openforcefield.topology.Molecule.are_isomorphic>` to match
  molecules.
- `PR #508 <https://github.com/openforcefield/openforcefield/pull/508>`_:
  In order to provide the same results for the same chemical species, regardless of input
  conformation, fractional bond order calculation methods now default to ignore input conformers
  and generate a new conformer of the molecule before running semiempirical calculations.
  Users can override this behavior by specifying the keyword argument
  ``use_conformers=molecule.conformers``
- `PR #544 <https://github.com/openforcefield/openforcefield/pull/544>`_: Raises 
  ``NotImplementedError`` when calling 
  :py:meth:`ParameterHandler.get_parameter   <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler.get_parameter>`, 
  which is not yet implemented, but would previously silently return ``None``.
- `PR #551 <https://github.com/openforcefield/openforcefield/pull/551>`_: Implemented the
  :py:meth:`ParameterHandler.get_parameter   <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler.get_parameter>` function.

Tests added
"""""""""""
- `PR #558 <https://github.com/openforcefield/openforcefield/pull/558>`_: Adds tests ensuring
  that the new Topology particle indexing system are properly implemented, and that TopologyVirtualSites
  reference the correct TopologyAtoms.
- `PR #469 <https://github.com/openforcefield/openforcefield/pull/469>`_: Added round-trip SMILES test
  to add coverage for :py:meth:`Molecule.from_smiles <openforcefield.topology.Molecule.from_smiles>`.
- `PR #469 <https://github.com/openforcefield/openforcefield/pull/469>`_: Added tests for unique atom
  naming behavior in  :py:meth:`Topology.to_openmm <openforcefield.topology.Topology.to_openmm>`, as
  well as tests of the ``ensure_unique_atom_names=False`` kwarg disabling this behavior.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_: Added tests for
  :py:meth:`Molecule.hill_formula <openforcefield.topology.Molecule.hill_formula>` and
  :py:meth:`Molecule.to_hill_formula <openforcefield.topology.Molecule.to_hill_formula>` for the
  various supported input types.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_: Added round-trip test for
  :py:meth:`Molecule.from_qcschema <openforcefield.topology.Molecule.from_qcschema>` and
  :py:meth:`Molecule.to_qcschema <openforcefield.topology.Molecule.to_qcschema>`.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_: Added tests for
  :py:meth:`Molecule.is_isomorphic_with <openforcefield.topology.Molecule.is_isomorphic_with>` and
  :py:meth:`Molecule.are_isomorphic <openforcefield.topology.Molecule.are_isomorphic>`
  with various levels of isomorphic graph matching.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_: Added toolkit dependent tests
  for :py:meth:`Molecule.canonical_order_atoms <openforcefield.topology.Molecule.canonical_order_atoms>`
  due to differences in the algorithms used.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_: Added a test for
  :py:meth:`Molecule.from_mapped_smiles <openforcefield.topology.Molecule.from_mapped_smiles>` using
  the molecule from issue #412 to ensure it is now fixed.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_: Added a test for
  :py:meth:`Molecule.remap <openforcefield.topology.Molecule.remap>`, this also checks for expected
  error when the mapping is not complete.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_: Added tests for
  :py:meth:`Molecule.from_pdb_and_smiles <openforcefield.topology.Molecule.from_pdb_and_smiles>`
  to check for a correct combination of smiles and PDB and incorrect combinations.
- `PR #509 <https://github.com/openforcefield/openforcefield/pull/509>`_: Added test for
  :py:meth:`Molecule.chemical_environment_matches <openforcefield.topology.Molecule.chemical_environment_matches>`
  to check that the complete set of matches is returned.
- `PR #509 <https://github.com/openforcefield/openforcefield/pull/509>`_: Added test for
  :py:meth:`Forcefield.create_openmm_system <openforcefield.typing.engines.smirnoff.forcefield.create_openmm_system>`
  to check that a protein system can be created.
- `PR #506 <https://github.com/openforcefield/openforcefield/pull/506>`_: Added a test for the molecule
   identified in issue #513 as loosing aromaticity when converted to rdkit.
- `PR #506 <https://github.com/openforcefield/openforcefield/pull/506>`_: Added a verity of toolkit dependent tests
   for identifying rotatable bonds while ignoring the user requested types.
- `PR #521 <https://github.com/openforcefield/openforcefield/pull/521>`_: Added toolkit independent round-trip InChI
  tests which add coverage for :py:meth:`Molecule.to_inchi <openforcefield.topology.Molecule.to_inchi>` and
  :py:meth:`Molecule.from_inchi <openforcefield.topology.Molecule.from_inchi>`. Also added coverage for bad inputs and
  :py:meth:`Molecule.to_inchikey <openforcefield.topology.Molecule.to_inchikey>`.
- `PR #529 <https://github.com/openforcefield/openforcefield/pull/529>`_: Added to XYZ file coverage tests.
- `PR #563 <https://github.com/openforcefield/openforcefield/pull/563>`_: Added `LibraryCharges` parameterization test
  for monatomic ions in ``test_forcefields/ion_charges.offxml``.
- `PR #543 <https://github.com/openforcefield/openforcefield/pull/543>`_: Added tests to assure that state enumeration can
  correctly find molecules tautomers, stereoisomers and protomers when possible.
- `PR #573 <https://github.com/openforcefield/openforcefield/pull/573>`_: Added test for ``quacpac`` error output
  for ``quacpac`` failure in ``Molecule.compute_partial_charges_am1bcc``.
- `PR #579 <https://github.com/openforcefield/openforcefield/pull/579>`_: Adds regression tests to ensure RDKit can be
  be used to write multi-model PDB files.


Bugfixes
""""""""
- `PR #558 <https://github.com/openforcefield/openforcefield/pull/558>`_: Fixes a bug where
  :py:meth:`TopologyVirtualSite.atoms <openforcefield.topology.TopologyVirtualSite.atoms>` would
  not correctly apply ``TopologyMolecule`` atom ordering on top of the reference molecule ordering,
  in cases where the same molecule appears multiple times, but in a different order, in the same Topology.
- `Issue #460 <https://github.com/openforcefield/openforcefield/issues/460>`_: Creates unique atom
  names in :py:meth:`Topology.to_openmm <openforcefield.topology.Topology.to_openmm>` if the existing
  ones are not unique. The lack of unique atom names had been causing problems in workflows involving
  downstream tools that expect unique atom names.
- `Issue #448 <https://github.com/openforcefield/openforcefield/issues/448>`_: We can now make molecules
  from mapped smiles using :py:meth:`Molecule.from_mapped_smiles <openforcefield.topology.Molecule.from_mapped_smiles>`
  where the order will correspond to the indeing used in the smiles.
  Molecules can also be re-indexed at any time using the
  :py:meth:`Molecule.remap <openforcefield.topology.Molecule.remap>`.
- `Issue #462 <https://github.com/openforcefield/openforcefield/issues/462>`_: We can now instance the
  :py:class:`Molecule <openforcefield.topology.Molecule>` from a QCArchive entry record instance or dictionary
  representation.
- `Issue #412 <https://github.com/openforcefield/openforcefield/issues/412>`_: We can now instance the
  :py:class:`Molecule <openforcefield.topology.Molecule>` using
  :py:meth:`Molecule.from_mapped_smiles <openforcefield.topology.Molecule.from_mapped_smiles>`. This resolves
  an issue caused by RDKit considering atom map indices to be a distinguishing feature of an atom, which led
  to erroneous definition of chirality (as otherwise symmetric substituents would be seen as different).
  We anticipate that this will reduce the number of times you need to
  type ``allow_undefined_stereo=True`` when processing molecules that do not actually contain stereochemistrty.
- `Issue #513 <https://github.com/openforcefield/openforcefield/issues/513>`_: The
   :py:meth:`Molecule.to_rdkit <openforcefield.topology.Molecule.to_rdkit>` now re-sets the aromaticity model
   after sanitizing the molecule.
- `Issue #500 <https://github.com/openforcefield/openforcefield/issues/500>`_: The
   :py:meth:`Molecule.find_rotatable_bonds <openforcefield.topology.Molecule.find_rotatable_bonds>` has been added
   which returns a list of rotatable :py:class:`Bond <openforcefield.topology.Bond>` instances for the molecule.
- `Issue #491 <https://github.com/openforcefield/openforcefield/issues/491>`_: We can now parse large molecules without hitting a match limit cap.
- `Issue #474 <https://github.com/openforcefield/openforcefield/issues/474>`_: We can now  convert molecules to InChI and
   InChIKey and from InChI.
- `Issue #523 <https://github.com/openforcefield/openforcefield/issues/523>`_: The
   :py:meth:`Molecule.to_file <openforcefield.topology.Molecule.to_file>` can now correctly write to `MOL` files in
   line with the support file type list.
- `Issue #568 <https://github.com/openforcefield/openforcefield/issues/568>`_: The
  :py:meth:`Molecule.to_file <openforcefield.topology.Molecule.to_file>` can now correctly write multi-model PDB files
  when using the RDKit backend toolkit.

Example added
"""""""""""""
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_: Adds an example notebook
  `QCarchive_interface.ipynb <https://github.com/openforcefield/openforcefield/blob/master/examples/QCArchive_interface/QCarchive_interface.ipynb>`_
  which shows users how to instance the :py:class:`Molecule <openforcefield.topology.Molecule>` from
  a QCArchive entry level record and calculate the energy using RDKit through QCEngine.

0.6.0 - Library Charges
-----------------------

This release adds support for a new SMIRKS-based charge assignment method,
`Library Charges <https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html#librarycharges-library-charges-for-polymeric-residues-and-special-solvent-models>`_.
The addition of more charge assignment methods opens the door for new types of
experimentation, but also introduces several complex behaviors and failure modes.
Accordingly, we have made changes
to the charge assignment infrastructure to check for cases when partial charges do
not sum to the formal charge of the molecule, or when no charge assignment method is able
to generate charges for a molecule. More detailed explanation of the new errors that may be raised and
keywords for overriding them are in the "Behavior Changed" section below.


With this release, we update ``test_forcefields/tip3p.offxml`` to be a working example of assigning LibraryCharges.
However, we do not provide any force field files to assign protein residue ``LibraryCharges``.
If you are interested in translating an existing protein FF to SMIRNOFF format or developing a new one, please
feel free to contact us on the `Issue tracker <https://github.com/openforcefield/openforcefield/issues>`_ or open a
`Pull Request <https://github.com/openforcefield/openforcefield/pulls>`_.


New features
""""""""""""
- `PR #433 <https://github.com/openforcefield/openforcefield/pull/433>`_: Closes
  `Issue #25 <https://github.com/openforcefield/openforcefield/issues/25>`_ by adding
  initial support for the
  `LibraryCharges tag in the SMIRNOFF specification <https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html#librarycharges-library-charges-for-polymeric-residues-and-special-solvent-models>`_
  using
  :py:class:`LibraryChargeHandler <openforcefield.typing.engines.smirnoff.parameters.LibraryChargeHandler>`.
  For a molecule to have charges assigned using Library Charges, all of its atoms must be covered by
  at least one ``LibraryCharge``. If an atom is covered by multiple ``LibraryCharge`` s, then the last
  ``LibraryCharge`` matched will be applied (per the hierarchy rules in the SMIRNOFF format).

  This functionality is thus able to apply per-residue charges similar to those in traditional
  protein force fields. At this time, there is no concept of "residues" or "fragments" during
  parametrization, so it is not possible to assign charges to `some` atoms in a molecule using
  ``LibraryCharge`` s, but calculate charges for other atoms in the same molecule using a different
  method. To assign charges to a protein, LibraryCharges SMARTS must be provided for
  the residues and protonation states in the molecule, as well as for any capping groups
  and post-translational modifications that are present.

  It is valid for ``LibraryCharge`` SMARTS to `partially` overlap one another. For example, a molecule
  consisting of atoms ``A-B-C`` connected by single bonds could be matched by a SMIRNOFF
  ``LibraryCharges`` section containing two ``LibraryCharge`` SMARTS: ``A-B`` and ``B-C``. If
  listed in that order, the molecule would be assigned the ``A`` charge from the ``A-B`` ``LibraryCharge``
  element and the ``B`` and ``C`` charges from the ``B-C`` element. In testing, these types of
  partial overlaps were found to frequently be sources of undesired behavior, so it is recommended
  that users define whole-molecule ``LibraryCharge`` SMARTS whenever possible.

- `PR #455 <https://github.com/openforcefield/openforcefield/pull/455>`_: Addresses
  `Issue #393 <https://github.com/openforcefield/openforcefield/issues/393>`_ by adding
  :py:meth:`ParameterHandler.attribute_is_cosmetic <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler.attribute_is_cosmetic>`
  and
  :py:meth:`ParameterType.attribute_is_cosmetic <openforcefield.typing.engines.smirnoff.parameters.ParameterType.attribute_is_cosmetic>`,
  which return True if the provided attribute name is defined for the queried object
  but does not correspond to an allowed value in the SMIRNOFF spec.

Behavior changed
""""""""""""""""
- `PR #433 <https://github.com/openforcefield/openforcefield/pull/433>`_: If a molecule
  can not be assigned charges by any charge-assignment method, an
  ``openforcefield.typing.engines.smirnoff.parameters.UnassignedMoleculeChargeException``
  will be raised. Previously, creating a system without either ``ToolkitAM1BCCHandler`` or
  the ``charge_from_molecules`` keyword argument to ``ForceField.create_openmm_system`` would
  produce a system where the molecule has zero charge on all atoms. However, given that we
  will soon be adding more options for charge assignment, it is important that
  failures not be silent. Molecules with zero charge can still be produced by setting the
  ``Molecule.partial_charges`` array to be all zeroes, and including the molecule in the
  ``charge_from_molecules`` keyword argument to ``create_openmm_system``.
- `PR #433 <https://github.com/openforcefield/openforcefield/pull/433>`_: Due to risks
  introduced by permitting charge assignment using partially-overlapping ``LibraryCharge`` s,
  the toolkit will now raise a
  ``openforcefield.typing.engines.smirnoff.parameters.NonIntegralMoleculeChargeException``
  if the sum of partial charges on a molecule are found to be more than 0.01 elementary charge units
  different than the molecule's formal charge. This exception can be overridden by providing
  the ``allow_nonintegral_charges=True`` keyword argument to ``ForceField.create_openmm_system``.




Tests added
"""""""""""
- `PR #430 <https://github.com/openforcefield/openforcefield/pull/430>`_: Added test for
  Wiberg Bond Order implemented in OpenEye Toolkits. Molecules taken from
  DOI:10.5281/zenodo.3405489 . Added by Sukanya Sasmal.
- `PR #569 <https://github.com/openforcefield/openforcefield/pull/569>`_: Added round-trip tests for more serialization formats (dict, YAML, TOML, JSON, BSON, messagepack, pickle). Note that some are unsupported, but the tests raise the appropriate error.


Bugfixes
""""""""
- `PR #431 <https://github.com/openforcefield/openforcefield/pull/431>`_: Fixes an issue
  where ``ToolkitWrapper`` objects would improperly search for functionality in the
  ``GLOBAL_TOOLKIT_REGISTRY``, even though a specific ``ToolkitRegistry`` was requested for an
  operation.
- `PR #439 <https://github.com/openforcefield/openforcefield/pull/439>`_: Fixes
  `Issue #438 <https://github.com/openforcefield/openforcefield/issues/438>`_, by replacing
  call to NetworkX ``Graph.node`` with call to ``Graph.nodes``, per
  `2.4 migration guide <https://networkx.github.io/documentation/stable/release/release_2.4.html>`_.

Files modified
""""""""""""""
- `PR #433 <https://github.com/openforcefield/openforcefield/pull/433>`_: Updates
  the previously-nonfunctional ``test_forcefields/tip3p.offxml`` to a functional state
  by updating it to the SMIRNOFF
  0.3 specification, and specifying atomic charges using the ``LibraryCharges`` tag.


0.5.1 - Adding the parameter coverage example notebook
------------------------------------------------------

This release contains a new notebook example,
`check_parameter_coverage.ipynb <https://github.com/openforcefield/openforcefield/blob/master/examples/check_dataset_parameter_coverage/check_parameter_coverage.ipynb>`_,
which loads sets of molecules, checks whether they are parameterizable,
and generates reports of chemical motifs that are not.
It also fixes several simple issues, improves warnings and docstring text,
and removes unused files.

The parameter coverage example notebook goes hand-in-hand with the
release candidate of our initial force field,
`openff-1.0.0-RC1.offxml <https://github.com/openforcefield/openforcefields>`_
, which will be temporarily available until the official force
field release is made in October.
Our goal in publishing this notebook alongside our first major refitting is to allow interested
users to check whether there is parameter coverage for their molecules of interest.
If the force field is unable to parameterize a molecule, this notebook will generate
reports of the specific chemistry that is not covered. We understand that many organizations
in our field have restrictions about sharing specific molecules, and the outputs from this
notebook can easily be cropped to communicate unparameterizable chemistry without revealing
the full structure.

The force field release candidate is in our new refit force field package,
`openforcefields <https://github.com/openforcefield/openforcefields>`_.
This package is now a part of the Open Force Field Toolkit conda recipe, along with the original
`smirnoff99Frosst <https://github.com/openforcefield/smirnoff99Frosst>`_ line of force fields.

Once the ``openforcefields`` conda package is installed, you can load the release candidate using:

``ff = ForceField('openff-1.0.0-RC1.offxml')``

The release candidate will be removed when the official force field,
``openff-1.0.0.offxml``, is released in early October.

Complete details about this release are below.

Example added
"""""""""""""
- `PR #419 <https://github.com/openforcefield/openforcefield/pull/419>`_: Adds
  an example notebook
  `check_parameter_coverage.ipynb <https://github.com/openforcefield/openforcefield/blob/master/examples/check_dataset_parameter_coverage/check_parameter_coverage.ipynb>`_
  which shows how to use the toolkit to check a molecule
  dataset for missing parameter coverage, and provides functionality to output
  tagged SMILES and 2D drawings of the unparameterizable chemistry.


New features
""""""""""""
- `PR #419 <https://github.com/openforcefield/openforcefield/pull/419>`_: Unassigned
  valence parameter exceptions now include a list of tuples of
  :py:class:`TopologyAtom <openforcefield.topology.TopologyAtom>`
  which were unable to be parameterized (``exception.unassigned_topology_atom_tuples``)
  and the class of the
  :py:class:`ParameterHandler <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler>`
  that raised the exception (``exception.handler_class``).
- `PR #425 <https://github.com/openforcefield/openforcefield/pull/425>`_: Implements
  Trevor Gokey's suggestion from
  `Issue #411 <https://github.com/openforcefield/openforcefield/issues/411>`_, which
  enables pickling of
  :py:class:`ForceFields <openforcefield.typing.engines.smirnoff.forcefield.ForceField>`
  and
  :py:class:`ParameterHandlers <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler>`.
  Note that, while XML representations of ``ForceField``s are stable and conform to the SMIRNOFF
  specification, the pickled ``ForceField``s that this functionality enables are not guaranteed
  to be compatible with future toolkit versions.

Improved documentation and warnings
"""""""""""""""""""""""""""""""""""
- `PR #425 <https://github.com/openforcefield/openforcefield/pull/425>`_: Addresses
  `Issue #410 <https://github.com/openforcefield/openforcefield/issues/410>`_, by explicitly
  having toolkit warnings print ``Warning:`` at the beginning of each warning, and adding
  clearer language to the warning produced when the OpenEye Toolkits can not be loaded.
- `PR #425 <https://github.com/openforcefield/openforcefield/pull/425>`_: Addresses
  `Issue #421 <https://github.com/openforcefield/openforcefield/issues/421>`_ by
  adding type/shape information to all Molecule partial charge and conformer docstrings.
- `PR #425 <https://github.com/openforcefield/openforcefield/pull/425>`_: Addresses
  `Issue #407 <https://github.com/openforcefield/openforcefield/issues/421>`_ by
  providing a more extensive explanation of why we don't use RDKit's mol2 parser
  for molecule input.

Bugfixes
""""""""
- `PR #419 <https://github.com/openforcefield/openforcefield/pull/419>`_: Fixes
  `Issue #417 <https://github.com/openforcefield/openforcefield/issues/417>`_ and
  `Issue #418 <https://github.com/openforcefield/openforcefield/issues/418>`_, where
  :py:meth:`RDKitToolkitWrapper.from_file <openforcefield.utils.toolkits.RDKitToolkitWrapper.from_file>`
  would disregard the ``allow_undefined_stereo`` kwarg and skip the first molecule
  when reading a SMILES file.


Files removed
"""""""""""""
- `PR #425 <https://github.com/openforcefield/openforcefield/pull/425>`_: Addresses
  `Issue #424 <https://github.com/openforcefield/openforcefield/issues/424>`_ by
  deleting the unused files ``openforcefield/typing/engines/smirnoff/gbsaforces.py``
  and ``openforcefield/tests/test_smirnoff.py``. ``gbsaforces.py`` was only used internally
  and ``test_smirnoff.py`` tested unsupported functionality from before the 0.2.0 release.




0.5.0 - GBSA support and quality-of-life improvements
-----------------------------------------------------

This release adds support for the
`GBSA tag in the SMIRNOFF specification <https://open-forcefield-toolkit.readthedocs.io/en/0.5.0/smirnoff.html#gbsa>`_.
Currently, the ``HCT``, ``OBC1``, and ``OBC2`` models (corresponding to AMBER keywords
``igb=1``, ``2``, and ``5``, respectively) are supported, with the ``OBC2`` implementation being
the most flexible. Unfortunately, systems produced
using these keywords are not yet transferable to other simulation packages via ParmEd, so users are restricted
to using OpenMM to simulate systems with GBSA.

OFFXML files containing GBSA parameter definitions are available,
and can be loaded in addition to existing parameter sets (for example, with the command
``ForceField('test_forcefields/smirnoff99Frosst.offxml', 'test_forcefields/GBSA_OBC1-1.0.offxml')``).
A manifest of new SMIRNOFF-format GBSA files is below.


Several other user-facing improvements have been added, including easier access to indexed attributes,
which are now accessible as ``torsion.k1``, ``torsion.k2``, etc. (the previous access method
``torsion.k`` still works as well). More details of the new features and several bugfixes are listed below.

New features
""""""""""""
- `PR #363 <https://github.com/openforcefield/openforcefield/pull/363>`_: Implements
  :py:class:`GBSAHandler <openforcefield.typing.engines.smirnoff.parameters.GBSAHandler>`,
  which supports the
  `GBSA tag in the SMIRNOFF specification <https://open-forcefield-toolkit.readthedocs.io/en/0.5.0/smirnoff.html#gbsa>`_.
  Currently, only GBSAHandlers with ``gb_model="OBC2"`` support
  setting non-default values for the ``surface_area_penalty`` term (default ``5.4*calories/mole/angstroms**2``),
  though users can zero the SA term for ``OBC1`` and ``HCT`` models by setting ``sa_model="None"``.
  No model currently supports setting ``solvent_radius`` to any value other than ``1.4*angstroms``.
  Files containing experimental SMIRNOFF-format implementations of ``HCT``, ``OBC1``, and ``OBC2`` are
  included with this release (see below). Additional details of these models, including literature references,
  are available on the
  `SMIRNOFF specification page <https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html#supported-generalized-born-gb-models>`_.

    .. warning :: The current release of ParmEd
      `can not transfer GBSA models produced by the Open Force Field Toolkit
      to other simulation packages
      <https://github.com/ParmEd/ParmEd/blob/3.2.0/parmed/openmm/topsystem.py#L148-L150>`_.
      These GBSA forces are currently only computable using OpenMM.

- `PR #363 <https://github.com/openforcefield/openforcefield/pull/363>`_: When using
  :py:meth:`Topology.to_openmm() <openforcefield.topology.Topology.to_openmm>`, periodic
  box vectors are now transferred from the Open Force Field Toolkit Topology
  into the newly-created OpenMM Topology.
- `PR #377 <https://github.com/openforcefield/openforcefield/pull/377>`_: Single indexed parameters in
  :py:class:`ParameterHandler <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler>`
  and :py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>`
  can now be get/set through normal attribute syntax in addition to the list syntax.
- `PR #394 <https://github.com/openforcefield/openforcefield/pull/394>`_: Include element and atom name
  in error output when there are missing valence parameters during molecule parameterization.

Bugfixes
""""""""
- `PR #385 <https://github.com/openforcefield/openforcefield/pull/385>`_: Fixes
  `Issue #346 <https://github.com/openforcefield/openforcefield/issues/346>`_ by
  having :py:meth:`OpenEyeToolkitWrapper.compute_partial_charges_am1bcc <openforcefield.utils.toolkits.OpenEyeToolkitWrapper.compute_partial_charges_am1bcc>`
  fall back to using standard AM1-BCC if AM1-BCC ELF10 charge generation raises
  an error about "trans COOH conformers"
- `PR #399 <https://github.com/openforcefield/openforcefield/pull/399>`_: Fixes
  issue where
  :py:class:`ForceField <openforcefield.typing.engines.smirnoff.forcefield.ForceField>`
  constructor would ignore ``parameter_handler_classes`` kwarg.
- `PR #400 <https://github.com/openforcefield/openforcefield/pull/400>`_: Makes
  link-checking tests retry three times before failing.



Files added
"""""""""""
- `PR #363 <https://github.com/openforcefield/openforcefield/pull/363>`_: Adds
  ``test_forcefields/GBSA_HCT-1.0.offxml``, ``test_forcefields/GBSA_OBC1-1.0.offxml``,
  and ``test_forcefields/GBSA_OBC2-1.0.offxml``, which are experimental implementations
  of GBSA models. These are primarily used in validation tests against OpenMM's models, and
  their version numbers will increment if bugfixes are necessary.

0.4.1 - Bugfix Release
----------------------

This update fixes several toolkit bugs that have been reported by the community.
Details of these bugfixes are provided below.

It also refactors how
:py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>`
and
:py:class:`ParameterHandler <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler>`
store their attributes, by introducing
:py:class:`ParameterAttribute <openforcefield.typing.engines.smirnoff.parameters.ParameterAttribute>`
and
:py:class:`IndexedParameterAttribute <openforcefield.typing.engines.smirnoff.parameters.IndexedParameterAttribute>`.
These new attribute-handling classes provide a consistent backend which should simplify manipulation of parameters
and implementation of new handlers.

Bug fixes
"""""""""
- `PR #329 <https://github.com/openforcefield/openforcefield/pull/329>`_: Fixed a
  bug where the two
  :py:class:`BondType <openforcefield.typing.engines.smirnoff.parameters.BondHandler.BondType>`
  parameter attributes ``k`` and ``length`` were treated as indexed attributes. (``k`` and
  ``length`` values that correspond to specific bond orders will be indexed under
  ``k_bondorder1``, ``k_bondorder2``, etc when implemented in the future)
- `PR #329 <https://github.com/openforcefield/openforcefield/pull/329>`_: Fixed a
  bug that allowed setting indexed attributes to single values instead of strictly lists.
- `PR #370 <https://github.com/openforcefield/openforcefield/pull/370>`_: Fixed a
  bug in the API where
  :py:class:`BondHandler <openforcefield.typing.engines.smirnoff.parameters.BondHandler>`,
  :py:class:`ProperTorsionHandler <openforcefield.typing.engines.smirnoff.parameters.ProperTorsionHandler>`
  , and
  :py:class:`ImproperTorsionHandler <openforcefield.typing.engines.smirnoff.parameters.ImproperTorsionHandler>`
  exposed non-functional indexed parameters.
- `PR #351 <https://github.com/openforcefield/openforcefield/pull/351>`_: Fixes
  `Issue #344 <https://github.com/openforcefield/openforcefield/issues/344>`_,
  in which the main :py:class:`FrozenMolecule <openforcefield.topology.FrozenMolecule>`
  constructor and several other Molecule-construction functions ignored or did not
  expose the ``allow_undefined_stereo`` keyword argument.
- `PR #351 <https://github.com/openforcefield/openforcefield/pull/351>`_: Fixes
  a bug where a molecule which previously generated a SMILES using one cheminformatics toolkit
  returns the same SMILES, even though a different toolkit (which would generate
  a different SMILES for the molecule) is explicitly called.
- `PR #354 <https://github.com/openforcefield/openforcefield/pull/354>`_: Fixes
  the error message that is printed if an unexpected parameter attribute is found while loading
  data into a :py:class:`ForceField <openforcefield.typing.engines.smirnoff.forcefield.ForceField>`
  (now instructs users to specify ``allow_cosmetic_attributes`` instead of ``permit_cosmetic_attributes``)
- `PR #364 <https://github.com/openforcefield/openforcefield/pull/364>`_: Fixes
  `Issue #362 <https://github.com/openforcefield/openforcefield/issues/362>`_ by
  modifying
  :py:meth:`OpenEyeToolkitWrapper.from_smiles <openforcefield.utils.toolkits.OpenEyeToolkitWrapper.from_smiles>`
  and
  :py:meth:`RDKitToolkitWrapper.from_smiles <openforcefield.utils.toolkits.RDKitToolkitWrapper.from_smiles>`
  to make implicit hydrogens explicit before molecule creation. These functions also
  now raise an error if the optional keyword ``hydrogens_are_explicit=True`` but the
  SMILES are interpreted by the backend cheminformatic toolkit as having implicit
  hydrogens.
- `PR #371 <https://github.com/openforcefield/openforcefield/pull/371>`_: Fixes
  error when reading early SMIRNOFF 0.1 spec files enclosed by a top-level ``SMIRFF`` tag.

.. note ::
  The enclosing ``SMIRFF`` tag is present only in legacy files.
  Since developing a formal specification, the only acceptable top-level tag value in a SMIRNOFF data structure is
  ``SMIRNOFF``.

Code enhancements
"""""""""""""""""
- `PR #329 <https://github.com/openforcefield/openforcefield/pull/329>`_:
  :py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>`
  was refactored to improve its extensibility. It is now possible to create new parameter
  types by using the new descriptors
  :py:class:`ParameterAttribute <openforcefield.typing.engines.smirnoff.parameters.ParameterAttribute>`
  and
  :py:class:`IndexedParameterAttribute <openforcefield.typing.engines.smirnoff.parameters.IndexedParameterAttribute>`.
- `PR #357 <https://github.com/openforcefield/openforcefield/pull/357>`_: Addresses
  `Issue #356 <https://github.com/openforcefield/openforcefield/issues/356>`_ by raising
  an informative error message if a user attempts to load an OpenMM topology which
  is probably missing connectivity information.



Force fields added
""""""""""""""""""
- `PR #368 <https://github.com/openforcefield/openforcefield/pull/368>`_: Temporarily adds
  ``test_forcefields/smirnoff99frosst_experimental.offxml`` to address hierarchy problems, redundancies, SMIRKS
  pattern typos etc., as documented in `issue #367 <https://github.com/openforcefield/openforcefield/issues/367>`_.
  Will ultimately be propagated to an updated forcefield in the ``openforcefield/smirnoff99frosst`` repo.
- `PR #371 <https://github.com/openforcefield/openforcefield/pull/371>`_: Adds
  ``test_forcefields/smirff99Frosst_reference_0_1_spec.offxml``, a SMIRNOFF 0.1 spec file enclosed by the legacy
  ``SMIRFF`` tag. This file is used in backwards-compatibility testing.



0.4.0 - Performance optimizations and support for SMIRNOFF 0.3 specification
----------------------------------------------------------------------------

This update contains performance enhancements that significantly reduce the time to create OpenMM systems for topologies containing many molecules via :py:meth:`ForceField.create_openmm_system <openforcefield.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system>`.

This update also introduces the `SMIRNOFF 0.3 specification <https://open-forcefield-toolkit.readthedocs.io/en/0.4.0/smirnoff.html>`_.
The spec update is the result of discussions about how to handle the evolution of data and parameter types as further functional forms are added to the SMIRNOFF spec.


We provide methods to convert SMIRNOFF 0.1 and 0.2 forcefields written with the XML serialization (``.offxml``) to the SMIRNOFF 0.3 specification.
These methods are called automatically when loading a serialized SMIRNOFF data representation written in the 0.1 or 0.2 specification.
This functionality allows the toolkit to continue to read files containing SMIRNOFF 0.2 spec force fields, and also implements backwards-compatibility for SMIRNOFF 0.1 spec force fields.


.. warning :: The SMIRNOFF 0.1 spec did not contain fields for several energy-determining parameters that are exposed in later SMIRNOFF specs.
  Thus, when reading SMIRNOFF 0.1 spec data, the toolkit must make assumptions about the values that should be added for the newly-required fields.
  The values that are added include 1-2, 1-3 and 1-5 scaling factors, cutoffs, and long-range treatments for nonbonded interactions.
  Each assumption is printed as a warning during the conversion process.
  Please carefully review the warning messages to ensure that the conversion is providing your desired behavior.



`SMIRNOFF 0.3 specification updates <https://open-forcefield-toolkit.readthedocs.io/en/0.4.0/smirnoff.html>`_
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
* The SMIRNOFF 0.3 spec introduces versioning for each individual parameter section, allowing asynchronous updates to the features of each parameter class.
  The top-level ``SMIRNOFF`` tag, containing information like ``aromaticity_model``, ``Author``, and ``Date``, still has a version (currently 0.3).
  But, to allow for independent development of individual parameter types, each section (such as ``Bonds``, ``Angles``, etc) now has its own version as well (currently all 0.3).
* All units are now stored in expressions with their corresponding values. For example, distances are now stored as ``1.526*angstrom``, instead of storing the unit separately in the section header.
* The current allowed value of the ``potential`` field for ``ProperTorsions`` and ``ImproperTorsions`` tags is no longer ``charmm``, but is rather ``k*(1+cos(periodicity*theta-phase))``.
  It was pointed out to us that CHARMM-style torsions deviate from this formula when the periodicity of a torsion term is 0, and we do not intend to reproduce that behavior.
* SMIRNOFF spec documentation has been updated with tables of keywords and their defaults for each parameter section and parameter type.
  These tables will track the allowed keywords and default behavior as updated versions of individual parameter sections are released.

Performance improvements and bugfixes
"""""""""""""""""""""""""""""""""""""

* `PR #329 <https://github.com/openforcefield/openforcefield/pull/329>`_: Performance improvements when creating systems for topologies with many atoms.
* `PR #347 <https://github.com/openforcefield/openforcefield/pull/347>`_: Fixes bug in charge assignment that occurs when charges are read from file, and reference and charge molecules have different atom orderings.


New features
""""""""""""

* `PR #311 <https://github.com/openforcefield/openforcefield/pull/311>`_: Several new experimental functions.

  * Adds :py:meth:`convert_0_2_smirnoff_to_0_3 <openforcefield.utils.utils.convert_0_2_smirnoff_to_0_3>`, which takes a SMIRNOFF 0.2-spec data dict, and updates it to 0.3.
    This function is called automatically when creating a ``ForceField`` from a SMIRNOFF 0.2 spec OFFXML file.
  * Adds :py:meth:`convert_0_1_smirnoff_to_0_2 <openforcefield.utils.utils.convert_0_1_smirnoff_to_0_2>`, which takes a SMIRNOFF 0.1-spec data dict, and updates it to 0.2.
    This function is called automatically when creating a ``ForceField`` from a SMIRNOFF 0.1 spec OFFXML file.
  * NOTE: The format of the "SMIRNOFF data dict" above is likely to change significantly in the future.
    Users that require a stable serialized ForceField object should use the output of :py:meth:`ForceField.to_string('XML') <openforcefield.typing.engines.smirnoff.forcefield.ForceField.to_string>` instead.
  * Adds :py:class:`ParameterHandler <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler>` and :py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>` :py:meth:`add_cosmetic_attribute <openforcefield.typing.engines.smirnoff.parameters.ParameterType.add_cosmetic_attribute>` and :py:meth:`delete_cosmetic_attribute <openforcefield.typing.engines.smirnoff.parameters.ParameterType.delete_cosmetic_attribute>` functions.
    Once created, cosmetic attributes can be accessed and modified as attributes of the underlying object (eg. ``ParameterType.my_cosmetic_attrib = 'blue'``)
    These functions are experimental, and we are interested in feedback on how cosmetic attribute handling could be improved. (`See Issue #338 <https://github.com/openforcefield/openforcefield/issues/338>`_)
    Note that if a new cosmetic attribute is added to an object without using these functions, it will not be recognized by the toolkit and will not be written out during serialization.
  * Values for the top-level ``Author`` and ``Date`` tags are now kept during SMIRNOFF data I/O.
    If multiple data sources containing these fields are read, the values are concatenated using "AND" as a separator.


API-breaking changes
""""""""""""""""""""
* :py:meth:`ForceField.to_string <openforcefield.typing.engines.smirnoff.forcefield.ForceField.to_string>` and :py:meth:`ForceField.to_file <openforcefield.typing.engines.smirnoff.forcefield.ForceField.to_file>` have had the default value of their ``discard_cosmetic_attributes`` kwarg set to False.
* :py:class:`ParameterHandler <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler>` and :py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>` constructors now expect the ``version`` kwarg (per the SMIRNOFF spec change above)
  This requirement can be skipped by providing the kwarg ``skip_version_check=True``
* :py:class:`ParameterHandler <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler>` and :py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>` functions no longer handle ``X_unit`` attributes in SMIRNOFF data (per the SMIRNOFF spec change above).
* The scripts in ``utilities/convert_frosst`` are now deprecated.
  This functionality is important for provenance and will be migrated to the ``openforcefield/smirnoff99Frosst`` repository in the coming weeks.
* :py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>` ``._SMIRNOFF_ATTRIBS`` is now :py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>` ``._REQUIRED_SPEC_ATTRIBS``, to better parallel the structure of the ``ParameterHandler`` class.
* :py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>` ``._OPTIONAL_ATTRIBS`` is now :py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>` ``._OPTIONAL_SPEC_ATTRIBS``, to better parallel the structure of the ``ParameterHandler`` class.
* Added class-level dictionaries :py:class:`ParameterHandler <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler>` ``._DEFAULT_SPEC_ATTRIBS`` and :py:class:`ParameterType <openforcefield.typing.engines.smirnoff.parameters.ParameterType>` ``._DEFAULT_SPEC_ATTRIBS``.

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
