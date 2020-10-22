Release History
===============

Releases follow the ``major.minor.micro`` scheme recommended by `PEP440 <https://www.python.org/dev/peps/pep-0440/#final-releases>`_, where

* ``major`` increments denote a change that may break API compatibility with previous ``major`` releases
* ``minor`` increments add features but do not break API compatibility
* ``micro`` increments represent bugfix releases or improvements in documentation

0.8.0 - Virtual Sites
---------------------

This release implements the SMIRNOFF virtual site specification. The implementation enables support for models using off-site charges, including 4- and 5-point water models, in addition to lone pair modeling on various functional groups. The primary focus was on the ability to parameterize a system using virtual sites, and generating an OpenMM system with all virtual sites present and ready for evaluation. Support for formats other than OpenMM has not be implemented in this release, but may come with the appearance of the OpenFF system object. In addition to implementing the specification, the toolkit :py:class:`Molecule <openforcefield.topology.Molecule>` objects now allow the creation and manipulation of virtual sites.

**Major Feature: Support for the SMIRNOFF VirtualSite tag**

Virtual sites can be added to a System in two ways:

* `SMIRNOFF Force Fields can contain a VirtualSites tag <https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html#virtualsites-virtual-sites-for-off-atom-charges>`_ , specifying the addition of virtual sites according to SMARTS-based rules.
* Virtual sites can be added to a :py:class:`Molecule <openforcefield.topology.Molecule>`, and these will appear in the final OpenMM system if a virtual site handler is present in the :py:class:`ForceField <openforcefield.typing.engines.smirnoff.forcefield.ForceField>`.

Virtual sites are the first parameters which directly depend on 3D conformation, where the position of the virtual sites are based on vectors defined on the atoms that were matched during parameterization. Because of this, a virtual site matching the triplet of atoms 1-2-3 will define a point that is different from a triplet matching 3-2-1. This is similar to defining "right-handed" and "left-handed" coordinate systems. This subtlety interplays with two major concepts in force field development:

1) we sometimes want to define a single virtual site describing two points with the same parameters (distance, angle, etc.), such as 5-point water models
2) we have a match that produces multiple orderings of the atoms (e.g. if wildcards are present in the SMARTS pattern), and we only want one to be applied.

Case 1) is very useful for parameter optimization, where a single SMARTS-based parameter can be used to optimize both points, such as the angle defining the virtual points for a 5-point water model. Case 2) is the typical scenario for the nitrogen lone pair in ammonia, where only one point needs to be specified. We discuss a few more illustrative examples below. Beyond these attributes, the virtual site specification allows a policy for specifying how to handle exclusions in the OpenMM force evaluator. The current default is to add pairwise energy exclusions in the OpenMM system between a virtual site and all tagged atoms matched in its SMARTS (``exclusion_policy="parents",`` ). Currently defined are ``"none"``, ``"minimal"``, and ``"parents"``, where ``"minimal"`` specifies the single atom that the virtual site defines as the "origin". For water, for example, ``"minimal"`` would mean just the oxygen, whereas ``"parents"`` would mean all three atoms.

In order to give consistent and intended behavior, the specification was modified from its draft form in following manner: The ``"name"`` and ``"match"`` attributes have been added to each virtual site parameter type. These changes allow for

* specifying different virtual site types using the same atoms
* allowing two virtual sites with the same type and same atoms but different physical parameters to be added simultaneously
* allowing the ability to control whether the virtual site encodes one or multiple particles, based on the number of ways the matching atoms can be ordered.

The ``"name"`` attribute encodes whether the virtual site to be added should override an existing virtual site of the same type (e.g. hierarchy preference), or if this virtual site should be added in addition to the other existing virtual sites on the given atoms. This means that different virtual site types can share the same group of parent atoms and use the same name without overwriting each other (the default ``name`` is ``EP`` for all sites, which gives the expected hierarchical behavior used in other SMIRNOFF tags).

The ``"match"`` attribute accepts either ``"once"`` or ``"all_permutations"``, offering control for situations where a SMARTS pattern can possibly match the same group of atoms in different orders (either due to wildcards or local symmetry) and it is desired to either add just one or all of the possible virtual particles. The default value is ``"all_permutations",`` but for TrivalentLonePair it is always set to ``"once"``, regardless of what the file contains, since all orderings always place the particle in the exact same position.

The following cases exemplify our reasoning in implementing this behavior, and should draw caution to complex issues that may arise when designing virtual site parameters. Let us consider 4-, 5-, and 6-point water models:

* A 4-point water model with a ``DivalentLonePair``: This can be implemented by specifying ``match="once"``, ``outOfPlaneAngle="0*degree"``, and ``distance=-.15*angstrom"``. Since the SMIRKS pattern ``"[#1:1]-[#8X2:2]-[#2:3]"`` would match water twice and would create two particles in the exact same position if ``all_permutations`` was specified, we specify ``"once"`` to have only one particle generated. Although having two particles in the same position should not affect the physics if the proper exclusion policy is applied, it would effectively make the 4-point model just as expensive as 5-point models.

* A 5-point water model with a ``DivalentLonePair``: This can be implemented by using ``match="all_permutations"`` (unlike the 4-point model), ``outOfPlaneAngle="56.26*degree``, and ``distance=0.7*angstrom``, for example. Here the permutations will cause particles to be placed at ±56.26 degrees, and changing any of the physical quantities will affect *both* particles.

* A 6-point water model with both ``DivalentLonePair`` sites above. Since these two parameters look identical, it is unclear whether they should both be applied or if one should override the other. The toolkit never compares the physical numbers to determine equality as this can lead to instability during e.g. parameter fitting. To get this to work, we specify ``name="EP1"`` for the first parameter, and ``name="EP2"`` for the second parameter. This instructs the parameter handler keep them separate, and therefore both are applied. (If both had the same name, then the typical SMIRNOFF hierarchy rules are used, and only the last matched parameter would be applied.)

* Dinitrogen, ``N#N`` with a ``BondCharge`` virtual site. Since we want a ``BondCharge`` on both ends, we specify ``match="all_permutations"``.

* Formaldehyde, ``H2C=O``, with ``MonovalentLonePair`` virtual site(s) on the oxygen, with the aim of modeling both lone pairs. This one is subtle, since ``[#1:3]-[#6X3:2]=[#8X1:1]`` matches two unique groups of atoms (``1-3-4`` and ``2-3-4``). It is important to note in this situation that ``match="all_permutations"`` behaves exactly the same as ``match="once"``. Due to the anchoring hydrogens (``1`` and ``2``) being symmetric but opposite about the bond between ``3`` and ``4``, a single parameter does correctly place both lone pairs. A standing issue here is that the default exclusion policy (``parents``) will allow these two virtual sites to interact since they have different indexed atoms (parents), causing the energy to be different than the non-virtual site parameterization. In the future, the ``exclusion_policy="local"`` will account for this, and make virtual sites that share at least one "parent" atom not interact with each other. As a special note: when applying a ``MonovalentLonePair`` to a completely symmetric molecule, e.g. water, ``all_permutations`` can come into play, but this will apply two particles (one for each hydrogen).

Finally, the toolkit handles the organization of atoms and virtual sites in a specific manner. Virtual sites are expected to be added *after all molecules in the topology are present*. This is because the Open Force Field Toolkit organizes a topology by placing all atoms first, then all virtual sites last. This differs from the OpenMM Modeller object, for example, which interleaves the order of atoms and virtual sites in such a way that all particles of a molecule are contiguous. In addition, due to the fact that a virtual site may contain multiple particles coupled to single parameters, the toolkit makes a distinction between a virtual *site*, and a virtual *particle*. A virtual site may represent multiple virtual particles, so the total number of particles cannot be directly determined by simply summing the number of atoms and virtual sites in a molecule. This is taken into account, however, and the :py:class:`Molecule <openforcefield.topology.Molecule>` and :py:class:`Topology <openforcefield.topology.Topology>` classes now implement ``particle`` iterators.


**Minor Feature: Support for the 0.4 ChargeIncrementModel tag**

To allow for more convenient fitting of ``ChargeIncrement`` parameters, it is now possible to specify one less
``charge_increment`` value than there are tagged atoms in a ``ChargeIncrement``'s ``smirks``. The missing
``charge_increment`` value will be calculated at parameterization-time to make the sum of
the charge contributions from a ``ChargeIncrement`` parameter equal to zero.
Since this change allows for force fields that are incompatible with
the previous specification, this new style of ``ChargeIncrement`` must specify a ``ChargeIncrementModel``
section version of ``0.4``. All ``0.3``-compatible ``ChargeIncrement`` parameters are compatible with
the ``0.4`` ``ChargeIncrementModel`` specification.

More details and examples of this change are available in `The ChargeIncrementModel tag in the SMIRNOFF specification <https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html#chargeincrementmodel-small-molecule-and-fragment-charges>`_


New features
""""""""""""
- `PR #726 <https://github.com/openforcefield/openforcefield/pull/726>`_: Adds support for the 0.4
  ChargeIncrementModel spec, allowing for the specification of one fewer ``charge_increment`` values
  than there are tagged atoms in the ``smirks``, and automatically assigning the final atom an offsetting charge.
- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Adds support for the ``VirtualSites`` tag in the SMIRNOFF specification

- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Adds ``replace`` and ``all_permutations`` kwarg to

  - :py:meth:`Molecule.add_bond_charge_virtual_site <openforcefield.topology.Molecule.add_bond_charge_virtual_site>`
  - :py:meth:`Molecule.add_monovalent_lone_pair_virtual_site <openforcefield.topology.Molecule.add_monovalent_lone_pair_virtual_site>`
  - :py:meth:`Molecule.add_divalent_lone_pair_virtual_site <openforcefield.topology.Molecule.add_divalent_lone_pair_virtual_site>`
  - :py:meth:`Molecule.add_trivalent_lone_pair_virtual_site <openforcefield.topology.Molecule.add_trivalent_lone_pair_virtual_site>`

- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Adds ``orientations`` to

  - :py:class:`BondChargeVirtualSite <openforcefield.topology.BondChargeVirtualSite>`
  - :py:class:`MonovalentLonePairVirtualSite <openforcefield.topology.MonovalentLonePairVirtualSite>`
  - :py:class:`DivalentLonePairVirtualSite <openforcefield.topology.DivalentLonePairVirtualSite>`
  - :py:class:`TrivalentLonePairVirtualSite <openforcefield.topology.TrivalentLonePairVirtualSite>`

- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Adds

  - :py:class:`VirtualParticle <openforcefield.topology.VirtualParticle>`
  - :py:class:`TopologyVirtualParticle <openforcefield.topology.TopologyVirtualParticle>`
  - :py:meth:`BondChargeVirtualSite.get_openmm_virtual_site <openforcefield.topology.BondChargeVirtualSite.get_openmm_virtual_site>`
  - :py:meth:`MonovalentLonePairVirtualSite.get_openmm_virtual_site <openforcefield.topology.MonovalentLonePairVirtualSite.get_openmm_virtual_site>`
  - :py:meth:`DivalentLonePairVirtualSite.get_openmm_virtual_site <openforcefield.topology.DivalentLonePairVirtualSite.get_openmm_virtual_site>`
  - :py:meth:`TrivalentLonePairVirtualSite.get_openmm_virtual_site <openforcefield.topology.TrivalentLonePairVirtualSite.get_openmm_virtual_site>`
  - :py:meth:`ValenceDict.key_transform <openforcefield.topology.ValenceDict.key_transform>`
  - :py:meth:`ValenceDict.index_of <openforcefield.topology.ValenceDict.index_of>`
  - :py:meth:`ImproperDict.key_transform <openforcefield.topology.ImproperDict.key_transform>`
  - :py:meth:`ImproperDict.index_of <openforcefield.topology.ImproperDict.index_of>`

- `PR #705 <https://github.com/openforcefield/openforcefield/pull/705>`_: Adds interpolation
  based on fractional bond orders for harmonic bonds. This includes interpolation for both
  the force constant ``k`` and/or equilibrium bond distance ``length``. This is accompanied by a
  bump in the ``<Bonds>`` section of the SMIRNOFF spec (but not the entire spec).
- `PR #718 <https://github.com/openforcefield/openforcefield/pull/718>`_: Adds ``.rings`` and
  ``.n_rings`` to :py:class:`Molecule <openforcefield.topology.Molecule>` and ``.is_in_ring``
  to :py:class:`Atom <openforcefield.topology.Atom>` and
  :py:class:`Bond <openforcefield.topology.Bond>`

Bugfixes
"""""""""
- `PR #682 <https://github.com/openforcefield/openforcefield/pull/682>`_: Catches failures in
  :py:meth:`Molecule.from_iupac <openforcefield.topology.Molecule.from_iupac>` instead of silently
  failing.
- `PR #743 <https://github.com/openforcefield/openforcefield/pull/743>`_: Prevents the non-bonded
  (vdW) cutoff from silently falling back to the OpenMM default of 1 nm in
  :py:meth:`Forcefield.create_openmm_system
  <openforcefield.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system>` and instead
  sets its to the value specified by the force field.
- `PR #737 <https://github.com/openforcefield/openforcefield/pull/737>`_: Prevents OpenEye from
  incidentally being used in the conformer generation step of
  :py:class:`AmberToolsToolkitWrapper.assign_fractional_bond_orders
  <openforcefield.utils.toolkits.AmberToolsToolkitWrapper.assign_fractional_bond_orders>`.

Behavior changed
""""""""""""""""
- `PR #705 <https://github.com/openforcefield/openforcefield/pull/705>`_: Changes the default values
  in the ``<Bonds>`` section of the SMIRNOFF spec to ``fractional_bondorder_method="AM1-Wiberg"``
  and ``potential="(k/2)*(r-length)^2"``, which is backwards-compatible with and equivalent to
  ``potential="harmonic"``.

Examples added
""""""""""""""
- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Adds a virtual site example notebook to run
  an OpenMM simulation with virtual sites, and compares positions and potential energy of TIP5P water between OpenFF
  and OpenMM forcefields.

API-breaking changes
""""""""""""""""""""
- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Methods

  - :py:meth:`Molecule.add_bond_charge_virtual_site <openforcefield.topology.Molecule.add_bond_charge_virtual_site>`
  - :py:meth:`Molecule.add_monovalent_lone_pair_virtual_site <openforcefield.topology.Molecule.add_monovalent_lone_pair_virtual_site>`
  - :py:meth:`Molecule.add_divalent_lone_pair_virtual_site <openforcefield.topology.Molecule.add_divalent_lone_pair_virtual_site>`
  - :py:meth:`Molecule.add_trivalent_lone_pair_virtual_site <openforcefield.topology.Molecule.add_trivalent_lone_pair_virtual_site>`
  now only accept a list of atoms, not a list of integers, to define to parent atoms

- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Removes
  :py:meth:`VirtualParticle.molecule_particle_index <openforcefield.topology.VirtualParticle.molecule_particle_index>`

- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Removes ``outOfPlaneAngle`` from

  - :py:class:`DivalentLonePairVirtualSite <openforcefield.topology.DivalentLonePairVirtualSite>`
  - :py:class:`TrivalentLonePairVirtualSite <openforcefield.topology.TrivalentLonePairVirtualSite>`

- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Removes ``inPlaneAngle`` from
  :py:class:`TrivalentLonePairVirtualSite <openforcefield.topology.TrivalentLonePairVirtualSite>`

- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Removes ``weights`` from

  - :py:class:`BondChargeVirtualSite <openforcefield.topology.BondChargeVirtualSite>`
  - :py:class:`MonovalentLonePairVirtualSite <openforcefield.topology.MonovalentLonePairVirtualSite>`
  - :py:class:`DivalentLonePairVirtualSite <openforcefield.topology.DivalentLonePairVirtualSite>`
  - :py:class:`TrivalentLonePairVirtualSite <openforcefield.topology.TrivalentLonePairVirtualSite>`

Tests added
"""""""""""

- `PR #548 <https://github.com/openforcefield/openforcefield/pull/548>`_: Adds test for 

  - The virtual site parameter handler
  - TIP5P water dimer energy and positions
  - Adds tests to for virtual site/particle indexing/counting


0.7.2 - Bugfix and minor feature release
----------------------------------------

New features
""""""""""""
- `PR #662 <https://github.com/openforcefield/openforcefield/pull/662>`_: Adds ``.aromaticity_model``
  of :py:class:`ForceField <openforcefield.typing.engines.smirnoff.forcefield.ForceField>` and ``.TAGNAME``
  of :py:class:`ParameterHandler <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler>` as
  public attributes.
- `PR #667 <https://github.com/openforcefield/openforcefield/pull/667>`_ and
  `PR #681 <https://github.com/openforcefield/openforcefield/pull/681>`_ linted the codebase with
  ``black`` and ``isort``, respectively.
- `PR #675 <https://github.com/openforcefield/openforcefield/pull/675>`_ adds
  ``.toolkit_version`` to
  :py:class:`ToolkitWrapper <openforcefield.utils.toolkits.ToolkitWrapper>` and
  ``.registered_toolkit_versions`` to
  :py:class:`ToolkitRegistry <openforcefield.utils.toolkits.ToolkitRegistry>`.
- `PR #696 <https://github.com/openforcefield/openforcefield/pull/696>`_ Exposes a setter for
  :py:class:`ForceField.aromaticity_model <openforcefield.typing.engines.smirnoff.forcefield.ForceField>`
- `PR #685 <https://github.com/openforcefield/openforcefield/pull/685>`_ Adds a custom ``__hash__``
  function to
  :py:class:`ForceField <openforcefield.typing.engines.smirnoff.forcefield.ForceField>`


Behavior changed
""""""""""""""""
- `PR #684 <https://github.com/openforcefield/openforcefield/pull/684>`_: Changes
  :py:class:`ToolkitRegistry <openforcefield.utils.toolkits.ToolkitRegistry>` to return an empty
  registry when initialized with no arguments, i.e. ``ToolkitRegistry()`` and makes the
  ``register_imported_toolkit_wrappers`` argument private.
- `PR #711 <https://github.com/openforcefield/openforcefield/pull/711>`_: The
  setter for :py:class:`Topology.box_vectors <openforcefield.topology.Topology>`
  now infers box vectors (a 3x3 matrix) when box lengths
  (a 3x1 array) are passed, assuming an orthogonal box.
- `PR #649 <https://github.com/openforcefield/openforcefield/pull/648>`_: Makes SMARTS
  searches stereochemistry-specific (if stereo is specified in the SMARTS) for both OpenEye
  and RDKit backends. Also ensures molecule
  aromaticity is re-perceived according to the ForceField's specified
  aromaticity model, which may overwrite user-specified aromaticity on the ``Molecule``
- `PR #648 <https://github.com/openforcefield/openforcefield/pull/648>`_: Removes the
  ``utils.structure`` module, which was deprecated in 0.2.0.
- `PR #670 <https://github.com/openforcefield/openforcefield/pull/670>`_: Makes the
  :py:class:`Topology <openforcefield.topology.Topology>` returned by ``create_openmm_system``
  contain the partial charges and partial bond orders (if any) assigned during parameterization.
- `PR #675 <https://github.com/openforcefield/openforcefield/pull/675>`_ changes the
  exception raised when no ``antechamber`` executable is found from ``IOError`` to
  :py:class:`AntechamberNotFoundError <openforcefield.utils.toolkits.AntechamberNotFoundError>`
- `PR #696 <https://github.com/openforcefield/openforcefield/pull/696>`_ Adds an
  ``aromaticity_model`` keyword argument to the
  :py:class:`ForceField <openforcefield.typing.engines.smirnoff.forcefield.ForceField>`
  constructor, which defaults to ``DEFAULT_AROMATICITY_MODEL``.

Bugfixes
"""""""""
- `PR #715 <https://github.com/openforcefield/openforcefield/pull/715>`_: Closes issue `Issue #475
  <https://github.com/openforcefield/openforcefield/issues/475>`_ writing a "PDB" file using OE backend rearranges
  the order of the atoms by pushing the hydrogens to the bottom.
- `PR #649 <https://github.com/openforcefield/openforcefield/pull/648>`_: Prevents 2020 OE
  toolkit from issuing a warning caused by doing stereo-specific smarts searches on certain
  structures.
- `PR #724 <https://github.com/openforcefield/openforcefield/pull/724>`_: Closes issue `Issue #502
  <https://github.com/openforcefield/openforcefield/issues/502>`_ Adding a utility function Topology.to_file() to 
  write topology and positions to a "PDB" file using openmm backend for pdb file write.

Tests added
"""""""""""
- `PR #694 <https://github.com/openforcefield/openforcefield/pull/694>`_: Adds automated testing
  to code snippets in docs.
- `PR #715 <https://github.com/openforcefield/openforcefield/pull/715>`_: Adds tests for pdb file writes using OE
  backend.
- `PR #724 <https://github.com/openforcefield/openforcefield/pull/724>`_: Adds tests for the utility function Topology.to_file().
  

0.7.1 - OETK2020 Compatibility and Minor Update
-----------------------------------------------

This is the first of our patch releases on our new planned monthly release schedule.

Detailed release notes are below, but the major new features of this release are updates for
compatibility with the new 2020 OpenEye Toolkits release, the
``get_available_force_fields`` function, and the disregarding of pyrimidal nitrogen stereochemistry
in molecule isomorphism checks.

Behavior changed
""""""""""""""""
- `PR #646 <https://github.com/openforcefield/openforcefield/pull/646>`_: Checking for
  :py:class:`Molecule <openforcefield.topology.Molecule>`
  equality using the ``==`` operator now disregards all pyrimidal nitrogen stereochemistry
  by default. To re-enable, use
  :py:class:`Molecule.{is|are}_isomorphic <openforcefield.topology.Molecule>`
  with the ``strip_pyrimidal_n_atom_stereo=False`` keyword argument.
- `PR #646 <https://github.com/openforcefield/openforcefield/pull/646>`_: Adds
  an optional ``toolkit_registry`` keyword argument to
  :py:class:`Molecule.are_isomorphic <openforcefield.topology.Molecule>`,
  which identifies the toolkit that should be used to search for pyrimidal nitrogens.


Bugfixes
""""""""
- `PR #647 <https://github.com/openforcefield/openforcefield/pull/647>`_: Updates
  :py:class:`OpenEyeToolkitWrapper <openforcefield.utils.toolkits.OpenEyeToolkitWrapper>`
  for 2020.0.4 OpenEye Toolkit behavior/API changes.
- `PR #646 <https://github.com/openforcefield/openforcefield/pull/646>`_: Fixes a bug where
  :py:class:`Molecule.chemical_environment_matches <openforcefield.topology.Molecule>`
  was not able to accept a :py:class:`ChemicalEnvironment <openforcefield.typing.chemistry.ChemicalEnvironment>` object
  as a query.
- `PR #634 <https://github.com/openforcefield/openforcefield/pull/634>`_: Fixes a bug in which calling
  :py:class:`RDKitToolkitWrapper.from_file <openforcefield.utils.toolkits.RDKitToolkitWrapper>` directly
  would not load files correctly if passed lowercase ``file_format``. Note that this bug did not occur when calling
  :py:class:`Molecule.from_file <openforcefield.topology.Molecule>`.
- `PR #631 <https://github.com/openforcefield/openforcefield/pull/631>`_: Fixes a bug in which calling
  :py:class:`unit_to_string <openforcefield.utils.utils.unit_to_string>` returned
  ``None`` when the unit is dimensionless. Now ``"dimensionless"`` is returned.
- `PR #630 <https://github.com/openforcefield/openforcefield/pull/630>`_: Closes issue `Issue #629 
  <https://github.com/openforcefield/openforcefield/issues/629>`_ in which the wrong exception is raised when
  attempting to instantiate a :py:class:`ForceField <openforcefield.typing.engines.smirnoff.forcefield.ForceField>`
  from an unparsable string.

New features
""""""""""""
- `PR #632 <https://github.com/openforcefield/openforcefield/pull/632>`_: Adds
  :py:class:`ForceField.registered_parameter_handlers <openforcefield.typing.engines.smirnoff.forcefield.ForceField>`
- `PR #614 <https://github.com/openforcefield/openforcefield/pull/614>`_: Adds 
  :py:class:`ToolkitRegistry.deregister_toolkit <openforcefield.utils.toolkits.ToolkitRegistry>`
  to de-register registered toolkits, which can include toolkit wrappers loaded into ``GLOBAL_TOOLKIT_REGISTRY``
  by default.
- `PR #656 <https://github.com/openforcefield/openforcefield/pull/656>`_: Adds
  a new allowed ``am1elf10`` option to the OpenEye implementation of
  :py:class:`assign_partial_charges <openforcefield.utils.toolkits.OpenEyeToolkitWrapper>` which
  calculates the average partial charges at the AM1 level of theory using conformers selected using the ELF10 method.
- `PR #643 <https://github.com/openforcefield/openforcefield/pull/643>`_: Adds
  :py:class:`openforcefield.typing.engines.smirnoff.forcefield.get_available_force_fields <openforcefield.typing.engines.smirnoff.forcefield.get_available_force_fields>`,
  which returns paths to the files of force fields available through entry point plugins.


0.7.0 - Charge Increment Model, Proper Torsion interpolation, and new Molecule methods
--------------------------------------------------------------------------------------

This is a relatively large release, motivated by the idea that changing existing functionality is bad
so we shouldn't do it too often, but when we do change things we should do it all at once.

Here's a brief rundown of what changed, migration tips, and how to find more details in the full release notes below:

* To provide more consistent partial charges for a given molecule, existing conformers are now disregarded by default
  by ``Molecule.assign_partial_charges``. Instead, new conformers are generated for use in semiempirical calculations.
  Search for ``use_conformers``.
* Formal charges are now always returned as ``simtk.unit.Quantity`` objects, with units of elementary charge.
  To convert them to integers, use ``from simtk import unit`` and
  ``atom.formal_charge.value_in_unit(unit.elementary_charge)`` or
  ``mol.total_charge.value_in_unit(unit.elementary_charge)``.
  Search ``atom.formal_charge``.
* The OpenFF Toolkit now automatically reads and writes partial charges in SDF files. Search for
  ``atom.dprop.PartialCharges``.
* The OpenFF Toolkit now has different behavior for handling multi-molecule and multi-conformer SDF files. Search
  ``multi-conformer``.
* The OpenFF Toolkit now distinguishes between partial charges that are all-zero and partial charges that are unknown.
  Search ``partial_charges = None``.
* ``Topology.to_openmm`` now assigns unique atoms names by default. Search ``ensure_unique_atom_names``.
* Molecule equality checks are now done by graph comparison instead of SMILES comparison.
  Search ``Molecule.are_isomorphic``.
* The ``ChemicalEnvironment`` module was almost entirely removed, as it is an outdated duplicate of some Chemper
  functionality. Search ``ChemicalEnvironment``.
* ``TopologyMolecule.topology_particle_start_index`` has been removed from the ``TopologyMolecule`` API, since atoms
  and virtualsites are no longer contiguous in the ``Topology`` particle indexing system. Search
  ``topology_particle_start_index``.
* ``compute_wiberg_bond_orders`` has been renamed to ``assign_fractional_bond_orders``.

There are also a number of new features, such as:

* Support for ``ChargeIncrementModel`` sections in force fields.
* Support for ``ProperTorsion`` ``k`` interpolation in force fields using fractional bond orders.
* Support for AM1-Mulliken, Gasteiger, and other charge methods using the new ``assign_partial_charges`` methods.
* Support for AM1-Wiberg bond order calculation using either the OpenEye or RDKit/AmberTools backends and the
  ``assign_fractional_bond_orders`` methods.
* Initial (limited) interoperability with QCArchive, via ``Molecule.to_qcschema`` and ``from_qcschema``.
* A ``Molecule.visualize`` method.
* Several additional ``Molecule`` methods, including state enumeration and mapped SMILES creation.

**Major Feature: Support for the SMIRNOFF ChargeIncrementModel tag**

`The ChargeIncrementModel tag in the SMIRNOFF specification <https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html#chargeincrementmodel-small-molecule-and-fragment-charges>`_
provides analagous functionality to AM1-BCC, except that instead of AM1-Mulliken charges, a number of different charge
methods can be called, and instead of a fixed library of two-atom charge corrections, an arbitrary number of
SMIRKS-based, N-atom charge corrections can be defined in the SMIRNOFF format.

The initial implementation of the SMIRNOFF ``ChargeIncrementModel`` tag accepts keywords for ``version``,
``partial_charge_method``, and ``number_of_conformers``. ``partial_charge_method`` can be any string, and it is
up to the ``ToolkitWrapper``'s ``compute_partial_charges`` methods to understand what they mean. For
geometry-independent ``partial_charge_method`` choices, ``number_of_conformers`` should be set to zero.

SMIRKS-based parameter application for ``ChargeIncrement`` parameters is different than other SMIRNOFF sections.
The initial implementation of ``ChargeIncrementModelHandler`` follows these rules:

* an atom can be subject to many ``ChargeIncrement`` parameters, which combine additively.
* a ``ChargeIncrement`` that matches a set of atoms is overwritten only if another ``ChargeIncrement``
  matches the same group of atoms, regardless of order. This overriding follows the normal SMIRNOFF hierarchy.

To give a concise example, what if a molecule ``A-B(-C)-D`` were being parametrized, and the force field
defined ``ChargeIncrement`` SMIRKS in the following order?

1) ``[A:1]-[B:2]``
2) ``[B:1]-[A:2]``
3) ``[A:1]-[B:2]-[C:3]``
4) ``[*:1]-[B:2](-[*:3])-[*:4]``
5) ``[D:1]-[B:2](-[*:3])-[*:4]``

In the case above, the ChargeIncrement from parameters 1 and 4 would NOT be applied to the molecule,
since another parameter matching the same set of atoms is specified further down in the parameter hierarchy
(despite those subsequent matches being in a different order).

Ultimately, the ChargeIncrement contributions from parameters 2, 3, and 5 would be summed and applied.

It's also important to identify a behavior that these rules were written to *avoid*: if not for the
"regardless of order" clause in the second rule, parameters 4 and 5 could actually have been applied six and two times,
respectively (due to symmetry in the SMIRKS and the use of wildcards). This situation could also arise as a result
of molecular symmetry. For example, a methyl group could match the SMIRKS ``[C:1]([H:2])([H:3])([H:4])`` six ways
(with different orderings of the three hydrogen atoms), but the user would almost certainly not intend for the charge
increments to be applied six times. The "regardless of order" clause was added specifically to address this.

In short, the first time a group of atoms becomes involved in a ``ChargeIncrement`` together, the System gains a new
parameter "slot". Only another ``ChargeIncrement`` which applies to the exact same group of atoms (in any order) can
take over the "slot", pushing the original ``ChargeIncrement`` out.

**Major Feature: Support for ProperTorsion k value interpolation**

`Chaya Stern's work <https://chayast.github.io/fragmenter-manuscript/>`_
showed that we may be able to produce higher-quality proper torsion parameters by taking into
account the "partial bond order" of the torsion's central bond. We now have the machinery to compute AM1-Wiberg
partial bond orders for entire molecules using the ``assign_fractional_bond_orders`` methods of either  ``OpenEyeToolkitWrapper`` or ``AmberToolsToolkitWrapper``. The thought is that, if some simple electron population analysis shows
that a certain aromatic bond's order is 1.53, maybe rotations about that bond can be described well by interpolating
53% of the way between the single and double bond k values.

Full details of how to define a torsion-interpolating SMIRNOFF force fields are available in
`the ProperTorsions section of the SMIRNOFF specification <https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html#fractional-torsion-bond-orders>`_.

Behavior changed
""""""""""""""""
- `PR #508 <https://github.com/openforcefield/openforcefield/pull/508>`_:
  In order to provide the same results for the same chemical species, regardless of input
  conformation,
  :py:class:`Molecule <openforcefield.topology.Molecule>`
  ``assign_partial_charges``, ``compute_partial_charges_am1bcc``, and
  ``assign_fractional_bond_orders`` methods now default to ignore input conformers
  and generate new conformer(s) of the molecule before running semiempirical calculations.
  Users can override this behavior by specifying the keyword argument
  ``use_conformers=molecule.conformers``.
- `PR #281 <https://github.com/openforcefield/openforcefield/pull/281>`_: Closes
  `Issue #250 <https://github.com/openforcefield/openforcefield/issues/250>`_
  by adding support for partial charge I/O in SDF. The partial charges are stored as a property in the
  SDF molecule block under the tag ``<atom.dprop.PartialCharge>``.
- `PR #281 <https://github.com/openforcefield/openforcefield/pull/281>`_: If a
  :py:class:`Molecule <openforcefield.topology.Molecule>`'s
  ``partial_charges`` attribute is set to ``None`` (the default value), calling ``to_openeye`` will
  now produce a OE molecule with partial charges set to ``nan``. This would previously produce an OE
  molecule with partial charges of 0.0, which was a loss of information, since it wouldn't be clear
  whether the original OFFMol's partial charges were REALLY all-zero as opposed to ``None``. OpenEye toolkit
  wrapper methods such as ``from_smiles`` and ``from_file`` now produce OFFMols with
  ``partial_charges = None`` when appropriate (previously these would produce OFFMols with
  all-zero charges, for the same reasoning as above).
- `PR #281 <https://github.com/openforcefield/openforcefield/pull/281>`_:
  :py:class:`Molecule <openforcefield.topology.Molecule>`
  ``to_rdkit``
  now sets partial charges on the RDAtom's ``PartialCharges`` property (this was previously set
  on the ``partial_charges`` property). If the
  :py:class:`Molecule <openforcefield.topology.Molecule>`'s partial_charges attribute is ``None``, this property
  will not be defined on the RDAtoms.
- `PR #281 <https://github.com/openforcefield/openforcefield/pull/281>`_:
  Enforce the behavior during SDF I/O that a SDF may contain multiple
  `molecules`, but that the OFF Toolkit
  does not assume that it contains multiple `conformers of the same molecule`. This is an
  important distinction, since otherwise there is ambiguity around whether properties of one
  entry in a SDF are shared among several molecule blocks or not, or how to resolve conflicts if properties
  are defined differently for several "conformers" of chemically-identical species (More info
  `here <https://docs.eyesopen.com/toolkits/python/oechemtk/oemol.html#dude-where-s-my-sd-data>`_).
  If the user requests the OFF Toolkit to write a multi-conformer
  :py:class:`Molecule <openforcefield.topology.Molecule>` to SDF, only the first conformer will be written.
  For more fine-grained control of writing properties, conformers, and partial charges, consider
  using ``Molecule.to_rdkit`` or ``Molecule.to_openeye`` and using the functionality offered by
  those packages.
- `PR #281 <https://github.com/openforcefield/openforcefield/pull/281>`_: Due to different
  constraints placed on the data types allowed by external toolkits, we make our best effort to
  preserve :py:class:`Molecule <openforcefield.topology.Molecule>`
  ``properties`` when converting molecules to other packages, but users should be aware that
  no guarantee of data integrity is made. The only data format for keys and values in the property dict that
  we will try to support through a roundtrip to another toolkit's Molecule object is ``string``.
- `PR #574 <https://github.com/openforcefield/openforcefield/pull/574>`_: Removed check that all
  partial charges are zero after assignment by ``quacpac`` when AM1BCC used for charge assignment.
  This check fails erroneously for cases in which the partial charge assignments are correctly all zero,
  such as for ``N#N``. It is also an unnecessary check given that ``quacpac`` will reliably indicate when
  it has failed to assign charges.
- `PR #597 <https://github.com/openforcefield/openforcefield/pull/597>`_: Energy-minimized sample systems
  with Parsley 1.1.0.
- `PR #558 <https://github.com/openforcefield/openforcefield/pull/558>`_: The
  :py:class:`Topology <openforcefield.topology.Topology>`
  particle indexing system now orders :py:class:`TopologyVirtualSites <openforcefield.topology.TopologyVirtualSite>`
  after all atoms.
- `PR #469 <https://github.com/openforcefield/openforcefield/pull/469>`_:
  When running :py:meth:`Topology.to_openmm <openforcefield.topology.Topology.to_openmm>`, unique atom names
  are generated if the provided atom names are not unique (overriding any existing atom names). This
  uniqueness extends only to atoms in the same molecule. To disable this behavior, set the kwarg
  ``ensure_unique_atom_names=False``.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_:
  :py:meth:`Molecule.__eq__ <openforcefield.topology.Molecule.__eq__>` now uses the new
  :py:meth:`Molecule.are_isomorphic <openforcefield.topology.Molecule.are_isomorphic>` to perform the
  similarity checking.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_:
  The :py:meth:`Topology.from_openmm <openforcefield.topology.Topology.from_openmm>` and
  :py:meth:`Topology.add_molecule <openforcefield.topology.Topology.add_molecule>` methods now use the
  :py:meth:`Molecule.are_isomorphic <openforcefield.topology.Molecule.are_isomorphic>` method to match
  molecules.
- `PR #551 <https://github.com/openforcefield/openforcefield/pull/551>`_: Implemented the
  :py:meth:`ParameterHandler.get_parameter <openforcefield.typing.engines.smirnoff.parameters.ParameterHandler.get_parameter>`
  function (would previously return ``None``).

API-breaking changes
""""""""""""""""""""
- `PR #471 <https://github.com/openforcefield/openforcefield/pull/471>`_: Closes
  `Issue #465 <https://github.com/openforcefield/openforcefield/issues/465>`_.
  ``atom.formal_charge`` and ``molecule.total_charge`` now return ``simtk.unit.Quantity`` objects
  instead of integers. To preserve backward compatibility, the setter for ``atom.formal_charge``
  can accept either a ``simtk.unit.Quantity`` or an integer.
- `PR #601 <https://github.com/openforcefield/openforcefield/pull/601>`_: Removes
  almost all of the previous
  :py:class:`ChemicalEnvironment <openforcefield.typing.chemistry.ChemicalEnvironment>`
  API, since this entire module was simply copied from
  `Chemper <https://github.com/MobleyLab/chemper>`_ several years ago and has fallen behind on updates.
  Currently only
  :py:meth:`ChemicalEnvironment.get_type <openforcefield.typing.chemistry.ChemicalEnvironment.get_type>`,
  :py:meth:`ChemicalEnvironment.validate <openforcefield.typing.chemistry.ChemicalEnvironment.validate>`,
  and an equivalent classmethod
  :py:meth:`ChemicalEnvironment.validate_smirks <openforcefield.typing.chemistry.ChemicalEnvironment.validate_smirks>`
  remain. Also, please comment on
  `this GitHub issue <https://github.com/MobleyLab/chemper/issues/90>`_ if you HAVE been using
  the previous extra functionality in this module and would like us to prioritize creation of a Chemper
  conda package.
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
- `PR #595 <https://github.com/openforcefield/openforcefield/pull/595>`_: Removed functions
  ``openforcefield.utils.utils.temporary_directory`` and
  ``openforcefield.utils.utils.temporary_cd`` and replaced their behavior with
  ``tempfile.TemporaryDirectory()``.

New features
""""""""""""
- `PR #471 <https://github.com/openforcefield/openforcefield/pull/471>`_: Closes
  `Issue #208 <https://github.com/openforcefield/openforcefield/issues/208>`_
  by implementing support for the
  ``ChargeIncrementModel`` tag in the `SMIRNOFF specification <https://open-forcefield-toolkit.readthedocs.io/en/latest/smirnoff.html#chargeincrementmodel-small-molecule-and-fragment-charges>`_.
- `PR #471 <https://github.com/openforcefield/openforcefield/pull/471>`_: Implements
  ``Molecule.assign_partial_charges``, which calls one of the newly-implemented
  ``OpenEyeToolkitWrapper.assign_partial_charges``, and
  ``AmberToolsToolkitWrapper.assign_partial_charges``. ``strict_n_conformers`` is a
  optional boolean keyword argument indicating whether an ``IncorrectNumConformersError`` should be raised if an invalid
  number of conformers is supplied during partial charge calculation. For example, if two conformers are
  supplied, but ``partial_charge_method="AM1BCC"`` is also set, then there is no clear use for
  the second conformer. The previous behavior in this case was to raise a warning, and to preserve that
  behavior, ``strict_n_conformers`` defaults to a value of ``False``.
- `PR #471 <https://github.com/openforcefield/openforcefield/pull/471>`_: Adds
  keyword argument ``raise_exception_types`` (default: ``[Exception]``) to
  :py:meth:`ToolkitRegistry.call <openforcefield.utils.toolkits.ToolkitRegistry.call>`.
  The default value will provide the previous OpenFF Toolkit behavior, which is that the first ToolkitWrapper
  that can provide the requested method is called, and it either returns on success or raises an exception. This new
  keyword argument allows the ToolkitRegistry to *ignore* certain exceptions, but treat others as fatal.
  If ``raise_exception_types = []``, the ToolkitRegistry will attempt to call each ToolkitWrapper that provides the
  requested method and if none succeeds, a single ``ValueError`` will be raised, with text listing the
  errors that were raised by each ToolkitWrapper.
- `PR #601 <https://github.com/openforcefield/openforcefield/pull/601>`_: Adds
  :py:meth:`RDKitToolkitWrapper.get_tagged_smarts_connectivity <openforcefield.utils.toolkits.RDKitToolkitWrapper.get_tagged_smarts_connectivity>`
  and
  :py:meth:`OpenEyeToolkitWrapper.get_tagged_smarts_connectivity <openforcefield.utils.toolkits.OpenEyeToolkitWrapper.get_tagged_smarts_connectivity>`,
  which allow the use of either toolkit for smirks/tagged smarts validation.
- `PR #600 <https://github.com/openforcefield/openforcefield/pull/600>`_:
  Adds :py:meth:`ForceField.__getitem__ <openforcefield.typing.engines.smirnoff.forcefield.ForceField.__getitem__>`
  to look up ``ParameterHandler`` objects based on their string names.
- `PR #508 <https://github.com/openforcefield/openforcefield/pull/508>`_:
  Adds :py:meth:`AmberToolsToolkitWrapper.assign_fractional_bond_orders <openforcefield.utils.toolkits.AmberToolsToolkitWrapper.assign_wiberg_bond_orders>`.
- `PR #469 <https://github.com/openforcefield/openforcefield/pull/469>`_: The
  :py:class:`Molecule <openforcefield.topology.Molecule>` class adds
  :py:meth:`Molecule.has_unique_atom_names <openforcefield.topology.Molecule.has_unique_atom_names>`
  and :py:meth:`Molecule.has_unique_atom_names <openforcefield.topology.Molecule.generate_unique_atom_names>`.
- `PR #472 <https://github.com/openforcefield/openforcefield/pull/472>`_:
  Adds to the :py:class:`Molecule <openforcefield.topology.Molecule>` class
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
      .. note::
         The to_qcschema method accepts an extras dictionary which is passed into the validated qcelemental.models.Molecule
         object.
- `PR #506 <https://github.com/openforcefield/openforcefield/pull/506>`_:
  The :py:class:`Molecule <openforcefield.topology.Molecule>` class adds
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
  Adds ``test_forcefields/ion_charges.offxml``, giving ``LibraryCharges`` for monatomic ions.
- `PR #543 <https://github.com/openforcefield/openforcefield/pull/543>`_:
  Adds 3 new methods to the :py:class:`Molecule <openforcefield.topology.Molecule>` class which allow the enumeration of molecule
  states. These are :py:meth:`Molecule.enumerate_tautomers <openforcefield.topology.Molecule.enumerate_tautomers>`,
  :py:meth:`Molecule.enumerate_stereoisomers <openforcefield.topology.Molecule.enumerate_stereoisomers>`,
  :py:meth:`Molecule.enumerate_protomers <openforcefield.topology.Molecule.enumerate_protomers>`
      .. warning::
         Enumerate protomers is currently only available through the OpenEye toolkit.
- `PR #573 <https://github.com/openforcefield/openforcefield/pull/573>`_:
  Adds ``quacpac`` error output to ``quacpac`` failure in ``Molecule.compute_partial_charges_am1bcc``.
- `PR #560 <https://github.com/openforcefield/openforcefield/issues/560>`_: Added visualization method to the the Molecule class.
- `PR #620 <https://github.com/openforcefield/openforcefield/pull/620>`_: Added the ability to register parameter handlers via entry point plugins. This functionality is accessible by initializing a ``ForceField`` with the ``load_plugins=True`` keyword argument. 
- `PR #582 <https://github.com/openforcefield/openforcefield/pull/582>`_: Added fractional bond order interpolation
  Adds `return_topology` kwarg to
  :py:meth:`Forcefield.create_openmm_system <openforcefield.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system>`,
  which returns the processed topology along with the system when ``True`` (default ``False``).

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
  :py:meth:`Forcefield.create_openmm_system <openforcefield.typing.engines.smirnoff.forcefield.ForceField.create_openmm_system>`
  to check that a protein system can be created.
- `PR #506 <https://github.com/openforcefield/openforcefield/pull/506>`_: Added a test for the molecule
  identified in issue #513 as losing aromaticity when converted to rdkit.
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
- `PR #582 <https://github.com/openforcefield/openforcefield/pull/582>`_: Added fractional bond order interpolation tests,
  tests for :py:class:`ValidatedDict <openforcefield.utils.collections.ValidatedDict>`.


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
  :py:meth:`Molecule.to_file <openforcefield.topology.Molecule.to_file>` method can now correctly write to ``MOL``
  files, in line with the supported file type list.
- `Issue #568 <https://github.com/openforcefield/openforcefield/issues/568>`_: The
  :py:meth:`Molecule.to_file <openforcefield.topology.Molecule.to_file>` can now correctly write multi-model PDB files
  when using the RDKit backend toolkit.


Examples added
""""""""""""""
- `PR #591 <https://github.com/openforcefield/openforcefield/pull/591>`_ and
  `PR #533 <https://github.com/openforcefield/openforcefield/pull/533>`_: Adds an
  `example notebook and utility to compute conformer energies <https://github.com/openforcefield/openforcefield/blob/master/examples/conformer_energies>`_.
  This example is made to be reverse-compatible with the 0.6.0 OpenFF Toolkit release.
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
