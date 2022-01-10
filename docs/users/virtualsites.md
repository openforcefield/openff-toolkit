(virtualsites)=

# Virtual sites

The Open Force Field Toolkit supports the SMIRNOFF virtual site specification
for models using off-site charges, including 4- and 5-point water models, in
addition to lone pair modelling on various functional groups. The primary focus
is on the ability to parameterize a system using virtual sites, and generate an
OpenMM system with all virtual sites present and ready for evaluation. Support
for formats other than OpenMM has not yet been implemented, but may come with
the appearance of the OpenFF system object. In addition to implementing the
specification, the toolkit [`Molecule`](openff.toolkit.topology.Molecule) object allows the creation and manipulation of virtual sites.

## Support for the SMIRNOFF VirtualSite tag

Virtual sites can be added to a System in two ways:

* SMIRNOFF Force Fields can contain a [VirtualSites tag](https://openforcefield.github.io/standards/standards/smirnoff/#virtualsites-virtual-sites-for-off-atom-charges),
  specifying the addition of virtual sites according to SMARTS-based rules.
* Virtual sites can be added to a [`Molecule`](openff.toolkit.topology.Molecule), and
  these will appear in the final OpenMM system if a virtual site handler is present
  in the[`ForceField`](openff.toolkit.typing.engines.smirnoff.forcefield.ForceField).

Virtual sites directly depend on 3D conformation, because the position of the
virtual sites depends on vectors defined by the atoms specified during
parameterization. Because of this, a virtual site matching the triplet of atoms
1-2-3 will define a point that is different from a triplet matching 3-2-1. This
is similar to defining "right-handed" and "left-handed" coordinate systems.
This subtlety plays into two major concepts in force field development:

1. We sometimes want to define a single virtual site describing two points with the
   same parameters (distance, angle, etc.) via symmetry, such as 5-point water
   models.
2. We have a match that produces multiple orderings of the atoms (e.g. if wildcards
   are present in the SMARTS pattern), and we only want one to be applied.

Case (1) is very useful for parameter optimization, where a single SMARTS-based
parameter can be used to optimize both points, such as the angle defining the
virtual points for a 5-point water model. Case (2) is the typical scenario for
the nitrogen lone pair in ammonia, where only one point needs to be specified.
We discuss a few more illustrative examples below. Beyond these attributes, the
virtual site specification allows a policy for specifying how to handle
exclusions in the OpenMM force evaluator. The current default is to add
pairwise energy exclusions in the OpenMM system between a virtual site and all
tagged atoms matched in its SMARTS(`exclusion_policy="parents",` ). Currently
defined are `"none"`, `"minimal"`, and `"parents"`, where `"minimal"` specifies
the single atom that the virtual site defines as the "origin". For water, for
example, `"minimal"` would mean just the oxygen, whereas `"parents"` would mean
all three atoms.

In order to give consistent and intended behavior, the specification was
modified from its draft form in following manner: The `"name"` and `"match"`
attributes have been added to each virtual site parameter type. These changes
allow for

* specifying different virtual site types using the same atoms
* allowing two virtual sites with the same type and same atoms but different
  physical parameters to be added simultaneously
* allowing the ability to control whether the virtual site encodes one or
  multiple particles, based on the number of ways the matching atoms can be
  ordered.

The `"name"` attribute encodes whether the virtual site to be added should
override an existing virtual site of the same type (e.g. hierarchy preference),
or if this virtual site should be added in addition to the other existing
virtual sites on the given atoms. This means that different virtual site types
can share the same group of parent atoms and use the same name without
overwriting each other (the default `name` is `EP` for all sites, which gives
the expected hierarchical behavior used in other SMIRNOFF tags).

The `"match"` attribute accepts either `"once"` or `"all_permutations"`,
offering control for situations where a SMARTS pattern can possibly match the
same group of atoms in different orders(either due to wildcards or local
symmetry) and it is desired to either add just one or all of the possible
virtual particles. The default value is `"all_permutations",` but for
`TrivalentLonePair` it is always set to `"once"`, regardless of what the file
contains, since all orderings always place the particle in the exact same
position.

The following cases exemplify our reasoning in implementing this behavior, and
should draw caution to complex issues that may arise when designing virtual
site parameters. Let us consider 4-, 5-, and 6-point water models:

* A 4-point water model with a `DivalentLonePair`: This can be implemented by
  specifying `match="once"`, `outOfPlaneAngle="0*degree"`, and
  `distance=-.15*angstrom"`. Since the SMIRKS pattern `"[#1:1]-[#8X2:2]-
  [#2:3]"` would match water twice and would create two particles in the exact
  same position if `all_permutations` was specified, we specify `"once"` to
  have only one particle generated. Although having two particles in the same
  position should not affect the physics if the proper exclusion policy is
  applied, it would effectively make the 4-point model just as expensive as
  5-point models.

* A 5-point water model with a `DivalentLonePair`: This can be implemented by
  using `match="all_permutations"` (unlike the 4-point model),
  `outOfPlaneAngle="56.26*degree`, and `distance=0.7*angstrom`, for example.
  Here the permutations will cause particles to be placed at ±56.26 degrees,
  and changing any of the physical quantities will affect *both* particles.

* A 6-point water model with both `DivalentLonePair` sites above. Since these
  two parameters look identical, it is unclear whether they should both be
  applied or if one should override the other. The toolkit never compares the
  physical numbers to determine equality as this can lead to instability during
  e.g. parameter fitting. To get this to work, we specify `name="EP1"` for the
  first parameter, and `name="EP2"` for the second parameter. This instructs
  the parameter handler keep them separate, and therefore both are applied.
  (If both had the same name, then the typical SMIRNOFF hierarchy rules are
  used, and only the last matched parameter would be applied.)

* Dinitrogen, `N#N` with a `BondCharge` virtual site. Since we want a
  `BondCharge` on both ends, we specify `match="all_permutations"`.

* Formaldehyde, `H2C=O`, with `MonovalentLonePair` virtual site(s) on the
  oxygen, with the aim of modeling both lone pairs. This one is subtle, since `
  [#1:3]-[#6X3:2]=[#8X1:1]` matches two unique groups of atoms (`1-3-4` and
  `2-3-4`). It is important to note in this situation that
  `match="all_permutations"` behaves exactly the same as `match="once"`. Due to
  the anchoring hydrogens (`1` and `2`) being symmetric but opposite about the
  bond between `3` and `4`, a single parameter does correctly place both lone
  pairs. A standing issue here is that the default exclusion policy
  (`parents`) will allow these two virtual sites to interact since they have
  different indexed atoms (parents), causing the energy to be different than
  the non-virtual site parameterization. In the future, the
  `exclusion_policy="local"` will account for this, and make virtual sites that
  share at least one "parent" atom not interact with each other. As a special
  note: when applying a `MonovalentLonePair` to a completely symmetric
  molecule, e.g. water, `all_permutations` can come into play, but this will
  apply two particles (one for each hydrogen).

Finally, the toolkit handles the organization of atoms and virtual sites in a
specific manner. Virtual sites are expected to be added *after all molecules in
the topology are present*. This is because the Open Force Field Toolkit
organizes a topology by placing all atoms first, then all virtual sites last.
This differs from the OpenMM Modeller object, for example, which interleaves
the order of atoms and virtual sites in such a way that all particles of a
molecule are contiguous. In addition, due to the fact that a virtual site may
contain multiple particles coupled to single parameters, the toolkit makes a
distinction between a virtual *site*, and a virtual *particle*. A virtual site
may represent multiple virtual particles, so the total number of particles
cannot be directly determined by simply summing the number of atoms and virtual
sites in a molecule. This is taken into account, however, and the 
[`Molecule`](openff.toolkit.topology.Molecule) and
[`Topology`](openff.toolkit.topology.Topology) classes both implement `particle`
iterators.

