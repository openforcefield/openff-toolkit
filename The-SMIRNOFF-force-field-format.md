# The SMIRks Native Open Force Field (SMIRNOFF) v0.2

The SMIRNOFF format is based on the [`OpenMM`](http://openmm.org) [`ForceField`](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField) class and provides an XML format for encoding force fields based on [SMIRKS](http://www.daylight.com/dayhtml/doc/theory/theory.smirks.html)-based chemical perception.
While designed for [`OpenMM`](http://openmm.org), parameters encoded in this format can be applied to systems and then these systems converted via [`ParmEd`](http://parmed.github.io/ParmEd) and [`InterMol`](https://github.com/shirtsgroup/InterMol) for simulations in a variety of other simulation packages.

## Basic structure

The SMIRNOFF format provides XML `ffxml` files that are parseable by the `ForceField` class of the `openforcefield.typing.smirnoff` module.
These encode parameters for a force field based on a SMIRKS-based specification of the chemical environment the parameters are to be applied to.
The file has tags corresponding to OpenMM force terms (`BondForce`, `AngleForce`, `TorsionForce`, etc., as discussed in more detail below); these specify units used for the different constants provided for individual force terms, for example (see the [AlkEthOH example ffxml](https://github.com/openforcefield/openforcefield/blob/master/openforcefield/data/forcefield/Frosst_AlkEthOH.ffxml)):
```XML
<AngleForce angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
   ...
</AngleForce>     
```
which introduces following `Angle` terms which will use units of degrees for the angle and kilocalories per mole per square radian for the force constant.

Under each of these force terms, there are tags for individual parameter lines such as these:
```XML
<AngleForce angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
   <Angle smirks="[a,A:1]-[#6X4:2]-[a,A:3]" angle="109.50" k="100.0"/>
   <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.50" k="70.0"/>
</AngleForce>     
```
The first of these specifies the `[a,A:1]-[#6X4:2]-[a,A:3]` SMIRKS pattern for an angle, with a tetravalent carbon at the center with single bonds to two atoms of any type.
Atoms are labeled 1, 2, and 3, with 2 being the central atom.
Equilibrium angle values are provided, along with force constants (with units as given above).

**SMIRNOFF parameters are hierarchical** in that parameters which come later in a file override those which come earlier if they match the same pattern.
 This can be seen in the example above, where the first line provides a generic angle parameter for any tetravalent carbon (single bond) angle, and the second line overrides this for the specific case of a hydrogen-(tetravalent carbon)-hydrogen angle.

This hierarchical structure means that a typical parameter file will tend to have generic parameters early in the section for each force type, with more specialized parameters assigned later.

**Technical note**: Because this is an XML format, certain special characters that occur in valid SMIRKS patterns (such as ampersand symbols `&`) must be treated specially to process.
See [this entry](https://en.wikipedia.org/wiki/List_of_XML_and_HTML_character_entity_references) for more details.

## Units

If not specified by a `*_unit` attribute, the [standard OpenMM unit system](http://docs.openmm.org/latest/userguide/theory.html#units) is assumed:
* distances are specified in `nanometers`
* mass is specified in `daltons`
* time is specified in `picoseconds`
* charge is specified in `elementary_charge`,
* temperature is specified in `kelvin`
* quantity is specified in `moles`
* angles are specified in `radians`

Allowed values for units are given in [`simtk.unit`](https://github.com/pandegroup/openmm/blob/master/wrappers/python/simtk/unit/unit_definitions.py).

## Functional forms, etc.

**Functional form**: The SMIRNOFF format specifies parameters; once specified, these are processed by the SMIRNOFF `ForceField` class and used to assign parameters to OpenMM Forces.
This means that specific forces are generally implemented as discussed in the [OpenMM Documentation](http://docs.openmm.org/latest/userguide/theory.html), see especially [Section 19 on Standard Forces](http://docs.openmm.org/latest/userguide/theory.html#standard-forces) for functional forms.
In some cases, typically for consistency with the AMBER force field philosophy motivating some of the authors, we do some manipulation of parameters from these files as discussed below in "Parameter sections".

**Charges**: In keeping with the AMBER force field philosophy, especially as implemented in small molecule force fields such as [GAFF](http://ambermd.org/antechamber/gaff.html), [GAFF2](https://mulan.swmed.edu/group/gaff.php), and [parm@Frosst](http://www.ccl.net/cca/data/parm_at_Frosst/), we at least initially treat partial charges as something to be obtained separately from the rest of the force field (BAIT, which includes bonds, angles, impropers, and torsions; and vdW terms), typically via QM calculations or a method such as the [AM1-BCC](https://dx.doi.org/10.1002/jcc.10128) approach from Christopher Bayly, thus, for system setup we provide the option of specifying a charging method, though charges are not normally specified in the FFXML itself.
For other force fields with "library"-style charges, we could introduce a new section providing specific charges via SMIRKS pattern, though this is not yet implemented. Bond charge corrections, however, are supported as we discuss below.

## Parameter sections

For this section it will help to have on hand an example SMIRNOFF file, such as that the [AlkEthOH example ffxml](https://github.com/openforcefield/openforcefield/blob/master/openforcefield/data/forcefield/Frosst_AlkEthOH.ffxml) or the larger prototype [SMIRNOFF99Frosst ffxml](https://github.com/openforcefield/SMIRNOFF99Frosst/blob/master/SMIRNOFF99Frosst.ffxml).

Before getting in to individual sections, it's worth noting that the XML parser ignores attributes in the XML that it does not understand, so providing a parameter line for an angle that specifies (for example) a second force constant `k2` will lead to no effect.

### `<StericsForce>`

Nonbonded sterics parameters, which model repulsive and dispersive forces, are specified via the `StericsForce` tag with sub-tags for individual `Atom` entries, such as:
```XML
<StericsForce potential="Lennard-Jones" scale12="0.0" scale13="0.0" scale14="0.5" sigma_unit="angstroms" epsilon_unit="kilocalories_per_mole">
   <Atom smirks="[#1:1]" rmin_half="1.4870" epsilon="0.0157"/>
   <Atom smirks="[#1:1]-[#6]" rmin_half="1.4870" epsilon="0.0157"/>
   ...
</StericsForce>
```
Attributes in the `<StericsForce>` tag specify the scaling terms applied to the energies of 1-2 (`scale12`, default: 0), 1-3 (`scale13`, default: 0), and 1-4 (`scale14`, default: 0.5) interactions.

Currently, only `potential="Lennard-Jones"` is supported:
```
U(r) = 4*epsilon*((sigma/r)^12 - (sigma/r)^6)
```
Later revisions will add support for additional potential types (e.g., Buckingham exp-6) and the ability to support arbitrary algebraic functional forms.
If the `potential` attribute is omitted, it defaults to `Lennard-Jones`.

For compatibility, the size property of an atom can be specified either via providing the `sigma` attribute, such as `sigma="1.3"`, or via the `r_0/2` (`rmin/2`) values used in AMBER force fields (here denoted `rmin_half` as in the example above).
The two are related by `r0 = 2^(1/6)*sigma` and conversion is done internally in `ForceField` into the `sigma` values used in OpenMM.
`epsilon` denotes the well depth.

### `<ElectrostaticsForce>`

Electrostatic interactions are specified via the `<ElectrostaticsForce>` tag.
```XML
<ElectrostaticsForce method="PME" scale12="0.0" scale13="0.0" scale14="0.833333"/>
```
The `method` attribute specifies the manner in which electrostatic interactions are to be computed:
* `PME` - [particle mesh Ewald](http://docs.openmm.org/latest/userguide/theory.html#coulomb-interaction-with-particle-mesh-ewald) should be used (DEFAULT)
* `reaction-field` - [reaction-field electrostatics](http://docs.openmm.org/latest/userguide/theory.html#coulomb-interaction-with-cutoff) should be used
* `Coulomb` - direct Coulomb interactions (with no reaction-field attenuation) should be used

The interaction scaling parameters applied to atoms connected by a few bonds are
* `scale12` (default: 0) specifies the scaling applied to 1-2 bonds
* `scale13` (default: 0) specifies the scaling applied to 1-3 bonds
* `scale14` (default: 0.833333) specifies the scaling applied to 1-4 bonds

Currently, no child tags are used because the charge model is specified via different means (currently library charges or BCCs).

**QUESTION:** Should we also also require the `cutoff` (or lack thereof) to be specified?

**QUESTION:** Should we also also require the `switch` (or lack thereof) to be specified, since this can also have a large effect on the forces near the cutoff?

### `<BondForce>`

Bond parameters are specified via a `<BondForce>...</BondForce>` block, with individual `<Bond/>` tags containing attributes specifying the equilibrium bond length (`length`) and force constant (`k`) values for specific bonds.
For example:
```XML
<BondForce potential="harmonic" length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
   <Bond smirks="[#6X4:1]-[#6X4:2]" length="1.526" k="620.0"/>
   <Bond smirks="[#6X4:1]-[#1:2]" length="1.090" k="680.0"/>
   ...
</BondForce>
```
Currently, only `potential="harmonic"` is supported, where we utilize the standard harmonic functional form:
```
U(r) = (k/2)*(r-length)^2
```
Later revisions will add support for additional potential types and the ability to support arbitrary algebraic functional forms.
If the `potential` attribute is omitted, it defaults to `harmonic`.

**Note that AMBER and CHARMM define a modified functional form**, such that `U(r) = k*(r-length)^2`, so that force constants would need to be multiplied by two in order to be used in the SMIRNOFF format.

Constrained bonds are handled by a separate `<Constraints>` tag, which can either specify constraint distances or draw them from equilibrium distances specified in `<BondForce>`.

### `<AngleForce>`

Angle parameters are specified via an `<AngleForce>...</AngleForce>` block, with individual `<Angle/>` tags containing attributes specifying the equilibrium angle (`angle`) and force constant (`k`), as in this example:
```XML
<AngleForce potential="harmonic" angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
   <Angle smirks="[a,A:1]-[#6X4:2]-[a,A:3]" angle="109.50" k="100.0"/>
   <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.50" k="70.0"/>
   ...
</AngleForce>
```
Currently, only `potential="harmonic"` is supported, where we utilize the standard harmonic functional form:
```
U(r) = (k/2)*(theta-angle)^2
```
Later revisions will add support for additional potential types and the ability to support arbitrary algebraic functional forms.
If the `potential` attribute is omitted, it defaults to `harmonic`.

**Note that AMBER and CHARMM define a modified functional form**, such that `U(r) = k*(theta-angle)^2`, so that force constants would need to be multiplied by two in order to be used in the SMIRNOFF format.

### `<ProperTorsionForce>``

Proper torsions are specified via a `<ProperTorsionForce>...</ProperTorsionForce>` block, with individual `<Proper/>` tags containing attributes specifying the periodicity (`periodicityN`), phase (`phaseN`), barrier height (`kN`), and optional torsion multiplicity to divide the barrier height by (`idivfN`):
```XML
<ProperTorsionForce potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
   <Proper smirks="[a,A:1]-[#6X4:2]-[#6X4:3]-[a,A:4]" idivf1="9" periodicity1="3" phase1="0.0" k1="1.40"/>
   <Proper smirks="[#6X4:1]-[#6X4:2]-[#8X2:3]-[#6X4:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="0.383" idivf2="1" periodicity2="2" phase2="180.0" k2="0.1"/>
   ...
</ProperTorsionForce>
```
Here, child `Proper` tags specify at least `k1`, `phase1`, and `periodicity1` attributes for the corresponding parameters of the first force term applied to this torsion.
However, additional values are allowed in the form `kN`, `phaseN`, and `periodicityN`, where all `N` values must be consecutive (e.g., it is impermissible to specify `k1` and `k3` values without a `k2` value) but `N` can go as high as necessary.

Currently, only `potential="charmm"` is supported, where we utilize the functional form:
```
U = \sum_{i=1}^N k_i * (1 + cos(periodicity_i * phi - phase_i))
```
**Note that AMBER defines a modified functional form**, such that `U = \sum_{i=1}^N (k_i/2) * (1 + cos(periodicity_i * phi - phase_i))`, so that barrier heights would need to be divided by two in order to be used in the SMIRNOFF format.
If the `potential` attribute is omitted, it defaults to `charmm`.

Optionally, an `idivfN` attribute may be specified for each torsional term (for easier compatibility with AMBER files); this specifies a numerical value (in AMBER, always an integer) which is used as a divisor for the barrier height when assigning the torsion; i.e., a torsion with `idivf1="9"` is assigned a barrier height `k1` that is 1/9th the specified value.
If `idivfN` is not specified, the barrier height is applied as stated.

In the future, we may switch to a model where torsional barriers are [automatically divided by the number of torsions along a bond](https://github.com/openforcefield/smarty/issues/131), effectively resulting in the torsional barrier being the average of all barriers applied to that bond, rather than the current model where barriers are summed.
(Barrier heights would need to be increased accordingly.)
This would result in better handling of some cases where a small change in a molecule (such as a change in tautomer) could currently (as in AMBER) result in a dramatically different applied barrier height because of a change in the number of torsions passing through that bond.
The averaging approach would make it easier to avoid this problem without requiring as many different torsional terms.

### `<ImproperTorsionForce>`

Improper torsions are specified via an `<ImproperTorsionForce>...</ImproperTorsionForce>` block, with individual `<Improper/>` tags containing attributes that specify the same properties as `<ProperTorsionForce>`:
```XML
<ImproperTorsionForce potential="charmm" phase_unit="degrees" k_unit="kilocalories_per_mole">
    <Improper smirks="[*:1]~[#6X3:2](=[#7X2,#7X3+1:3])~[#7:4]" k1="10.5" periodicity1="2" phase1="180."/>
    ...
</ImproperTorsionForce>
```
Currently, only `potential="charmm"` is supported, where we utilize the functional form:
```
U = \sum_{i=1}^N k_i * (1 + cos(periodicity_i * phi - phase_i))
```
**Note that AMBER defines a modified functional form**, such that `U = \sum_{i=1}^N (k_i/2) * (1 + cos(periodicity_i * phi - phase_i))`, so that barrier heights would need to be divided by two in order to be used in the SMIRNOFF format.
If the `potential` attribute is omitted, it defaults to `charmm`.

**Improper torsion definitions deviate profoundly from AMBER handling of impropers** in two ways.
First, to eliminate ambiguity, we treat impropers as a [trefoil](https://upload.wikimedia.org/wikipedia/commons/thumb/5/5c/Trefoil_knot_left.svg/2000px-Trefoil_knot_left.svg.png) and apply the same set of parameters to all six paths around the trefoil.
*Because of this, all barrier heights are divided by six before we apply them*, for consistency with AMBER force fields.
Second, the *second* atom in an improper (in the example above, the trivalent carbon) is the central atom in the trefoil.

### `<GBSAForce>`

Generalized-Born surface area (GBSA) implicit solvent parameters are specified via a `<GBSAForce>...</GBSAForce>` using `<Atom/>` tags with GBSA model specific attributes:
```XML
 <GBSAForce gb_model="OBC1" solvent_dielectric="78.5" solute_dielectric="1" radius_units="nanometers" sa_model="ACE" surface_area_penalty="5.4*calories/mole/angstroms**2" solvent_radius="1.4*angstroms">
   <Atom smirks="[#1:1]" radius="0.12" scale="0.85"/>
   <Atom smirks="[#1:1]~[#6]" radius="0.13" scale="0.85"/>
   <Atom smirks="[#1:1]~[#8]" radius="0.08" scale="0.85"/>
   <Atom smirks="[#1:1]~[#16]" radius="0.08" scale="0.85"/>
   <Atom smirks="[#6:1]" radius="0.22" scale="0.72"/>
   <Atom smirks="[#7:1]" radius="0.155" scale="0.79"/>
   <Atom smirks="[#8:1]" radius="0.15" scale="0.85"/>
   <Atom smirks="[#9:1]" radius="0.15" scale="0.88"/>
   <Atom smirks="[#14:1]" radius="0.21" scale="0.8"/>
   <Atom smirks="[#15:1]" radius="0.185" scale="0.86"/>
   <Atom smirks="[#16:1]" radius="0.18" scale="0.96"/>
   <Atom smirks="[#17:1]" radius="0.17" scale="0.8"/>
 </GBSAForce>
```

#### Generalized Born (GB) model
In the `<GBSAForce/>` tag, `gb_model` selects which GB model is used.
Currently, this can be selected from a subset of the [GBSA models available in OpenMM's `simtk.openmm.app`](http://docs.openmm.org/latest/userguide/application.html#amber-implicit-solvent:
* `HCT`: [Hawkins-Cramer-Truhlar](http://docs.openmm.org/latest/userguide/zbibliography.html#hawkins1995) (corresponding to `igb=1` in AMBER): requires `[radius, scale]`
* `OBC1`: [Onufriev-Bashford-Case](http://docs.openmm.org/latest/userguide/zbibliography.html#onufriev2004) using the GB(OBC)I parameters (corresponding to `igb=2` in AMBER): requires `[radius, scale]`
* `OBC2`: [Onufriev-Bashford-Case](http://docs.openmm.org/latest/userguide/zbibliography.html#onufriev2004) using the GB(OBC)II parameters (corresponding to `igb=5` in AMBER): requires `[radius, scale]`

If the `gb_model` attribute is omitted, it defaults to `OBC1`.

Each GB model can possess several attributes, which may be unitless (`1.0`, `78.5`) or unit-bearing quantities (`1.4*angstrom`, `5.4*calories/mole/angstrom**2`).
The attributes `solvent_dielectric` and `solute_dielectric` specify solvent and solute dielectric constants used by the GB model.
In this example, `radius` and `scale` are per-particle parameters of the `OBC1` GB model supported by OpenMM.
Units are for these per-particle parameters (such as `radius_units`) optionally specified in the `<GBSAForce/>` tag.

#### Surface area (SA) penalty model

The `sa_model` attribute specifies the solvent-accessible surface area model ("SA" part of GBSA) if one should be included; if omitted, no SA term is included.

Currently, only the [analytical continuum electrostatics (ACE) model](http://docs.openmm.org/latest/userguide/theory.html#surface-area-term), designated `ACE`, can be specified, but there are plans to add more models in the future, such as the Gaussian solvation energy component of [EEF1](https://www.ncbi.nlm.nih.gov/pubmed/10223287).
If `sa_model` is not specified, it defaults to `ACE`.

The `ACE` model permits two additional parameters to be specified:
* The `surface_area_penalty` attribute specifies the surface area penalty for the `ACE` model. (Default: `5.4*calories/mole/angstroms**2`)
* The `solvent_radius` attribute specifies the solvent radius. (Default: `1.4*angstroms`)

### `<BondChargeCorrections>`

Bond charge corrections similar to those used in the [AM1-BCC](https://dx.doi.org/10.1002/jcc.10128) charge model from Christopher Bayly can be applied via a `<BondChargeCorrections>...</BondChargeCorrections>` section with child tags specifying specific `<BondChargeCorrection/>` terms.
Here is an example not intended for actual use:
```XML
<BondChargeCorrections method="AM1/CM2" increment_unit="elementary_charge">
  <BondChargeCorrection smirks="[#6X4:1]-[#6X3a:2]" increment="+0.0073"/>
  <BondChargeCorrection smirks="[#6X4:1]-[#6X3a:2]-[#7]" increment="-0.0943"/>
  <BondChargeCorrection smirks="[#6X4:1]-[#8:2]" increment="+0.0718"/>
</BondChargeCorrections>
```
The `method="AM1/CM2"` attribute specifies that a CM2 population analysis should be applied to an AM1 wavefunction to compute initial partial charges, while `elementary_charge` specifies the units used for the `increment` attribute in the `<BondChargeCorrection/>` tags.
The charge correction `increment` (in units of proton charge) will be applied on top of this by subtracting `increment` from the atom tagged as 1 and adding it to the atom tagged as 2.

### `<Constraints>`

Distance constraints can be specified through a `<Constraints/>` block, which can constrain bonds to their equilibrium lengths or specify an interatomic constraint distance.
Two atoms must be tagged in the `smirks` attribute of each `<Constraint/>` record.

To constrain two atoms to their equilibrium bond length, it is critical that a `<BondForce/>` record be specified for those atoms:
```XML
<Constraints>
  <!-- constrain all bonds to hydrogen to their equilibrium bond length -->
  <Constraint smirks="[#1:1]-[*:2]" />
</Constraints>
```
Note that the two atoms must be bonded in the specified `Topology` for the equilibrium bond length to be used.

To specify the constraint distance, or constrain two atoms that are not directly bonded (such as the hydrogens in rigid water models), specify the `distance` attribute (and optional `distance_unit` attribute for the `<Constraints/>` tag):
```XML
<Constraints distance_unit="angstroms">
  <!-- constrain water O-H bond to equilibrium bond length (overrides earlier constraint) -->
  <Constraint smirks="[#1:1]-[#8X2H2:2]-[#1]" distance="0.9572"/>
  <!-- constrain water H...H, calculating equilibrium length from H-O-H equilibrium angle and H-O equilibrium bond lengths -->
  <Constraint smirks="[#1:1]-[#8X2H2]-[#1:2]" distance="1.8532"/>
</Constraints>
```
Typical molecular simulation practice is to constrain all bonds to hydrogen to their equilibrium bond lengths and enforce rigid TIP3P geometry on water molecules:
```XML
<Constraints distance_unit="angstroms">
  <!-- constrain all bonds to hydrogen to their equilibrium bond length -->
  <Constraint smirks="[#1:1]-[*:2]" />
  <!-- TIP3P rigid water -->
  <Constraint smirks="[#1:1]-[#8X2H2:2]-[#1]" distance="0.9572"/>
  <Constraint smirks="[#1:1]-[#8X2H2]-[#1:2]" distance="1.8532"/>
</Constraints>
```

## Advanced features

Standard usage is expected to rely primarily on the features documented above and potentially new features. However, some advanced features are also currently supported.

### Versioning

The SMIRNOFF forcefield format supports versioning via the `version` attribute to the root `<SMIRNOFF>` tag, e.g.:
```XML
<SMIRNOFF version="0.2">
...
</SMIRNOFF>
```
The version format is `x.y`, where `x` denotes the major version and `y` denotes the minor version.
SMIRNOFF versions are guaranteed to be backward-compatible within the *same major version number series*, but it is possible major version increments will break backwards-compatibility.

### Partial bond orders

Partial bond orders can be used to allow interpolation of parameters. For example, these parameters:
```XML
<BondForce potential="harmonic" length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
    <Bond smirks="[#6X3:1]-[#6X3:2]" k="820.0" length="1.45"/>
    <Bond smirks="[#6X3:1]:[#6X3:2]" k="938.0" length="1.40"/>
    <Bond smirks="[#6X3:1]=[#6X3:2]" k="1098.0" length="1.35"/>
    ...
```
can be replaced by a single parameter line by first invoking the `fractional_bondorder_method` attribute to specify a method for computing the fractional bond order and `fractional_bondorder_interpolation` for specifying the procedure for interpolating parameters between specified integral bond orders:
```XML
<BondForce potential="harmonic" length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2" fractional_bondorder_method="Wiberg" fractional_bondorder_interpolation="linear">
    <Bond smirks="[#6X3:1]!#[#6X3:2]" k_bondorder1="820.0" k_bondorder2="1098" length_bondorder1="1.45" length_bondorder2="1.35"/>
    ...
```
This allows specification of force constants and lengths for bond orders 1 and 2, and then interpolation between those based on the partial bond order.
* `fractional_bondorder_method` defaults to `none`, but the `Wiberg` method is supported.
* `fractional_bondorder_interpolation` defaults to `linear`, which is the only supported scheme for now.

### Aromaticity models

Before conduct SMIRKS substructure searches, molecules are prepared using one of the supported aromaticity models, with the default model (`MDL`) used unless otherwise requested.
Alternate aromaticity models can be requested by the force field, such as
```
<SMIRNOFF version="0.2" aromaticity_model="MDL">
```
used by SMIRNOFF99Frosst (a choice by Christopher Bayly to simplify handling of certain heteroaromatic compounds).

The only aromaticity model currently widely supported (by both the [OpenEye toolkit](https://docs.eyesopen.com/toolkits/python/oechemtk/aromaticity.html) and [RDKit](http://www.rdkit.org/docs/RDKit_Book.html)) is the `MDL` model.

### Future advanced features

At present, the SMIRNOFF format basically defaults to AMBER- or OpenMM-style decisions on many issues.
For example, AMBER-style (Lorentz-Berthelot) combining rules are used, and the AMBER force field functional form.
Angles potentials are assumed to be harmonic.
However, we have plans to support other combination rules, functional forms, and angle potentials.
In keeping with the above, whole-force field decisions will be handled as attributes of the SMIRNOFF tag.
For example, alternate combination rules or functional forms might be handled as follows:
* Geometric mean combining rules: `<StericsForce potential="Lennard-Jones" combining_rule="geometric_mean">`
* A Halgren buffered 14-7 potential for vdW could be handled as `<CustomForce NonbondedForm="Halgren-buffered-14-7">`
* Selection of the type of angle force applied would be handled in a similar manner, via an `<AngleForce potential="harmonic-cos">` tag or similar. [This feature is being planned.](https://github.com/openforcefield/smarty/issues/179)
* The ability to specify arbitrary functional forms, such as `<StericsForce potential="A*exp(-B*r)+C*r^(-6)>"`

### Additional plans for future development

See the [`openforcefield` GitHub issue tracker](https://github.com/openforcefield/openforcefield/issues) for a more thorough list, though some major areas are highlighted here:
* Exploring how use of partial bond order can simplify the environments needed
* Implementing partial bond order use for other parameter types
* Modifications discussed above to handling of torsions
* Possible modifications to make it easier to support other force field families
* [Off-atom charges](https://github.com/openforcefield/smarty/issues/132)
* Additional functional forms for vdW interactions

## Use for parameterization of systems

A relatively extensive set of examples is available under [`examples/`](https://github.com/openforcefield/openforcefield/tree/master/examples). Basic usage works as follows in python, however:

```python
# Read in a small molecule using the OpenEye toolkit
import oechem
mol = oechem.OEGraphMol()
ifs = oechem.oemolistream(mol_filename)
flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
ifs.SetFlavor(oechem.OEFormat_MOL2, flavor)
oechem.OEReadMolecule(ifs, mol)
oechem.OETriposAtomNames(mol)

# Load forcefield
from openforcefield.typing import smirnoff
forcefield = smirnoff.ForceField('forcefield/Frosst_AlkEthOH_parmAtFrosst.ffxml')

# Generate an openforcefield Topology and
from openforcefield.topology import Molecule, Topology
topology = Topology()
topology.add_molecule(mol)

# Create an OpenMM System from the openforcefield Topology
from simtk import openmm, unit
system = forcefield.createSystem(topology)
```
This example can essentially trivially be extended to handle the case of beginning from a SMILES string rather than a `.mol2` file.

The `SMIRNOFF_simulation` example in the examples directory shows how to extend the example above to simulate this molecule in the gas phase.

`createSystem()` can also handle a system consisting of a mixture of molecules; we've tested it on cyclohexane/ethanol and propane/methanol/butanol mixtures for example.
As input it is necessary to provide a Topology file representing the system, and a list of OpenEye molecules for the components of that Topology.
So, for example, one can read a PDB file describing a mixture and provide OpenEye molecules for the components (generated by the Mobleylab's [SolvationToolkit](https://github.com/MobleyLab/SolvationToolkit), for example) and create a system from that.

`createSystem()` allows the user to specify a choice of charge model, among other properties. Consult its help in python for more information.

One important note is that the OpenEye molecules currently must have atom names, hence the [`OETriposAtomNames`](https://docs.eyesopen.com/toolkits/python/oechemtk/OEChemFunctions/OETriposAtomNames.html) above.

**Use with protein force fields**: While SMIRNOFF format force fields can cover a wide range of biological systems, an initial focus is on small molecule force fields, meaning that users may have considerable interest in combining SMIRNOFF small molecule parameters to systems in combination with traditional biopolymer parameters from conventional force fields, such as the AMBER family of protein/nucleic acid force fields.
Thus, we provide an example of setting up a mixed system in [`examples/mixedFF_structure`](examples/mixedFF_structure), where an AMBER family force field is used for a protein and smirnoff99Frosst for a small molecule.

### The optional `id` attribute and other XML attributes

In general, additional optional XML attributes can be specified and will be ignored by `ForceField` unless they are specifically handled by the parser (and specified in this document).

One attribute we have found helpful in parameter file development is the `id` attribute for a specific parameter line, and we *recommend* that SMIRNOFF forcefields utilize this as effectively a parameter serial number, such as in:
```XML
 <Bond smirks="[#6X3:1]-[#6X3:2]" id="b5" k="820.0" length="1.45"/>
```
Some functionality in `ForceField`, such as `ForceField.labelMolecules`, looks for the `id` attribute.
Without this attribute, there is no way to uniquely identify a specific parameter line in the XML file without referring to it by its smirks string, and since some smirks strings can become long and relatively unwieldy (especially for torsions) this provides a more human- and search-friendly way of referring to specific sets of parameters.

The `parent_id` attribute is also frequently used to denote parameters from which the current parameter is derived in some manner.

### A remark about parameter availability

`ForceField` will currently raise an exception if any parameters are missing where expected for your system---i.e. if a bond is assigned no parameters, an exception will be raised.
However, use of generic parameters (i.e. `[*:1]~[*:2]` for a bond) in your FFXML will result in parameters being assigned everywhere, bypassing this exception.
So use generics sparingly unless it is your intention to provide generics that should be used.

## Version history

### 0.2

Overhaul of draft specification along with `ForceField` refactor:
* Aromaticity model now defaults to `MDL`, and aromaticity model names drop OpenEye-specific prefixes
* Added a description of default units if `*_unit` attributes are omitted.
* Forces were renamed to be more general:
    * `<NonbondedForce>` was renamed to `<StericsForce>`
    * `<HarmonicBondForce>` was renamed to `<BondForce>`
    * `<HarmonicAngleForce>` was renamed to `<AngleForce>`
* `<PeriodicTorsionForce>` was split into `<ProperTorsionForce>` and `<ImproperTorsionForce>`
* `<StericsForce>` now specifies 1-2, 1-3, and 1-4 scaling factors via `scale12` (default: 0), `scale13` (default: 0), and `scale14` (default: 0.5) attributes. Coulomb scaling parameters have been removed from `StericsForce`.
* Added the `<ElectrostaticsForce>` tag to separately specify 1-2, 1-3, and 1-4 scaling factors for electrostatics, as well as the method used to compute electrostatics (`PME`, `reaction-field`, `Coulomb`) since this has a huge effect on the energetics of the system.
* Made it clear that `Constraint` entries do not have to be between bonded atoms.
* The `potential` attribute was added to most forces to allow flexibility in extending forces to additional functional forms (or algebraic expressions) in the future. `potential` defaults to the current recommended scheme if omitted.
* `GBSAForce` now has defaults specified for `gb_method` and `sa_method`
* Changes to how fractional bond orders are handled:
    * Use of fractional bond order is now are specified at the force tag level, rather than the root level
    * The fractional bond order method is specified via the `fractional_bondorder_method` attribute
    * The fractional bond order interpolation scheme is specified via the `fractional_bondorder_interpolation`
* Section heading names were cleaned up.
* Example was updated to reflect use of the new `openforcefield.topology.Topology` class
* Eliminated "Requirements" section, since it specified requirements for the software, rather than described an aspect of the SMIRNOFF specification

### 0.1

Initial draft specification.
