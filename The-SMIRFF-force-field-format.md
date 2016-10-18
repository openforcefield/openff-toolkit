# The Smirks Force Field (SMIRFF) format and its use

The SMIRFF format is based on the [`OpenMM`](http://openmm.org) [`ForceField`](http://docs.openmm.org/7.0.0/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField) class and provides an XML format for encoding force fields based on [SMIRKS](http://www.daylight.com/dayhtml/doc/theory/theory.smirks.html)-based chemical perception. While designed for [`OpenMM`](http://openmm.org), parameters encoded in this format can be applied to systems and then these systems converted via [`ParmEd`](http://parmed.github.io/ParmEd) and [`InterMol`](https://github.com/shirtsgroup/InterMol) for simulations in a variety of other simulation packages.

## Basic structure

The SMIRFF format provides XML `ffxml` files that are parseable by the `ForceField` class of the `smarty.forcefield` module. These encode parameters for a force field based on a SMIRKS-based specification of the chemical environment the parameters are to be applied to. The file has tags corresponding to OpenMM force terms (`HarmonicBondForce`, `HarmonicAngleForce`, `PeriodicTorsionForce`, etc., as discussed in more detail below); these specify units used for the different constants provided for individual force terms, for example (see the [AlkEthOH example ffxml](https://github.com/open-forcefield-group/smarty/blob/master/smarty/data/forcefield/Frosst_AlkEtOH.ffxml)):
```XML
   <HarmonicAngleForce angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
```       
which introduces following `Angle` terms which will use units of degrees for the angle and kilocalories per mole per square radian for the force constant.

Under each of these force terms, there are tags for individual parameter lines such as these:
```XML
   <Angle smirks="[a,A:1]-[#6X4:2]-[a,A:3]" angle="109.50" k="100.0" id="a0001" parent_id="a0001"/>
   <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.50" k="70.0" id="a0002" parent_id="a0001"/>
```
The first of these specifies the `[a,A:1]-[#6X4:2]-[a,A:3]` SMIRKS pattern for an angle, with a tetravalent carbon at the center with single bonds to two atoms of any type. Atoms are labeled 1, 2, and 3, with 2 being the central atom. Equilibrium angle values are provided, along with force constants (with units as given above). Additional tags `id` and `parent_id` are essentially "informational", with `id` essentially serving as a serial number for the parameter (it can be used to report on which parameters are assigned to a given molecule, for example) and `parent_id` ultimately being helpful in expressing relationships between parameters and tracking *parameter hierarchies*, though we do very little with this yet.

**SMIRFF parameters are hierarchical** in that parameters which come later in a file override those which come earlier if they match the same pattern. This can be seen in the example above, where the first line provides a generic angle parameter for any tetravalent carbon (single bond) angle, and the second line overrides this for the specific case of a hydrogen-(tetravalent carbon)-hydrogen angle.

This hierarchical structure means that a typical parameter file will tend to have generic parameters early in the section for each force type, with more specialized parameters assigned later.

**Technical note**: Because this is an XML format, certain special characters that occur in valid SMIRKS patterns (such as ampersand symbols `&`) must be treated specially to process.

## Functional forms, etc.

**Functional form**: The SMIRFF format specifies parameters; once specified, these are processed by the SMIRFF `ForceField` class and used to assign parameters to OpenMM Forces. This means that specific forces are generally implemented as discussed in the [OpenMM Documentation](http://docs.openmm.org/7.0.0/userguide/theory.html), see especially [Section 19 on Standard Forces](http://docs.openmm.org/7.0.0/userguide/theory.html#standard-forces) for functional forms. In some cases, typically for consistency with the AMBER force field philosophy motivating some of the authors, we do some manipulation of parameters from these files as discussed below in "Parameter sections".

**Charges**: In keeping with the AMBER force field philosophy, especially as implemented in small molecule force fields such as [GAFF](http://ambermd.org/antechamber/gaff.html), [GAFF2](https://mulan.swmed.edu/group/gaff.php), and [parm@Frosst](http://www.ccl.net/cca/data/parm_at_Frosst/), we at least initially treat partial charges as something to be obtained separately from the rest of the force field (bonds, angles, impropers, torsions [BAIT] and vdW terms), typically via QM calculations or a method such as Bayly's [AM1-BCC](https://dx.doi.org/10.1002/jcc.10128) approach, thus, for system setup we provide the option of specifying a charging method, though charges are not normally specified in the FFXML itself. For other force fields with "library"-style charges, we could introduce a new section providing specific charges via SMIRKS pattern, though this is not yet implemented. Bond charge corrections, however, are supported as we discuss below.

## Parameter sections

For this section it will help to have on hand an example SMIRFF file, such as that the [AlkEthOH example ffxml](https://github.com/open-forcefield-group/smarty/blob/master/smarty/data/forcefield/Frosst_AlkEtOH.ffxml) or the larger prototype [SMIRFF99Frosst ffxml](https://github.com/open-forcefield-group/smirff99Frosst/blob/master/smirff99Frosst.ffxml).

Before getting in to individual sections, it's worth noting that the XML parser ignores attributes in the XML that it does not understand, so providing a parameter line for an angle that specifies (for example) a second force constant `k2` will lead to no effect.

### NONBONDED PARAMETERS (E.G. LENNARD-JONES)

Nonbonded parameters (currently, Lennard-Jones parameters) are specified via the [`NonbondedForce`](http://docs.openmm.org/7.0.0/userguide/theory.html#nonbondedforce) tag with sub-tags for individual `Atom` entries, such as:
```XML
<NonbondedForce coulomb14scale="0.833333" lj14scale="0.5" sigma_unit="angstroms" epsilon_unit="kilocalories_per_mole">
   <Atom smirks="[#1:1]" rmin_half="1.4870" epsilon="0.0157" id="n0001" parent_id="n0001"/>
   <Atom smirks="[$([#1]-[#6]):1]" rmin_half="1.4870" epsilon="0.0157" id="n0002" parent_id="n0001"/>
   ...
</NonbondedForce>
Scaling terms for 1-4 interactions should be specified in attributes for the `NonbondedForce` tag, along with units.

For compatibility, the size property of an atom can be specified either via providing the `sigma` attribute, such as `sigma="1.3"`, or via the `r_0/2` (`rmin/2`) values used in AMBER force fields (here denoted `rmin_half` as in the example above). The two are related by `r0 = 2^(1/6)*sigma` and conversion is done internally in `ForceField` into the `sigma` values used in OpenMM. `epsilon` denotes the well depth.

### BOND PARAMETERS

Bond parameters are specified via the [`HarmonicBondForce`](http://docs.openmm.org/7.0.0/userguide/theory.html#harmonicbondforce) tag with individual `Bond` tags providing equilibrium bond length `length` and force constant `k` values for specific bonds, for example:
```XML
<HarmonicBondForce length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
   <Bond smirks="[#6X4:1]-[#6X4:2]" length="1.526" k="620.0" id="b0001" parent_id="b0001"/>
   <Bond smirks="[#6X4:1]-[#1:2]" length="1.090" k="680.0" id="b0002" parent_id="b0002"/>
...
</HarmonicBondForce>
```

**AMBER functional forms define the force constant `k` in a manner that differs by a factor of two---we do not use that convention here, electing to use the standard harmonic definition `U(r) = (k/2)*(r-length)^2` instead.** Thus, comparing a SMIRFF file to a corresponding AMBER parameter file or `.frcmod` will make it appear that force constants here are twice as large as they ought to be.

### ANGLE PARAMETERS

Angle parameters are specified via the [`HarmonicAngleForce`](http://docs.openmm.org/7.0.0/userguide/theory.html#harmonicangleforce) tag with individual `Angle` tags providing parameters (equilibrium angle `angle` and force constant `k`), as in this example:
```XML
<HarmonicAngleForce angle_unit="degrees" k_unit="kilocalories_per_mole/radian**2">
   <Angle smirks="[a,A:1]-[#6X4:2]-[a,A:3]" angle="109.50" k="100.0" id="a0001" parent_id="a0001"/>
   <Angle smirks="[#1:1]-[#6X4:2]-[#1:3]" angle="109.50" k="70.0" id="a0002" parent_id="a0001"/>
...
</HarmonicAngleForce>
```

**AMBER functional forms drop the factor of 2 in the angle energy term, which we elect not to do here.** Thus, comparing a SMIRFF file to a corresponding AMBER parameter file or .frcmod will make it appear that force constants here are twice as large as they ought to be.

### PROPER TORSIONS

Torsions are implemented as a [`PeriodicTorsionForce`](http://docs.openmm.org/7.0.0/userguide/theory.html#periodictorsionforce) tag with child tags for `Proper` (discussed here) and `Improper` (discussed below) parameters, for example:
```
<PeriodicTorsionForce phase_unit="degrees" k_unit="kilocalories_per_mole">
   <Proper smirks="[a,A:1]-[#6X4:2]-[#6X4:3]-[a,A:4]" idivf1="9" periodicity1="3" phase1="0.0" k1="1.40" id="t0001" parent_id="t0001"/>
...
   <Proper smirks="[#6X4:1]-[#6X4:2]-[#8X2:3]-[#6X4:4]" idivf1="1" periodicity1="3" phase1="0.0" k1="0.383" idivf2="1" periodicity2="2" phase2="180.0" k2="0.1" id="t0008" parent_id="t0003"/>
...
</PeriodicTorsionForce>
```

Here, child `Proper` tags specify at least `k1`, `phase1`, and `periodicity1` attributes for the corresponding parameters of the first force term applied to this torsion. However, additional values are allowed in the form kN, phaseN, and periodicityN, where all N values must be consecutive (e.g. no `k1` and `k3` values without a `k2` value) but `N` can go as high as necessary.

Optionally, an `idivfN` attribute may be specified for each torsional term (for easier compatibility with AMBER files); this specifies a numerical value (in AMBER, always an integer) which is used as a divisor for the barrier height when assigning the torsion; i.e., a torsion with `idivf1="9"` is assigned a barrier height `k1` that is 1/9th the specified value. If `idivfN` is not specified, the barrier height is applied as stated.

In the future, we may switch to a model where torsional barriers are [automatically divided by the number of torsions along a bond](https://github.com/open-forcefield-group/smarty/issues/131), effectively resulting in the torsional barrier being the average of all barriers applied to that bond, rather than the current model where barriers are summed. (Barrier heights would need to be increased accordingly.) This would result in better handling of some cases where a small change in a molecule (such as a change in tautomer) could currently (as in AMBER) result in a dramatically different applied barrier height because of a change in the number of torsions passing through that bond. The averaging approach would make it easier to avoid this problem without requiring as many different torsional terms.

### IMPROPER TORSIONS

Impropers are applied in the same manner as proper torsions, via `PeriodicTorsionForce`, but with the `Improper` tag, as in:
```XML:
<Improper smirks="[*:1]~[#6X3:2](=[#7X2,#7X3+1:3])~[#7:4]" id="i4" k1="10.5" periodicity1="2" phase1="180."/>
...
```

**Improper torsions deviate profoundly from AMBER handling of impropers** in two ways. First, to eliminate ambiguity, we treat impropers as a [trefoil](https://upload.wikimedia.org/wikipedia/commons/thumb/5/5c/Trefoil_knot_left.svg/2000px-Trefoil_knot_left.svg.png) and apply the same set of parameters to all six paths around the trefoil. *Because of this, all barrier heights are divided by six before we apply them*, for consistency with AMBER force fields. Second, the *second* atom in an improper (in the example above, the trivalent carbon) is the central atom in the trefoil.

### SPECIAL SECTIONS

**Bond charge corrections**
Bond charge corrections (along the lines of Christopher Bayly's bond charge corrections in [AM1-BCC](https://dx.doi.org/10.1002/jcc.10128)) can be applied via a `BondChargeCorrections` tag with children specifying specific `BondChargeCorrection` terms. Here is an example not intended for actual use:
```XML
<BondChargeCorrections method="AM1" increment_unit="elementary_charge">
  <BondChargeCorrection smirks="[#6X4:1]-[#6X3a:2]" increment="+0.0073" id="c0001" parent_id="c0001"/>
  <BondChargeCorrection smirks="[#6X4:1]-[#6X3a:2]-[#7]" increment="-0.0943" id="c0002" parent_id="c0002"/>
  <BondChargeCorrection smirks="[#6X4:1]-[#8:2]" increment="+0.0718" id="c0003" parent_id="c0003"/>
</BondChargeCorrections>
```
The charge model specified must be a method understood by the OpenEye toolkits, and the charge correction `increment` (in units of proton charge) will be applied on top of this by subtracting `increment` from the atom tagged as 1 and adding it to the atom tagged as 2.

## Advanced features

Standard usage is expected to rely primarily on the features documented above and potentially new features. However, some advanced features are also currently supported.

### Partial bond orders
Partial bond orders can be used to allow interpolation of parameters. For example, these parameters:
```XML
<HarmonicBondForce length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2">
...
    <Bond smirks="[#6X3:1]-[#6X3:2]" id="b5" k="820.0" length="1.45"/>
    <Bond smirks="[#6X3:1]:[#6X3:2]" id="b6" k="938.0" length="1.40"/>
    <Bond smirks="[#6X3:1]=[#6X3:2]" id="b7" k="1098.0" length="1.35"/>
...
```
can be replaced by a single parameter line:
```XML
<HarmonicBondForce length_unit="angstroms" k_unit="kilocalories_per_mole/angstrom**2" fractional_bondorder="interpolate-linear">
...
    <Bond smirks="[#6X3:1]!#[#6X3:2]" id="b5" k_bondorder1="820.0" k_bondorder2="1098" length_bondorder1="1.45" length_bondorder2="1.35"/>
...
```
This allows specification of force constants and lengths for bond orders 1 and 2, and then interpolation between those based on the partial bond order. Currently the Wiberg bond order is used, which will be obtained automatically from an AM1-based charge calculation using the OpenEye toolkits if a beta version (or later) of the October 2016 toolkits is used.

Important usage notes:
* An interpolation scheme must be specified in the `HarmonicBondForce` attributes; currently only `interpolate-linear` is supported, though a spline interpolation may be preferable (this needs to be explored)
* If it is desired to use fractional bond orders, the introductory SMIRFF tag for the file must specify that the force field will use these via `<SMIRFF use_fractional_bondorder="True">` or similar. Otherwise, no partial bond orders will be obtained for possible later use.
* This feature is only implemented for bonds at present, though it needs to be extended to angles and torsions; possibly also it may have value for vdW parameters (which could vary depending on bond order) though this needs to be explored

### Aromaticity models

Before conduct SMIRKS substructure searches, molecules are prepared by applying one of OpenEye's aromaticity models, with the default model used unless otherwise requested. Alternate aromaticity models can be requested by the force field, such as
`<SMIRFF version="0.1" aromaticity_model="OEAroModel_MDL">` used by SMIRFF99Frosst (a choice by Christopher Bayly to simplify handling of certain heteroaromatic compounds). Any of the names of the [aromaticity models available in the OpenEye toolkit](https://docs.eyesopen.com/toolkits/python/oechemtk/aromaticity.html) can be used.

### Future plans

See the Issue tracker for a more thorough list, though some major areas are highlighted here:
* Exploring how use of partial bond order can simplify the environments needed
* Implementing partial bond order use for other parameter types
* Modifications discussed above to handling of torsions
* Possible modifications to make it easier to support other force field families
* [Off-atom charges](https://github.com/open-forcefield-group/smarty/issues/132)
* Additional functional forms for vdW interactions

## Use for parameterization of systems

A relatively extensive set of examples is available under [smarty/examples](https://github.com/open-forcefield-group/smarty/tree/master/examples). Basic usage works as follows in python, however:

```python
from simtk import openmm, unit
import numpy as np
import smarty
import oechem
mol = oechem.OEGraphMol()
ifs = oechem.oemolistream(mol_filename)
flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MOL2_Default | oechem.OEIFlavor_MOL2_Forcefield
ifs.SetFlavor(oechem.OEFormat_MOL2, flavor)
oechem.OEReadMolecule(ifs, mol)
oechem.OETriposAtomNames(mol)

# Load forcefield
forcefield = smarty.ForceField(smarty.get_data_filename('forcefield/Frosst_AlkEtOH_parmAtFrosst.ffxml'))

# Generate an OpenMM Topology and create an OpenMM System
topology = smarty.generateTopologyFromOEMol(mol)
system = forcefield.createSystem(topology, [mol])
```
This example can essentially trivially be extended to handle the case of beginning from a SMILES string rather than a `.mol2` file.

The SMIRFF_simulation example in the examples directory shows how to extend the example above to simulate this molecule in the gas phase.

`createSystem()` can also handle a system consisting of a mixture of molecules; we've tested it on cyclohexane/ethanol and propane/methanol/butanol mixtures for example. As input it is necessary to provide a Topology file representing the system, and a list of OpenEye molecules for the components of that Topology. So, for example, one can read a PDB file describing a mixture and provide OpenEye molecules for the components (generated by the Mobleylab's [SolvationToolkit](https://github.com/MobleyLab/SolvationToolkit), for example) and create a system from that.

`createSystem()` allows the user to specify a choice of charge model, among other properties. Consult its help in python for more information.

One important note is that the OpenEye molecules currently must have atom names, hence the [`OETriposAtomNames`](https://docs.eyesopen.com/toolkits/python/oechemtk/OEChemFunctions/OETriposAtomNames.html) above.


### A remark about parameter availability

`ForceField` will currently raise an exception if any parameters are missing where expected for your system---i.e. if a bond is assigned no parameters, an exception will be raised. However, use of generic parameters (i.e. `[*:1]~[*:2]` for a bond) in your FFXML will result in parameters being assigned everywhere, bypassing this exception. So use generics sparingly unless it is your intention to provide generics that should be used.


## Requirements

Requirements match those of the SMARTY package and will be provided in more detail on the main page for SMARTY in SMARTY's `README.md`. Currently, [OpenEye toolkits](http://www.eyesopen.com/toolkit-development) (free for academic use, but they require a license) are utilized for most of our chemistry. [OpenMM](http://openmm.org) is also required, as are a variety of other relatively standard python packages and other toolkits available via [`conda`](http://conda.pydata.org/docs/building/meta-yaml.html).

The easiest way to install SMARTY along with its dependencies is via `conda`:
```bash
conda config --add channels omnia
conda install --yes smarty
```
