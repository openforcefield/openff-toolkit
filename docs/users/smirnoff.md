# SMIRNOFF (SMIRks Native Open Force Field)

## The SMIRNOFF specification

The SMIRNOFF specification can be found in the OpenFF [standards repository].

## SMIRNOFF and the Toolkit

OpenFF releases all its force fields in SMIRNOFF format. SMIRNOFF is a format
developed by OpenFF; its specification can be found in our
[standards repository]. SMIRNOFF-format force fields are distributed as XML
files with the `.offxml` extension. Instead of using atom types like
traditional force field formats, SMIRNOFF associates parameters directly with
chemical groups using [SMARTS] and [SMIRKS], which are extensions of the
popular SMILES serialization format for molecules. SMIRNOFF goes to great
lengths to ensure reproducibility of results generated from its force fields.

The OpenFF Toolkit is the reference implementation of the SMIRNOFF spec. The
toolkit is responsible for reading and writing `.offxml` files, for
facilitating their modification, and for applying them to a molecular system in
order to produce an [`Interchange`] object. The OpenFF Interchange project then
takes over and is responsible for [producing input files and data] for actual
MD software. The toolkit strives to be backwards compatible with old versions
of the spec, but owing to the vagaries of the arrow of time cannot be forward
compatible. Trying to use an old version of the toolkit to load an `.OFFXML`
file created with a new version of the spec will lead to an error.

A simplified `.offxml` file for TIP3P water might look like this:

```xml
<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">
    <Author>The Open Force Field Initiative</Author>
    <Date>2021-08-16</Date>
    <Constraints version="0.3">
        <Constraint smirks="[#1:1]-[#8X2H2+0:2]-[#1]" id="c-tip3p-H-O" distance="0.9572 * angstrom"></Constraint>
        <Constraint smirks="[#1:1]-[#8X2H2+0]-[#1:2]" id="c-tip3p-H-O-H" distance="1.5139006545247014 * angstrom"></Constraint>
    </Constraints>
    <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1.0" cutoff="9.0 * angstrom" switch_width="1.0 * angstrom" method="cutoff">
        <Atom smirks="[#1]-[#8X2H2+0:1]-[#1]" epsilon="0.1521 * mole**-1 * kilocalorie" id="n-tip3p-O" sigma="3.1507 * angstrom"></Atom>
        <Atom smirks="[#1:1]-[#8X2H2+0]-[#1]" epsilon="0 * mole**-1 * kilocalorie" id="n-tip3p-H" sigma="1 * angstrom"></Atom>
    </vdW>
    <Electrostatics version="0.3" scale12="0.0" scale13="0.0" scale14="0.8333333333" scale15="1.0" cutoff="9.0 * angstrom" switch_width="0.0 * angstrom" method="PME"></Electrostatics>
    <LibraryCharges version="0.3">
        <LibraryCharge smirks="[#1]-[#8X2H2+0:1]-[#1]" charge1="-0.834 * elementary_charge" id="q-tip3p-O"></LibraryCharge>
        <LibraryCharge smirks="[#1:1]-[#8X2H2+0]-[#1]" charge1="0.417 * elementary_charge" id="q-tip3p-H"></LibraryCharge>
    </LibraryCharges>
</SMIRNOFF>
```

:::{note} TIP3P's geometry is specified entirely by constraints, but SMIRNOFF
   certainly supports a wide variety of bonded parameters and functional
   forms.
:::

Note that this format specifies not just the individual parameters, but also their
functional forms and units in very explicit terms. This both makes it easy to read
and means that the correct implementation of each force is specifically defined,
rather than being left up to the MD engine.

The complicated part is that each parameter is specified by a SMIRKS code. These
codes are SMARTS codes with an optional numerical index on some atoms given
after a colon. This indexing system comes from SMIRKS. Each parameter expects a
certain number of indexed atoms, and applies the force accordingly. Unindexed
atoms are used to match the chemistry, but forces are not applied to them.
SMARTS/SMIRKS codes are less intimidating than they look; `[#1]` matches any
Hydrogen atom (atomic number 1), while `[#8X2H2+0]` matches an oxygen atom
(atomic number 8) with some additional constraints. Dashes represent bonds. So
`[#1]-[#8X2H2+0:1]-[#1]` represents an oxygen atom indexed as 1 connected to
two unindexed hydrogen atoms. This system allows individual parameters to be as
general or as specific as needed.

[SMARTS]: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
[SMIRKS]: https://www.daylight.com/dayhtml/doc/theory/theory.smirks.html
[standards repository]: https://openforcefield.github.io/standards/standards/smirnoff/
[`Interchange`]: openff.interchange.Interchange
[producing input files and data]: openff.interchange:using/output

:::{hint} 
This page is not the SMIRNOFF spec; it has been moved to the
[standards repository].
:::