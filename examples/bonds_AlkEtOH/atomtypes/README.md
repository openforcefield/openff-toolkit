# Atom type SMARTS components

## Formats

### Initial types

A `basetypes` file specifies the initial atom types used to initialize the sampler.

Comments beginning with `%` are ignored throughout the file.
Each line has the format
```
<SMARTS> <typename>
```
where `<SMARTS>` is an [OpenEye SMARTS string](https://docs.eyesopen.com/toolkits/cpp/oechemtk/SMARTS.html) and `<typename>` is a human-readable typename associated with that atom type.

Atom type definitions are hierarchical, with the last match in the file taking precedence over earlier matches.

For example, we could use the elemental base types:
```
% atom types
H    hydrogen
C    carbon
N    nitrogen
O    oxygen
F    fluorine
P    phosphorous
S    sulfur
Cl   chlorine
Br   bromine
I    iodine
```

### Decorators

A `decorators` file contains a list of SMARTS

Comments beginning with `%` are ignored throughout the file.
Each line has the format
```
<SMARTS> <decoratorname>
```
where `<SMARTS>` is an [OpenEye SMARTS string](https://docs.eyesopen.com/toolkits/cpp/oechemtk/SMARTS.html) and `<decoratorname>` is a human-readable typename associated with that decorator.

The SMARTS component is ANDed together (using the `&` operator) with a parent atom type to create a new proposed child atom type.
The human-readable `<decoratorname>` is appended (with a space) to the parent name to keep a human-readable annotation of the proposed child atom type.

### Substitutions

It is often convenient to define various tokens that are substituted for more sophisticated SMARTS expressions.

% Substitution definitions
% Format:
% <SMARTS> <replacement-string>

Comments beginning with `%` are ignored throughout the file.
Each line has the format
```
<SMARTS> <substitution-name>
```
where `<SMARTS>` is an [OpenEye SMARTS string](https://docs.eyesopen.com/toolkits/cpp/oechemtk/SMARTS.html) and `<substitution-name>` is the token that will be substituted for this.

For example, we could define some elemental substitutions along with some substitutions for halogens:
```
% elements
[#9]    fluorine
[#17]   chlorine
[#35]   bromine
[#53]   iodine

% halogens
[$smallhals,$largehals]     halogen
[$fluorine,$chlorine]       smallhals
[$bromine,$iodine]          largehals
```

The [`OESmartsLexReplace`](http://docs.eyesopen.com/toolkits/python/oechemtk/OEChemFunctions/OESmartsLexReplace.html) function is used to implement these replacements.

## Manifest
* `basetypes-elemental.smarts` - basetypes file with elemental atom types - this is a good choice to begin with
* `basetypes.smarts` - basetypes file with more sophisticated atom types
* `decorators.smarts` - `decorators` file with a variety of decorators
* `decorators-simple.smarts` - minimal `decorators` file for testing
* `substitutions.smarts` - minimal `substitutions` file
