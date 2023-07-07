# GROMACS Data files

* `rb_torsion.top`
* `rb_torsion.gro`

[Source](https://github.com/shirtsgroup/InterMol/blob/v0.1.2/intermol/tests/gromacs/unit_tests/dihedral3_vacuum/dihedral3_vacuum.top)

Round-tripped through ParmEd to creata a monolithic file without using i.e. `[ dihedraltypes ]`
Manually modified `comb-rule` to 2 since that doesn't matter for testing RB torsions
