# Filtering Molecule Sets

There are a number of molecule sets that have been used in the open-forcefield effort. 
In an effort to test the SMIRFFparm99Frosst it became increasingly evident that some of these molecules sets include molecules outside what we want to currently test our tools. 

Molecule sets currently included in this directory are listed below, although when filtering is finalized, molecule sets should be moved to smart/data/molecules 

* `DrugBank.sdf` - [DrugBank Release Version 5.0.1](http://www.drugbank.ca/releases/latest)
* `DrugBank_atyped.oeb` - complete DrugBank molecule set with parm@Frosst atomtypes provided by Christopher Bayly
* `updated_DrugBank.mol2.gz` - created by calling `python filter_molecule_sets.py -f DrugBank.sdf -o updated_DrugBank.mol2.gz`
* `DrugBank_updated_tripos.mol2.gz` - created by calling `python filter_molecule_sets.py --input DrugBank_atyped.oeb --output DrugBank_updated_tripos.mol2.gz --repeats False --warnings False --heavy 100 --SMIRKS remove_smirks_simple.smarts --metals 0 --atoms elements_exclude.txt --type gg,Se1 --flavor tripos`
* `DrugBank_updated_ff.mol2.gz` - created by calling `python filter_molecule_sets.py --input DrugBank_atyped.oeb --output DrugBank_updated_ff.mol2.gz --repeats False --warnings False --heavy 100 --SMIRKS remove_smirks_simple.smarts --metals 0 --atoms elements_exclude.txt --type gg,Se1 --flavor ff`

The main script for this directory is `filter_molecule_sets.py` 

```
Usage:             Given a set of molecules filter for
            number of heavy atoms per molecule,
            number of metals atoms per molecule,
            inappropriate valence (carbons or nitrogens with 5+ bonds),
            SMIRKS patterns you do not wish to have included

    usage:
        filter_molecule_sets.py --input DrugBank_atyped.oeb --output updated_DrugBank.mol2.gz \
        --repeats False --warnings False --heavy 100 --SMIRKS remove_smikrs_example.smarts \
        --metals 0 --hydrogens True
    

Options:
  -h, --help            show this help message and exit
  -f INPUT_FILE, --input=INPUT_FILE
                        Name of initial molecule file
  -o OUTPUT_FILE, --output=OUTPUT_FILE
                        Name of file to save filtered Molecules
  -r REPEATS, --repeats=REPEATS
                        If True allows for multiple molecules with the same
                        isomeric SMILES string, default = False OPTIONAL
  -w WARNINGS, --warnings=WARNINGS
                        If True keeps molecules that result in warning while
                        loading, default = False, OPTIONAL
  -n HEAVY, --heavy=HEAVY
                        Maximum number of heavy atoms per molecule, default =
                        100, OPTIONAL
  -s SMIRKS_FILE, --SMIRKS=SMIRKS_FILE
                        If not None, the file of SMIRKS are parsed and
                        molecules that match that pattern are removed. The
                        file should have a single column with SMIRKS/SMARTS
                        strings where commented lines begin with a '%',
                        default = None, OPTIONAL
  -m METALS, --metals=METALS
                        Maximum number of metals per molecule, default = 0,
                        OPTIONAL
  -H HYDROGENS, --hydrogens=HYDROGENS
                        If True, hydrogens are added to the output molecules
  -a ATOMS, --atoms=ATOMS
			File name with atomic number of the elements that you
			do not want in your set of molecules, OPTIONAL
  -t ATOMTYPE, --atomtype=ATOMTYPE
			Atom types that you do not want in your set of molecules,
			OPTIONAL (usage: --atomtype gg,Se1)
  -y FLAVOR, --flavor=FLAVOR
			Choose between two different flavors for atom types, 
			default= tripos (options: tripos or ff). Use .oeb file
			for ff option.
```

elements_exclude.txt - File with element numbers that you do not want in your
set of molecules.
