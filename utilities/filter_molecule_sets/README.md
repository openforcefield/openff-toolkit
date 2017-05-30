# Filtering Molecule Sets

There are a number of molecule sets that have been used in the open-forcefield effort. 
In an effort to test the SMIRNOFFparm99Frosst it became increasingly evident that some of these molecules sets include molecules outside what we want to currently test our tools. 
Molecule sets and relevant scripts are listed below. 
For our purposes DrugBank refers to [DrugBank Release Version 5.0.1](http://www.drugbank.ca/releases/latest).
The molecules in this set where then typed with parm@frosst atom types by Christopher I. Bayly and stored as `openforcefield/data/molecules/DrugBank_atyped.oeb`

**Molecule Sets**
* `DrugBank_updated_tripos.mol2.gz` - created by calling 
```
python filter_molecule_sets.py --input DrugBank_atyped.oeb \
--output DrugBank_updated_tripos.mol2.gz --repeats False --warnings False \
--heavy 100 --SMIRKS remove_smirks_simple.smarts --metals 0 \
--atoms elements_exclude.txt --type gg,Se1 --flavor tripos
```
* `DrugBank_updated_ff.mol2.gz` - created by calling 
```
python filter_molecule_sets.py --input DrugBank_atyped.oeb \
--output DrugBank_updated_ff.mol2.gz --repeats False --warnings False \
--heavy 100 --SMIRKS remove_smirks_simple.smarts --metals 0 \
--atoms elements_exclude.txt --type gg,Se1 --flavor ff
```

**Python Scripts**
* `filter_molecule_sets.py` - This script was developed to filter unwanted molecules out of a larger set and was created to be as general as possible. See details below. 
* `coordinates_for_DrugBank.py` - This script specifically starts with `openforcefield/data/molecules/DrugBank_atyped.oeb` and generates 3D coordinates for as many molecules as possible. It generates `DrugBank_ff.mol2` and `DrugBank_tripos.mol2` that are available in the `data/molecules/` directory. Note that a time limit for conformer generation was created to prevent memory crashes on the UC Irvine cluster, 289 of 7133 initial molecules ran out of the designated time, oeomega failed to generate 3D coordinates for another 916 molecules.  
* `oeb-to-FF-and-tripos-mol2.py` - utility to convert oeb file with parm@Frosst atomtypes to two `*.mol2` files with tripos and parm@frosst atomtypes

**Input Files**
* `elements_exclude.txt` - File with element numbers that you do not want in your set of molecules.
* `remove_smirks_simple.txt` - Example input file, any molecules containing this SMIRKS pattern will be removed


 The main python script for this directory is `filter_molecule_sets.py` 

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

