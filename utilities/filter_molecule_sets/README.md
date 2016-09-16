# Filtering Molecule Sets

There are a number of molecule sets that have been used in the open-forcefield effort. 
In an effort to test the SMIRFFparm99Frosst it became increasingly evident that some of these molecules sets include molecules outside what we want to currently test our tools. 

The main script for this directory is `filter_molecule_sets.py` 

```
Usage:             Given a set of molecules filter for
            number of heavy atoms per molecule,
            number of metals atoms per molecule,
            inappropriate valence (carbons or nitrogens with 5+ bonds),
            SMIRKS patterns you do not wish to have included

    usage:
        filter_molecule_sets.py --input DrugBank.sdf --output updated_DrugBank.mol2.gz \
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

```
