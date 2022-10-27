# Molecule Files

This is a list of the molecule files provided here including a brief description when possible

* `AlkEthOH_tripos.tar.gz`: `*.crd`, `*.top`, and `*.mol2` (with Tripos SYBYL atom types) for all the molecules in
the AlkEthOH dataset. All molecules were taken from
https://github.com/openforcefield/open-forcefield-data/blob/master/Model-Systems/AlkEthOH_distrib/inputfiles.tar.gz.
All the files in the archive can be accessed through `openff.toolkit.tests.utils:get_alkethoh_file_path()`. The archive
has the following structure:
    - `AlkEthOH_rings_filt1`: A collection of cyclic alkanes, alcohols, and ethers.
    - `AlkEthOH_chain_filt1`: Acyclic alkanes, alcohols, and ethers.
    - `AlkEthOH_test_filt1`: A smaller test set composed of the first 42 molecules in `AlkEthOH_rings_filt1/`.
    - `convert_to_sybyl.py`: A script used to convert the original data from GAFF to Tripos SYBYL atom types.
* `AlkEthOH_test_filt1_*.mol2`: A test set of 42 alkanes, ethers, and alcohols from the AlkEthOH set. This is file
contains all the molecules in `AlkEthOH_tripos.tar.gz/AlkEthOH_test_filt1`.
* `*.crd`, `*.mol2`, and `*.top` for benzene
* `MiniDrugBank_*.mol2`: A subset of `DrugBank_*.oeb`. A version record is available on gitHub at [openforcefield/MiniDrugBank](https://github.com/openforcefield/MiniDrugBank)
**For all molecule sets**
        - `*` = `ff` for parm@frosst atom types
        - `*` = `tripos` for tripos atom types

