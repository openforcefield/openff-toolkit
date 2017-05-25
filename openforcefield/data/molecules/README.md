# Molecule Files

This is a list of the molecule files provided here
including a brieg description when possible

* `AlkEthOH-tripos.mol2.gz`: set of 42 alkanes, ethers, and alcohols
* `*.crd`, `*.mol2`, and `*.top` for the following AlkEthOH molecules
        - `AlkEthOH_c100`
        - `AlkEthOH_c1161`
        - `AlkEthOH_c1266`
        - `AlkEthOH_c38`
        - `AlkEthOH_r0`
        - `AlkEthOH_r118`
        - `AlkEthOH_r12`
* `*.crd`, `*.mol2`, and `*.top` for benzene
* `AlkEthOH_test_filt1_*.mol2`: A different set of 42 alkanes ethers and alcohols
* `zinc-subset-*.mol2.gz`: set of 7505 molecules
* `PhEthOH_pFrosstTyped.oeb`: 5082 molecules with phenyl groups in addition to alkanes, ethers, and alcohols
* `PhEthOH_pFrosstTyped_first200.oeb`: first 200 in PhEthOH set used for testing smarty and smirky
* `DrugBank_atyped.oeb`: 7133 parm@Frosst atom typed molecules from the [DrugBank Release Version 5.0.1](http://www.drugbank.ca/releases/latest)
* `DrugBank_*.mol2`: 5928 molecules with 3D coordinates generated with `utilities/filter_molecules/coordinates_for_DrugBank.py` from `DrugBank_atyped.oeb` 
* `MiniDrugBank_*.mol2`: A subset of `DrugBank_*.oeb`. A version record is available on gitHub at [open-forcefield-group/MiniDrugBank](https://github.com/open-forcefield-group/MiniDrugBank) 
**For all molecule sets**
        - `*` = `ff` for parm@frosst atom types 
        - `*` = `tripos` for tripos atom types 
          
