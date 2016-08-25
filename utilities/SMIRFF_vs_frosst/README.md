# Comparison of energies for molecules using parm@frosst parameters versus SMIRFF parameters

## Manifest
- `SMIRFF_vs_frosst.ipynb`: Main IPython/Jupyter notebook doing energy comparison, comparison of bonded energies, visualization of occurrences.
- `parm_at_Frosst.tgz`: Full archive from http://www.ccl.net/cca/data/parm_at_Frosst/ 
- `tleap_tools.py`: Utility tools for handling tleap conversion of parm@frosst parameters into AMBER prmtop and crd files, adapted from openmoltools
- `zinc-subset*.mol2.gz`: ZINC subset from CCL archive but in .mol2 format, with both TRIPOS and parm@frosst atom names, as provided by Christopher Bayly.
- `zinc_frosst_amber.tar.gz`: AMBER prmtop/crd files for zinc parm@frosst subset. Can be output by SMIRFF_vs_frosst.ipynb or loaded by it. 
