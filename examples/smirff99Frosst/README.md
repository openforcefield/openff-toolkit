# Manual conversion of amber-parm99+parm@Frosst to smirff format

Amber parm99 with a multitude of additional parameters in parm@Frosst can be effectively used on a wide spectrum of pharmaceutically relevant small molecules. An attempt is made here to convert these into smirff format.

## Manifest
Only a few of the files here are intended for use, most are intermediate files used to generate the final useful ones. The intermediate ones are here in case a forensic audit is needed for certain parameters.
The most useful files are:

smirffishFrcmod.parm99Frosst.txt - the key source file for the smirff ffxml file, this file has the format of an amber .frcmod file except the initial field giving bonded amber atomtypes is replaced with a smirks string for use with the smirff format for that parameter. The translation of this file into smirff99Frosst.ffxml is done in ipython by parameter_usage.ipynb.

smirff99Frosst.ffxml - this is the smirff format file containing the best efforts of C.I.Bayly during his sabbatical in summer 2016 to convert amber parm99.dat plus the parm@Frosst frcmod file into a single small-molecule smirff ffxml.

parameter_usage.ipynb - this ipython notebook originally written by David Mobley carries out several independent steps:
a) Reading in a smirffishFrcmod file (here smirffishFrcmod.parm99Frosst.txt) and converting it to smirff ffxml format.
b) Reading in a database of molecules and applying to them the smirff ffxml from part (a).
c) Analyzing the result of (b) to see which molecules get which parameters, including a depiction function which will depict a molecule highlighting the atoms associated with a specific parameter.

The supporting intermediate files consist mostly of the building up of smirffishFrcmod.parm99Frosst.txt by gradually adding Bonds, VdW, and Angles in that order. Then Impropers (a very few parameters) were added. With Torsions, the intermediate files were done first for parm99 and then parm@Frosst. Then the Torsions were combined into smirffishFrcmod.Torsions.parm99withParmFrosst.txt and harmonized. This file was then introduced into smirffishFrcmod.parm99Frosst.txt, replacing the parm99-only torsions. This was then iteratively improved based on how well the parm@Frosst zinc 7505 molecule set was parameterized in parameter_usage.ipynb. With the 2016sep06 version, there were only 17 molecules requiring thegeneric torsion, most of those having broken chemistry (e.g. pentavalent nitrogens or carbons).

