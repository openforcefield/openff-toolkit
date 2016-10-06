#!/usr/bin/env bash
smirky --molecules AlkEthOH_test_filt1_ff.mol2 \
    --typetag Torsion \
    --atomORbases atom_OR_bases.smarts \
    --atomORdecors atom_OR_decorators.smarts \
    --atomANDdecors atom_AND_decorators.smarts \
    --bondORbase bond_OR_bases.smarts \
    --bondANDdecors bond_AND_decorators.smarts \
    --atomOddsFile atom_odds_forTorsions.smarts \
    --initialtypes initial_Torsions.smarts \
    --bondOddsFile bond_odds_forTorsions.smarts \
    --substitutions substitutions.smarts \
    --smirff forcefield/Frosst_AlkEtOH.ffxml \
    --iterations $1 \
    --temperature $2 \
    --output $3 \
    --verbose True \
