#!/usr/bin/env bash
export TYPE=$1
smirky --molecules AlkEthOH_test_filt1_ff.mol2 \
    --typetag ${TYPE} \
    --initialtypes initial_${TYPE}.smarts \
    --smirff forcefield/Frosst_AlkEtOH.ffxml \
    --iterations 2 \
    --temperature 0.0 \
    --verbose True \
    --output ${TYPE}_answers
