#!/bin/sh
source ~/.condainit
conda activate scvi
bname=$1
ct=$2

#python HCA_ON/scripts/5_refine_major/16match_sanses_subclass_major.py $bname $ct
python HCA_ON/scripts/5_refine_major/16match_sanses_subclass.py $bname $ct
