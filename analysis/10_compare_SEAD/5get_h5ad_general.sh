#!/bin/sh
source ~/.condainit
conda activate scvi
cell=$1
fn0=$2
fn1=$3
subclass=$4
python HCA_ON/scripts/10_compare_SEAD/5get_h5ad_general.py $cell $fn0 $fn1 $subclass
