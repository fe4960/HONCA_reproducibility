#!/bin/sh

source ~/.condainit
conda activate scvi
h5ad=$1
cell=$2
l=$3
n=$4
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_hotspot/2_donor_100p_metacell_permut.py $h5ad $cell $l $n
