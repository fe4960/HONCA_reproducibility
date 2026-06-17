#!/bin/sh

source ~/.condainit
conda activate scvi
h5ad=$1
l=$2
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_hotspot/2_donor_100p_metacell_leiden.py $h5ad $l
