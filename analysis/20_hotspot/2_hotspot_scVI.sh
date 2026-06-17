#!/bin/sh
source ~/.condainit
conda activate hotspot

dir=$1
label=$2
hvg=$3
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_hotspot/2_hotspot_scVI.py $dir $label $hvg
