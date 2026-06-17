#!/bin/sh
source ~/.condainit
conda activate hotspot

dir=$1
label=$2

python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_hotspot/3_hotspot_heatmap.py $dir $label
