#!/bin/sh
source ~/.condainit
conda activate scvi
python HCA_ON/scripts/7_DEG/15GO_heatmap_gsea_general_major.py $1 $2 $3
