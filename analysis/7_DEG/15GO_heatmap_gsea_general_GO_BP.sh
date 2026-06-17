#!/bin/sh
source ~/.condainit
conda activate scvi
python HCA_ON/scripts/7_DEG/15GO_heatmap_gsea_general_GO_BP.py $1 $2 $3 $4 $5
