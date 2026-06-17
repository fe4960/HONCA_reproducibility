#!/bin/sh
source ~/.condainit
conda activate scvi

python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/6refine_scvi_plot_MALAT1.py $1 $2 $3 $4
