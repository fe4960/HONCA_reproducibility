#!/bin/sh
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/14GO_enrichR_sub.R 
source ~/.condainit
conda activate qtl
ct=$1
cell=$2
od=$3
dem=$4
Rscript --vanilla ${comd} $ct $cell
