#!/bin/sh
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/12GO_analysis_mash_ct_DEG.R
source ~/.condainit
conda activate qtl
main=$1
ct=$2
j=$3
dem=$4
Rscript --vanilla ${comd} $main $ct $j $dem
