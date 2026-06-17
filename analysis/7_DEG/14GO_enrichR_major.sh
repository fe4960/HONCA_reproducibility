#!/bin/sh
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/14GO_enrichR_major.R 
source ~/.condainit
conda activate qtl
main=$1
ct=$2
od=$3
dem=$4
Rscript --vanilla ${comd} $main $ct $od $dem
