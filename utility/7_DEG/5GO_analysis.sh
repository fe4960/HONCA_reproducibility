#!/bin/sh
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/5GO_analysis.R
source ~/.condainit
conda activate qtl
#main=$1
cell=$1
ct=$2
#j=$3
#dem=$4
Rscript --vanilla ${comd} $cell $ct
#Rscript --vanilla ${comd} $main $ct $j $dem
