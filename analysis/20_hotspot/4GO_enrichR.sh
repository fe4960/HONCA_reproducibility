#!/bin/sh
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/20_hotspot/4GO_enrichR.R
source ~/.condainit
conda activate qtl
main=$1
ct=$2
bg1=$3
t=$4
Rscript --vanilla ${comd} $main $ct $bg1 $t
