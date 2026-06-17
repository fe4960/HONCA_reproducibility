#!/bin/sh
source ~/.condainit
conda activate qtl

ct=$1
dir=$2
mkdir -p ${dir}
Rscript --vanilla HCA_ON/scripts/7_DEG/5filter_sample.R $ct $dir
