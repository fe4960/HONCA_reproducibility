#!/bin/sh
source ~/.condainit
conda activate qtl

ct=$1
dir=$2
mkdir -p ${dir}
Rscript --vanilla HCA_ON/scripts/7_DEG/3norm_gene_exp_batch.R $ct $dir  
