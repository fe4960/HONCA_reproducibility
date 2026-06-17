#!/bin/sh
source ~/.condainit
conda activate archr
ct=$1
meta=$2

Rscript --vanilla HCA_ON/scripts/7_DEG/19celltype_proportion.R $ct $meta 
