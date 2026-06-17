#!/bin/sh
source ~/.condainit
conda activate magma_celltype

sam=$1
label=$2
ct=$3
Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/14_GWAS_gene/2run_MAGMA.R $sam $label $ct
