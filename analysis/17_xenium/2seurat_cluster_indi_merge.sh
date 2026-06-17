#!/bin/sh
source ~/.condainit
conda activate spatial
err=HCA_ON/scripts/17_xenium/2seurat_cluster_indi_merge
mkdir $err
bc=$1
label=$2
Rscript --vanilla HCA_ON/scripts/17_xenium/2seurat_cluster_indi_merge.R $bc $label > ${err}/${label}.err 2> ${err}/${label}.out
