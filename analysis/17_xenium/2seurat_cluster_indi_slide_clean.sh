#!/bin/sh
source ~/.condainit
conda activate spatial1
err=HCA_ON/scripts/17_xenium/2seurat_cluster_indi_slide_clean
mkdir $err
bc=$1
label=$2
Rscript --vanilla HCA_ON/scripts/17_xenium/2seurat_cluster_indi_slide_clean.R $bc $label > ${err}/${label}.err 2> ${err}/${label}.out
