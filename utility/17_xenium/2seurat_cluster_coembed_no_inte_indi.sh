#!/bin/sh
source ~/.condainit
conda activate spatial

err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_no_inte_indi

mkdir $err

ref=$1
meta=$2
bc=$3
lab=$4
ds=$5
Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_no_inte_indi.R $ref $meta $bc $lab $ds > ${err}/${lab}.err 2> ${err}/${lab}.out
