#!/bin/sh
source ~/.condainit
conda activate spatial1


cut=$2

########
#####process xenium data

label=$1

f=$3
clu=$4


err="HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general"

Rscript --vanilla HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general.R   $label $cut $f $clu > ${err}/${label}_${cut}_${f}_${clu}.err 2> ${err}/${label}_${cut}_${f}_${clu}.out
