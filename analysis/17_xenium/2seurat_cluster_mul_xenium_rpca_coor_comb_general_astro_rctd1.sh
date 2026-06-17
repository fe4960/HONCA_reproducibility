#!/bin/sh
source ~/.condainit
conda activate spatial1


cut=$2

########
#####process xenium data

label=$1

f=$3
clu=$4
sc=$5
idx=$6
label1=$7
err="HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general_astro_rctd1"
#err="HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general_astro_rctd1_x"
mkdir $err
#Rscript --vanilla HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general_astro_rctd1.R   $label $cut $f $clu $sc $idx > ${err}/${label}_${cut}_${f}_${clu}_${sc}_${idx}.err 2> ${err}/${label}_${cut}_${f}_${clu}_${sc}_${idx}.out


R --vanilla --args   $label $cut $f $clu $sc $idx $label1 <  ${err}.R  > ${err}/${label}_${cut}_${f}_${clu}_${sc}_${idx}_${label1}.err 2> ${err}/${label}_${cut}_${f}_${clu}_${sc}_${idx}_${label1}.out
