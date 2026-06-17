#!/bin/sh
source ~/.condainit
conda activate spatial1


cut=$1 #50
#label="oligo"
label=$2
f=$3
#bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium"
#bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_w1_6"
bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_w1_8"

#Rscript --vanilla HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca.R $cut $label $bc $f


#Rscript --vanilla HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_new.R $cut $label $bc $f
 
R --vanilla --args $cut $label $bc $f < HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_new.R 

#> HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca.out 2> HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca.err 
