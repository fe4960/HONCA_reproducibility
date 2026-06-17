#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general_anno_RCTD_astro.R --no-save > HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general_anno_RCTD_astro.out 2> HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general_anno_RCTD_astro.err 
