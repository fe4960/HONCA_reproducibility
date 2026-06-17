#!/bin/sh
source ~/.condainit
conda activate spatial1


R < HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_all_niche.R --no-save > HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_all_niche.out 2> HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_all_niche.err
