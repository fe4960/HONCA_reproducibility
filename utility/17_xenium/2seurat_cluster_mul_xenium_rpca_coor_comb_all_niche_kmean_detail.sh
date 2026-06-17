#!/bin/sh
source ~/.condainit
conda activate spatial1


R < HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_all_niche_kmean_detail.R --no-save > HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_all_niche_kmean_detail.out 2> HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_all_niche_kmean_detail.err
