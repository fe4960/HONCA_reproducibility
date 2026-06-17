#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/2seurat_cluster_mul_coembed_rpca.R --no-save > HCA_ON/scripts/17_xenium/2seurat_cluster_mul_coembed_rpca.out 2> HCA_ON/scripts/17_xenium/2seurat_cluster_mul_coembed_rpca.err 
