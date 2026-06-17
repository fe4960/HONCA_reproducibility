#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/2seurat_cluster_mul_rest.R --no-save > HCA_ON/scripts/17_xenium/2seurat_cluster_mul_rest.out 2> HCA_ON/scripts/17_xenium/2seurat_cluster_mul_rest.err 
