#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_indi.R --no-save > HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_indi.out 2> HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_indi.err 
