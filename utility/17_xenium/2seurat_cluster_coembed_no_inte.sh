#!/bin/sh
source ~/.condainit
conda activate spatial

R < HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_no_inte.R --no-save > HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_no_inte.out 2> HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_no_inte.err 
