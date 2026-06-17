#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_no_inte_ref_xenium_update_plot2_new.R  --no-save > HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_no_inte_ref_xenium_update_plot2_new.out 2> HCA_ON/scripts/17_xenium/2seurat_cluster_coembed_no_inte_ref_xenium_update_plot2_new.err
