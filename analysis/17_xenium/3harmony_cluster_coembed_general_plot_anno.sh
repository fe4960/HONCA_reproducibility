#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_general_plot_anno.R --no-save > HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_general_plot_anno.out 2> HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_general_plot_anno.err 
