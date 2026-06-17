#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_general_plot.R --no-save > HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_general_plot.out 2> HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_general_plot.err 
