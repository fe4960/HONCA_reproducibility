#!/bin/sh
source ~/.condainit
conda activate spatial1
R < HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_prev_flt_plot.R --no-save > HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_prev_flt_plot.out 2> HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_prev_flt_plot.err
