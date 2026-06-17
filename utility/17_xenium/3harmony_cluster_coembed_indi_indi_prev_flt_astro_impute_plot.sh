#!/bin/sh
source ~/.condainit
conda activate spatial1

err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_indi_prev_flt_astro_impute_plot
#err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_indi_prev_flt_astro_impute_plot_ONHON


R  < ${err}.R --no-save > ${err}.err 2> ${err}.out
