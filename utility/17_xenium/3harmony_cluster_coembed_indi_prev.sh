#!/bin/sh
source ~/.condainit
conda activate spatial1

#err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_prev_plot
err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_prev_flt

mkdir $err

ref=$1
meta=$2
bc=$3
lab=$4
ds=$5
cut=$6
f=$7
sc=$8
#$rds $meta $bc $label $ds

#Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_prev.R $ref $lab $bc   > ${err}/${lab}.err 2> ${err}/${lab}.out

Rscript --vanilla ${err}.R $ref $meta $bc $lab $ds $cut $f $sc  > ${err}/${lab}.err 2> ${err}/${lab}.out
