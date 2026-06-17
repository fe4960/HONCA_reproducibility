#!/bin/sh
source ~/.condainit
conda activate spatial1


#err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_indi_prev_flt_astro_impute_ONHON

err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_indi_prev_flt_astro_impute_astro

mkdir $err

ref=$1
meta=$2
bc=$3
lab=$4
ds=$5
cut=$6
f=$7
sc=$8
i=$9
#$rds $meta $bc $label $ds

#Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_prev.R $ref $lab $bc   > ${err}/${lab}.err 2> ${err}/${lab}.out

R --vanilla  --args $ref $meta $bc $lab $ds $cut $f $sc $i < ${err}.R > ${err}/${lab}_${i}.err 2> ${err}/${lab}_${i}.out
