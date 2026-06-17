#!/bin/sh
source ~/.condainit
conda activate spatial1

err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi

mkdir $err

ref=$1
#meta=$2
bc=$3
lab=$2
#ds=$5

Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi.R $ref $lab $bc   > ${err}/${lab}.err 2> ${err}/${lab}.out

#Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi.R $ref $meta $bc $lab $ds > ${err}/${lab}.err 2> ${err}/${lab}.out
