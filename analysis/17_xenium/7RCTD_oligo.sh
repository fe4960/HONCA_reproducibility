#!/bin/sh
source ~/.condainit
conda activate spatial1

#err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/7RCTD_oligo
err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/7RCTD_oligo1

mkdir $err

ref=$1
bc=$3
lab=$2
meta=$4 
ds=$5
clu=$6
Rscript --vanilla ${err}.R $ref $lab $bc $meta $ds $clu  > ${err}/${lab}.err 2> ${err}/${lab}.out
######Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/7RCTD.R $ref $lab $bc $meta $ds $clu  > ${err}/${lab}.err 2> ${err}/${lab}.out
