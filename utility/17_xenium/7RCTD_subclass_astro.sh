#!/bin/sh
source ~/.condainit
conda activate spatial1

#err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/7RCTD_subclass_astro1
err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/7RCTD_subclass_astro

mkdir $err

ref=$1
bc=$3
lab=$2
meta=$4 
ds=$5
clu=$6
id=$7

R --vanilla --args $ref $lab $bc $meta $ds $clu $id < ${err}.R > ${err}/${lab}_${id}.err 2> ${err}/${lab}_${id}.out
#Rscript --vanilla ${err}.R $ref $lab $bc $meta $ds $clu $id > ${err}/${lab}_${id}.err 2> ${err}/${lab}_${id}.out
######Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/7RCTD.R $ref $lab $bc $meta $ds $clu  > ${err}/${lab}.err 2> ${err}/${lab}.out
