#!/bin/sh
source ~/.condainit
conda activate spatial1

#err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/7RCTD_subclass_oligo
err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/7RCTD_subclass_oligo1

mkdir $err

ref=$1
bc=$3
lab=$2
meta=$4 
ds=$5
clu=$6
id=$7
label=${lab}_${clu}_${id}

R --vanilla  --args $ref $lab $bc $meta $ds $clu $id < ${err}.R > ${err}/${label}.err 2> ${err}/${label}.out

#Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/7RCTD_subclass_oligo.R $ref $lab $bc $meta $ds $clu $id > ${err}/${label}.err 2> ${err}/${label}.out
######Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/7RCTD.R $ref $lab $bc $meta $ds $clu  > ${err}/${lab}.err 2> ${err}/${lab}.out
