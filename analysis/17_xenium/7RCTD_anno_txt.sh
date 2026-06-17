#!/bin/sh
source ~/.condainit 
conda activate spatial1
comd="HCA_ON/scripts/17_xenium/7RCTD_anno_txt"
mkdir $comd

Rscript --vanilla ${comd}.R $1 $2 $3 > ${comd}/${1}_${2}_${3}.out > ${comd}/${1}_${2}_${3}.err


