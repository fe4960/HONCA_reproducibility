#!/bin/sh
source ~/.condainit
conda activate spatial1

comd="HCA_ON/scripts/17_xenium/7RCTD_mk"

mkdir $comd

label=$1

cut=$2
fc=$3
sc=$4
idx=$5

Rscript --vanilla ${comd}.R $label $cut $fc $sc $idx > ${comd}/${label}_${cut}_${fc}_${sc}_${idx}.out  2> ${comd}/${label}_${cut}_${fc}_${sc}_${idx}.err
