#!/bin/sh
source ~/.condainit
conda activate scvi
dir=$1
h5ad=$2
list1=$3
list2=$4
label=$5
ct_list=$6
python HCA_ON/scripts/7_DEG/1gene_count_DEG.py ${dir} ${h5ad} ${list1} ${list2} ${label} ${ct_list}
