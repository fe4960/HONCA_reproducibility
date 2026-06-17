#!/bin/sh
source ~/.condainit
conda activate veloc

h5ad=$2
indir=$1
bname=$3
hvg=$4

python HCA_ON/scripts/6_velocyto/3scvelo_general.py $indir $h5ad $bname $hvg 
