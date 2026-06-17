#!/bin/sh

source ~/.condainit
conda activate scvi

h5ad=$1
res=$2
batch_key=$3
bname=$4
outdir=$5
rmlist=$6
mk=$7
save=$8
indir=$9
python HCA_ON/scripts/5_refine_major/5refine_resolution.py $h5ad $res ${batch_key} $bname $outdir $rmlist $mk $save $indir
