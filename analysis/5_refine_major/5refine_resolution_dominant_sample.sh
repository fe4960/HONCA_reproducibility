#!/bin/sh

source ~/.condainit
conda activate scvi

#h5ad=$1
#res=$2
#batch_key=$3
bname=$1
outdir=$2

python HCA_ON/scripts/5_refine_major/5refine_resolution_dominant_sample.py $bname $outdir
