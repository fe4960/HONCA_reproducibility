#!/bin/sh
source ~/.condainit
conda activate veloc

h5ad=$1
indir=$2


python HCA_ON/scripts/6_velocyto/2general_cell_loom.py $h5ad $indir 
