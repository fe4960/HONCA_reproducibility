#!/bin/sh
source ~/.condainit
conda activate scvi
python HCA_ON/scripts/11_spatial/4overlap_xenium_nsforst_bin_gene_general_major.py $1 $2 $3
#python HCA_ON/scripts/11_spatial/4overlap_xenium_nsforst_bin_gene_ononh.py $1 $2 $3
