#!/bin/sh
  
source ~/.condainit

conda activate scvi

python HCA_ON/scripts/7_DEG/1gene_count_per_donor_atlas_general.py $1 $2 $3 $4 $5
