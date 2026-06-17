#!/bin/sh
  
source ~/.condainit

conda activate scvi

python HCA_ON/scripts/7_DEG/1gene_count_per_donor_atlas_general_label.py $1 $2 $3 $4
