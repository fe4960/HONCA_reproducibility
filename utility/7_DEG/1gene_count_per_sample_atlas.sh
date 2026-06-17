#!/bin/sh
  
source ~/.condainit

conda activate scvi

python HCA_ON/scripts/7_DEG/1gene_count_per_sample_atlas_new.py
