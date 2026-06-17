#!/bin/sh
source ~/.condainit
conda activate scvi
python HCA_ON/scripts/20_finemap/12summarize_gene_upset.py
