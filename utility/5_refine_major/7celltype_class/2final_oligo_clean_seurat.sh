#!/bin/sh
source ~/.condainit
conda activate scvi
python HCA_ON/scripts/5_refine_major/7celltype_class/2final_oligo_clean_seurat.py
