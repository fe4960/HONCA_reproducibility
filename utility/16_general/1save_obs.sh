#!/bin/sh
source ~/.condainit
conda activate scvi

h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/oligo_opc_seurat_new"

python HCA_ON/scripts/16_general/1save_obs.py $h5ad
