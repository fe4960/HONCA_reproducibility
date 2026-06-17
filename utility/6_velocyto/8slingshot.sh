#!/bin/sh
comd=HCA_ON/scripts/6_velocyto/8slingshot
mkdir $comd
source ~/.condainit
conda activate slingshot
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/
s=oligo_opc_seurat_new_simple
Rscript --vanilla ${comd}.R $dir $s > ${comd}/${s}.out 2> ${comd}/${s}.err
