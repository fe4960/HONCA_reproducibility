#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"
d=$1 #"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
input=$2 #"Oligodendrocyte_subclass_seurat"
fn="${d}/${input}"
label=$3 #"subclass"
cellnum=$4 #5000
out="${d}/${input}_${label}_5k"
co=$5 #"t"

#############sh ${main}/eye_QTL/scripts/2_GWAS/7downsample1.sh $fn $label $cellnum $out $co 

#############sh ${main}/HCA_ON/scripts/16_general/1save_obs_general.sh $out

s=${input}_${label}_5k_simple

############sh ${main}/HCA_ON/scripts/16_general/2h5ad2h5seurat_general_final.sh $s $d

sh ${main}/eye_QTL/scripts/2_GWAS/MAGMA/2run_MAGMA.sh $s $label $d


sh ${main}/HCA_ON/scripts/14_GWAS_gene/4summarize_MAGMA_result_list_general.sh $s

perl ${main}/HCA_ON/scripts/14_GWAS_gene/4summarize_MAGMA_result.pl $s

sh HCA_ON/scripts/14_GWAS_gene/5dotplot_clean1.sh $s

