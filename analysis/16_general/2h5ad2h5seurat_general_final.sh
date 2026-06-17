#!/bin/sh
source ~/.condainit
conda activate seurat4_convert
#conda activate archr
main="/dfs3b/ruic20_lab/junw42"
s=$1 
d=$2
########ct="Oligodendrocyte"
#s="oligo_opc_seurat_new_simple"
##########d="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/"
#s1="oligo_opc_seurat_new"
#######s="Oligodendrocyte_subclass_seurat"
#######s1="Oligodendrocyte_subclass_seurat"
Rscript --vanilla ${main}/HCA_ON/scripts/16_general/2h5ad2h5seurat_general_final.R $s $d # $s1


