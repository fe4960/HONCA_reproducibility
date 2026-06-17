#!/bin/sh
comd="HCA_ON/scripts/14_GWAS_gene/MAGMA_pipeline_general"

mkdir $comd
ct="Oligodendrocyte_precursor_cell"
d="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/"
input=${ct}_subclass_seurat_cycling  #"Oligodendrocyte_subclass_seurat"
fn="${d}/${input}"
label="subclass"
cellnum=5000
#out="${d}/${input}_${label}_5k"
co="t"


##########sh sb.sh -c ${comd}.sh -e $comd -j $ct -p fg -m 60 -t 0-2 -- $d $input $label $cellnum $co


ct="Fibroblast"
d="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/"
input=${ct}_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location  #"Oligodendrocyte_subclass_seurat"
fn="${d}/${input}"
label="subclass"
cellnum=5000
#out="${d}/${input}_${label}_5k"
co="t"


#########sh sb.sh -c ${comd}.sh -e $comd -j $ct -p fg -m 60 -t 0-2 -- $d $input $label $cellnum $co


ct="Mural_cell"
d="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/"
input=${ct}_subclass_rmRPE  #"Oligodendrocyte_subclass_seurat"
fn="${d}/${input}"
label="subclass"
cellnum=5000
#out="${d}/${input}_${label}_5k"
co="t"


#########sh sb.sh -c ${comd}.sh -e $comd -j $ct -p fg -m 60 -t 0-2 -- $d $input $label $cellnum $co

#Endothelial_cell_subclass_rmRPE
ct="Endothelial_cell"
d="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/"
input=${ct}_subclass_sb_seurat_rmRPE  #"Oligodendrocyte_subclass_seurat"

#input=${ct}_subclass_rmRPE  #"Oligodendrocyte_subclass_seurat"
fn="${d}/${input}"
label="subclass"
cellnum=5000
#out="${d}/${input}_${label}_5k"
co="t"


sh sb.sh -c ${comd}.sh -e $comd -j $ct -p fg -m 60 -t 0-2 -- $d $input $label $cellnum $co
