#!/bin/sh
source ~/.condainit
conda activate scvi

ct="Oligodendrocyte_precursor_cell"
fn="subclass_seurat_cycling_rmRPE"
#python HCA_ON/scripts/5_refine_major/7celltype_class/1final_opc_clean_seurat_top_ranking_gene.py $ct $fn

fn=subclass_seurat_rmRPE
ct="Oligodendrocyte"
#python HCA_ON/scripts/5_refine_major/7celltype_class/1final_opc_clean_seurat_top_ranking_gene.py $ct $fn


fn=subclass_new5_clean
ct="Astrocyte"
#python HCA_ON/scripts/5_refine_major/7celltype_class/1final_opc_clean_seurat_top_ranking_gene.py $ct $fn


fn=subclass_clean_new_rename_rmRPE_seurat_v3_rename_location
ct="Fibroblast"
python HCA_ON/scripts/5_refine_major/7celltype_class/1final_opc_clean_seurat_top_ranking_gene.py $ct $fn
