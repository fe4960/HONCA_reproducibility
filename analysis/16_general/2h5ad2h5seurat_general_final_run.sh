#!/bin/sh
comd="HCA_ON/scripts/16_general/2h5ad2h5seurat_general_final"
mkdir $comd

s=Oligodendrocyte_subclass_seurat_simple

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"

#sh sb.sh -c ${comd}.sh -m 70 -t 0-1 -p fg -e $comd -j $s -- $s $dir


s=Oligodendrocyte_precursor_cell_subclass_seurat_simple

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"

#sh sb.sh -c ${comd}.sh -m 70 -t 0-1 -p fg -e $comd -j $s -- $s $dir

s=Astrocyte_subclass_new5_clean_simple

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"

#sh sb.sh -c ${comd}.sh -m 70 -t 0-1 -p fg -e $comd -j $s -- $s $dir

s=ONONH_all_rawcount_subset
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/
#sh sb.sh -c ${comd}.sh -m 70 -t 0-1 -p s -e $comd -j $s -- $s $dir

#s=Fibroblast_subclass_clean_new_simple
s=Fibroblast_subclass_clean_new_rename_rmRPE_seurat1_simple
dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
#sh sb.sh -c ${comd}.sh -m 70 -t 0-1 -p s -e $comd -j $s -- $s $dir

#s=Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_1_simple
s=Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_recluster_simple
dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
#sh sb.sh -c ${comd}.sh -m 50 -t 0-1 -p fg -e $comd -j $s -- $s $dir


s=ONONH_all_rawcount_5k_simple
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/
#sh sb.sh -c ${comd}.sh -m 100 -t 0-2 -p s -e $comd -j $s -- $s $dir


s=ONONH_all_rawcount_subset
#dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/
#dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/
#sh sb.sh -c ${comd}.sh -m 150 -t 0-2 -p fg -e $comd -j $s -- $s $dir
#/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_raw_normcount_only_2k_celltype_rmUnk
#s=ONONH_all_raw_normcount_only_2k_celltype_rmUnk_simple
#dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/
s=ONONH_all_rawcount_subset
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/





#dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/
#s=Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location1_simple

#dir=/dfs3b/ruic20_lab/junw42/DRG/data/scvi/raw/major/
#s=DRG_horse_ultima_merged_anno_full.clean.neuron_final_simple
dir=/dfs3b/ruic20_lab/junw42/HCA_TMCB/data/manuscript/MAGMA/
s=cb_fibro_subtype_02_rawcount_anno_clean_simple
sh sb.sh -c ${comd}.sh -m 20 -t 0-2 -p s -a ruic20_lab -e $comd -j $s -- $s $dir
