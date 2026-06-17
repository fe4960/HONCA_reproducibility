#!/bin/sh
comd=HCA_ON/scripts/16_general/1save_obs_general
mkdir $comd

#ct="Oligodendrocyte"
#h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/${ct}_subclass_seurat"

#sh sb.sh -m 70 -t 0-1 -p fg -c ${comd}.sh -e $comd -j ${ct}_subclass_seurat -- $h5ad


ct="Oligodendrocyte_precursor_cell"
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/${ct}_subclass_seurat"

#sh sb.sh -m 70 -t 0-1 -p fg -c ${comd}.sh -e $comd -j ${ct}_subclass_seurat -- $h5ad


ct="Astrocyte"
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/${ct}_subclass_new5_clean"

#sh sb.sh -m 70 -t 0-1 -p fg -c ${comd}.sh -e $comd -j ${ct}_subclass_new5_clean -- $h5ad 

ct="Fibroblast"

#h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/${ct}_subclass_clean_new"

h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat1"

#sh sb.sh -m 70 -t 0-1 -p fg -c ${comd}.sh -e $comd -j ${ct}_subclass_clean_new_rename_rmRPE_seurat1  -- $h5ad

#h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_1"
#h5ad=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_recluster

#sh sb.sh -m 50 -t 0-1 -p fg -c ${comd}.sh -e $comd -j ${ct}_subclass_clean_new_rename_rmRPE_seurat_v3_recluster  -- $h5ad

#sh sb.sh -m 70 -t 0-1 -p fg -c ${comd}.sh -e $comd -j ${ct}_subclass_clean_new_rename_rmRPE_seurat_v3_1  -- $h5ad



#sh sb.sh -m 70 -t 0-1 -p fg -c ${comd}.sh -e $comd -j ${ct}_subclass_clean_new  -- $h5ad


#ct="major"

#h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_all_rawcount_5k"

ct="major"
#h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_raw_normcount_only_5k"
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_raw_normcount_only_2k_celltype_rmUnk"
#sh sb.sh -m 50 -t 0-1 -p fg -c ${comd}.sh -e $comd -j ONONH_all_raw_normcount_only_2k_celltype_rmUnk  -- $h5ad

ct="major"
h5ad=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final/ONONH_all_rawcount_subset
#sh sb.sh -m 50 -t 0-1 -p fg -c ${comd}.sh -e $comd -j ONONH_all_raw_normcount_only_2k_celltype_rmUnk  -- $h5ad



ct=Astrocyte_subclass_new5_clean_5k
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_5k"

#sh sb.sh -m 50 -t 0-1 -p fg -c ${comd}.sh -e $comd -j $ct  -- $h5ad

ct=Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location1
h5ad=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location1

#ct=DRG_neuron
#h5ad=/dfs3b/ruic20_lab/junw42/DRG/data/scvi/raw/major/DRG_horse_ultima_merged_anno_full.clean.neuron_final

ct=CB_Fibro_new
h5ad=/dfs3b/ruic20_lab/junw42/HCA_TMCB/data/manuscript/MAGMA/cb_fibro_subtype_02_rawcount_anno_clean

sh sb.sh -m 20 -t 0-2 -p s -a ruic20_lab -c ${comd}.sh -e $comd -j $ct  -- $h5ad
