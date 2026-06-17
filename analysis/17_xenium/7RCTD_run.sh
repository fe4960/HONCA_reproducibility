#!/bin/sh
comd=HCA_ON/scripts/17_xenium/7RCTD
mkdir $comd
#rds="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_OLIGO_ref.combined_487_plus_4panel_cca.rds"
#rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_simple.rds
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_simple_sort.obs.gz

#label=astro #oligo
#ds=BCM_22_0784_ON_RNA
#ds="BCM_22_0698_ONH_RNA"
#label=oligo
#bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium"
#clu=subclass1




#rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_simple.rds
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_simple.obs.gz

#ds="BCM_22_0896_ON_RNA"
#label=oligo
#bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium"
#clu=subclass



#sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-20 -e $comd -j $label -- $rds $label $bc $meta $ds $clu


#rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_simple.rds
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_simple.obs.gz
rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_all_rawcount_5k_simple.rds
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_all_rawcount_5k_simple.obs.gz
ds="BCM_22_0896_ON_RNA"
label=major
bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium"
clu=majorclass

#comd=HCA_ON/scripts/17_xenium/7RCTD_major

sh sb.sh -c ${comd}.sh -m 50 -p s -t 0-20 -e $comd -j $label -- $rds $label $bc $meta $ds $clu
