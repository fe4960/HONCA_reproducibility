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


comd=HCA_ON/scripts/17_xenium/7RCTD_subclass_fibro

#rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_simple.rds
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_simple_sort.obs.gz

#rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_simple.rds
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_simple.obs.gz

rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_recluster_simple.rds
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_recluster_simple.obs.gz


ds="BCM_22_0698_ONH_RNA"
#label=Fibroblast
label=fibro
bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_w1_6"
clu=subclass1
label1="seurat_v3_recluster"


for i in {1..5}
do
sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-20 -e $comd -j ${label}_${label1} -- $rds $label $bc $meta $ds $clu $i $label1
done


#rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat1_simple.rds
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat1_simple.obs.gz

rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_simple.rds
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_simple.obs.gz


ds="BCM_22_0698_ONH_RNA"
#label=Fibroblast
label=fibro
bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_w1_6"
clu=subclass1
#label1="seurat1"
label1="seurat_v3"

for i in {1..5}
do
sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-20 -e $comd -j $label -- $rds $label $bc $meta $ds $clu $i $label1
done




#sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-20 -e $comd -j $label -- $rds $label $bc $meta $ds $clu

#comd=HCA_ON/scripts/17_xenium/7RCTD_subclass
#sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-20 -e $comd -j $label -- $rds $label $bc $meta $ds $clu
