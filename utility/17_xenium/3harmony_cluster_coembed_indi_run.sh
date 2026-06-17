#!/bin/sh
comd=HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi
mkdir $comd
#rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_simple.rds
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_simple.obs.gz
rds="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_OLIGO_ref.combined_487_plus_4panel_cca.rds"

#rds="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_OLIGO_ref.combined_5k_plus_panel.rds"
#bc=/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH/ONONH_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_oligo
label=oligo
bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium"
#ds=9000
#ds=BCM_22_0784_ON_RNA
sh sb.sh -c ${comd}.sh -m 50 -p fg -t 0-7 -e $comd -j $label -- $rds $label $bc

rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_simple.rds
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_simple.obs.gz
bc=/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH/ONONH_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_opc
label=opc

#sh sb.sh -c ${comd}.sh -m 100 -p fg -t 0-2 -e $comd -j $label -- $rds $meta $bc $label

rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_simple.rds
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_simple.obs.gz
bc=/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH/ONONH_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_fibro
label=fibro
#ds=15000
ds=BCM_22_0784_ON_RNA
#sh sb.sh -c ${comd}.sh -m 100 -p fg -t 0-2 -e $comd -j $label -- $rds $meta $bc $label $ds
