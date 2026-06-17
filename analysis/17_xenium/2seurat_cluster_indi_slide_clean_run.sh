#!/bin/sh

#HCA_ON/scripts/17_xenium/2seurat_cluster_indi.sh 

comd=HCA_ON/scripts/17_xenium/2seurat_cluster_indi_slide_clean
mkdir $comd
bc=/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH/ONONH_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_oligo
label=major
sh sb.sh -c ${comd}.sh -m 100 -p s -t 0-2 -e $comd -j $label -- $bc $label

rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_simple.rds
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_simple.obs.gz
bc=/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH/ONONH_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_opc
label=opc

#sh sb.sh -c ${comd}.sh -m 100 -p fg -t 0-2 -e $comd -j $label -- $rds $meta $bc $label

bc=/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH/ONONH_xenium_umap_p5_cutoff20_facet_ref_xenium_no_inte_fibro
label=fibro
#sh sb.sh -c ${comd}.sh -m 100 -p fg -t 0-2 -e $comd -j $label -- $bc $label

