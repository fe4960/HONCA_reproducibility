#!/bin/sh
source ~/.condainit
conda activate scvi

main="/dfs3b/ruic20_lab/junw42"

#ref="${main}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_simple.h5ad"
#ref="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.h5ad"
ref="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount.h5ad"
query="${main}/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/fibro_list"

#query="${main}/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/oligo_xenium_h5ad_list"
syst="xenium"
outdir="${main}/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi"
#label=major

label=fibro
cl=celltype1 #majorclass

sn_type="Fibroblast"
sn_col="majorclass"
st_type="fibro"
#dt="${main}/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi/fibro_subclass_subclass1"
mkdir $outdir
python ${main}/HCA_ON/scripts/17_xenium/10scVI_harmony_fibro.py $ref $query $syst $outdir $label ${sn_type} ${sn_col} 
