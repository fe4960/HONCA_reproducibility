#!/bin/sh
source ~/.condainit
conda activate scvi1.3

main="/dfs3b/ruic20_lab/junw42"

#ref="${main}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_simple.h5ad"
ref="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.h5ad"
query="${main}/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/major_list"

#query="${main}/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/oligo_xenium_h5ad_list"
syst="xenium"
outdir="${main}/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi"
label=major
cl=majorclass
mkdir $outdir
python ${main}/HCA_ON/scripts/17_xenium/3sysVI.py $ref $query $syst $outdir $label $cl
