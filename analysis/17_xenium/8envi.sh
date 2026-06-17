#!/bin/sh

source ~/.condainit
conda activate envi

sc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.h5ad"
output="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/"
target_tissue="ONONH_comb"
label="major"
cl="majorclass"

query_list="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/major_list"

python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/8envi.py $sc $output ${target_tissue} $label $cl ${query_list}
