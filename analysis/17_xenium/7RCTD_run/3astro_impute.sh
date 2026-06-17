#!/bin/sh
#comd=HCA_ON/scripts/17_xenium/7RCTD
comd=HCA_ON/scripts/17_xenium/7RCTD_subclass_astro_impute
mkdir $comd
#rds="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_OLIGO_ref.combined_487_plus_4panel_cca.rds"
rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple.rds
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple_sort_subclass3.obs.gz
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_celltype2.obs.gz
#####meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple_sort.obs.gz
label=astro
ds="BCM_22_0698_ONH_RNA"
bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_w1"

for clu in celltype1 celltype celltype2
do
for i in ONH_1 ONH_2 ONH_3 PP_1 PP_2
do
sh sb.sh -c ${comd}.sh -m 70 -p fg -t 0-20 -e $comd -j $label -- $rds $label $bc $meta $ds $clu $i
done
done
