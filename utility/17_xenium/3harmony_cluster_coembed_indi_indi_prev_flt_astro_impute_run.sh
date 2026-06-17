#!/bin/sh
comd=HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_indi_prev_flt_astro_impute_astro
mkdir $comd


#rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_simple.rds
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_celltype2.obs.gz
rds="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple.rds"
label=astro
bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_w1"
#ds="BCM_22_0698_ONH_RNA"
ds="BCM_22_0047_ON_ONH_RNA_ds"
#ds=BCM_22_0784_ON_RNA
#cut=20
#f=0
cut=20
f=5
sc="celltype2"
for i in ONH_1 ONH_2 ONH_3 PP_1 PP_2
do
echo $i
sh sb.sh -c ${comd}.sh -m 50 -p fg  -t 0-20 -e $comd -j ${label}_${i} -- $rds $meta $bc $label $ds $cut $f $sc $i
done



