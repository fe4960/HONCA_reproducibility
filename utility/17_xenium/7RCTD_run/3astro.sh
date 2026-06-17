#!/bin/sh
#comd=HCA_ON/scripts/17_xenium/7RCTD
comd=HCA_ON/scripts/17_xenium/7RCTD_subclass_astro
mkdir $comd
#rds="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_OLIGO_ref.combined_487_plus_4panel_cca.rds"

rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple.rds
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple_sort_subclass3.obs.gz

#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_celltype2.obs.gz
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_celltype3.obs.gz

#####meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple_sort.obs.gz
label=astro
#label=Astrocyte # astro #oligo
#ds=BCM_22_0784_ON_RNA
ds="BCM_22_0698_ONH_RNA"
#label=oligo

bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_w1_8"

###bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_w1"
#clu=subclass3
#clu=celltype3
#clu=subclass1
clu=celltype2




#rds=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_simple.rds
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_simple.obs.gz

#ds="BCM_22_0698_ONH_RNA"
#label=fibro
#bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium"
#clu=subclass
#-n hpc3-14-12
#for i in {1..5}
#for i in {1..5}
for i in 6 7
do
#sh sb.sh -c ${comd}.sh -m 70 -p s -n hpc3-l18-01  -t 0-20 -e $comd -j $label -- $rds $label $bc $meta $ds $clu $i
sh sb.sh -c ${comd}.sh -m 70 -p g -a ruic20_lab_gpu   -t 0-20 -e $comd -j $label -- $rds $label $bc $meta $ds $clu $i

#sh sb.sh -c ${comd}.sh -m 70 -p s -n hpc3-14-12  -t 0-20 -e $comd -j $label -- $rds $label $bc $meta $ds $clu $i
done
