#!/bin/sh
comd="/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/16_general/3seurat2h5ad_list"
mkdir $comd

#dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH/" 

####for s in ONONH_xenium_merge_cutoff20_oligo_1  # ONONH_xenium_merge_cutoff20_oligo_2
#####do
#s="ONONH_comb_xenium_data_list_comb_cutoff20"

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/"

for i in {1..5}
do
#echo "oligo"
s="ONONH_comb_xenium_merge_cutoff20_major_${i}"
sh sb.sh -c ${comd}.sh -m 70 -t 0-1 -p fg -e $comd -j $s -- $s $dir
done

for s in ONONH_xenium_merge_cutoff20_fibro_1  ONONH_xenium_merge_cutoff20_fibro_2
do
echo "fibro_1"
####sh sb.sh -c ${comd}.sh -m 70 -t 0-1 -p fg -e $comd -j $s -- $s $dir
done

##s=Astrocyte_subclass_new4_clean_simple

##dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"

#sh sb.sh -c ${comd}.sh -m 70 -t 0-1 -p fg -e $comd -j $s -- $s $dir

##s=ONONH_all_rawcount_subset
##dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/
#sh sb.sh -c ${comd}.sh -m 70 -t 0-1 -p s -e $comd -j $s -- $s $dir

##s=Fibroblast_subclass_clean_new_simple
##dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
##sh sb.sh -c ${comd}.sh -m 70 -t 0-1 -p s -e $comd -j $s -- $s $dir
