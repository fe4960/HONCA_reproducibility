#!/bin/sh
#main_dir=/storage/chenlab/Users/junwang/monkey/scripts/database/
#main_data_dir=/storage/chenlab/Users/junwang/monkey/data/database
main_dir=/dfs3b/ruic20_lab/junw42
#main_data_dir=
#export LD_LIBRARY_PATH=/storage/chen/Software/lib/usr/lib64:$LD_LIBRARY_PATH

#dir="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final/PeakCalls_bed"

dir=$1 #"/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/PeakCalls/"

file=$2 #"/storage/chenlab/Users/junwang/human_meta/data/cell_list"

#/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/DAR/

while read cell

do

############tail -n+2 ${dir}/${cell}-reproduciblePeaks.bed | sort -u | awk '{print $1":"$2"-"$3}' > ${dir}/${cell}_all_peak_hg38
#tail -n+2 ${dir}/${cell}_DAR_peak_FDR005 | sort -u | awk '{print $2":"$3"-"$4}' > ${dir}/${cell}_DAR_peak_FDR005_hg38

tail -n+2 ${dir}/${cell}_DAR_peak_FDR005 | sort -u | awk '{print $2":"$4"-"$5}' > ${dir}/${cell}_DAR_peak_FDR005_hg38


${main_dir}/software/liftOver -positions  ${dir}/${cell}_DAR_peak_FDR005_hg38  ${main_dir}/reference/ucsc_chain/hg38ToHg19.over.chain.gz  ${dir}/${cell}_DAR_peak_FDR005_hg38_hg19 ${dir}/${cell}_DAR_peak_FDR005_hg38_hg19_unmap


done < $file
