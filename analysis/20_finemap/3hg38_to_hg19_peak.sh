#!/bin/sh
##########main_dir=/storage/chenlab/Users/junwang/monkey/scripts/database/
##########main_data_dir=/storage/chenlab/Users/junwang/monkey/data/database
##########export LD_LIBRARY_PATH=/storage/chen/Software/lib/usr/lib64:$LD_LIBRARY_PATH

#dir="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final/PeakCalls_bed"

#dir="/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/" #"/storage/chenlab/Users/junwang/human_meta/data/proj4_clean1_final_major_full_max/PeakCalls/"
main_dir="/dfs3b/ruic20_lab/junw42"

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/"


#file=${dir}/peakset_clusters_major_full_peak_anno 
file=${dir}/peakset_clusters_major_full_peak_anno
tail -n+2 ${file} |  awk '{print $NF}' | sort -u  > ${file}_hg38
software/liftOver -positions  ${file}_hg38  ${main_dir}/reference/ucsc_chain/hg38ToHg19.over.chain.gz  ${file}_hg38_hg19 ${file}_hg38_hg19_unmap


perl human_meta/scripts/finemap/4match_hg38_hg19_coor.pl ${file}_hg38 ${file}_hg38_hg19 ${file}_hg38_hg19_unmap


file=${dir}/g2p_cA_final_major_full_max_LCRE_uniq_clean  #g2p_cA_final_major_full_max_LCRE_uniq_clean

tail -n+2 ${file} | awk '{print $NF}' | sort -u > ${file}_hg38

software/liftOver -positions  ${file}_hg38  ${main_dir}/reference/ucsc_chain/hg38ToHg19.over.chain.gz  ${file}_hg38_hg19 ${file}_hg38_hg19_unmap

perl human_meta/scripts/finemap/4match_hg38_hg19_coor.pl  ${file}_hg38 ${file}_hg38_hg19 ${file}_hg38_hg19_unmap


#done < $file
