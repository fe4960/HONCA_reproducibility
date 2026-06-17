#!/bin/sh

#dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE267301/"
#dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_SOX9"
dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_AD"

#dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC"
#file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE267301_list"
#file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/SRA_list_normal"
#file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/SRA_list_normal_sox9"
file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/SRA_list_AD"
while read line
do


mv ${dir}/${line}_2.fastq.gz  ${dir}/${line}_S1_L001_R1_001.fastq.gz
mv ${dir}/${line}_3.fastq.gz  ${dir}/${line}_S1_L001_R2_001.fastq.gz

	
#mv ${dir}/${line}_3.fastq.gz  ${dir}/${line}_S1_L001_R1_001.fastq.gz
#mv ${dir}/${line}_4.fastq.gz  ${dir}/${line}_S1_L001_R2_001.fastq.gz
done < $file

