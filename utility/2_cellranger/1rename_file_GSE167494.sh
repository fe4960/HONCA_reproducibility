#!/bin/sh

#dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_SOX9"
#file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_SOX9_list"


#dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC"
#file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_list"

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_AD"
file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_list_AD"


pre="none"
while read line
do
IFS="	" read -ra arr <<< "$line"
if [ $pre != ${arr[1]} ]; then
	n=1
	pre=${arr[1]}
else
	n=$(($n+1))
fi
#echo "${dir}/${arr[1]}_S1_L00${n}_R1_001.fastq.gz"
#mv ${dir}/${arr[0]}_S1_L00${n}_R1_001.fastq.gz ${dir}/${arr[1]}_S1_L00${n}_R1_001.fastq.gz
#mv ${dir}/${arr[0]}_S1_L00${n}_R2_001.fastq.gz ${dir}/${arr[1]}_S1_L00${n}_R2_001.fastq.gz

mv ${dir}/${arr[0]}_S1_L001_R1_001.fastq.gz ${dir}/${arr[1]}_S1_L00${n}_R1_001.fastq.gz
mv ${dir}/${arr[0]}_S1_L001_R2_001.fastq.gz ${dir}/${arr[1]}_S1_L00${n}_R2_001.fastq.gz


done < $file

