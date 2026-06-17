#!/bin/sh
file=/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/fastq_paths
rm ${file}_full
dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/fastq"
while read line
do
declare -a array=()
IFS="	" ARRA=($line)
IFS="," ARR=(${ARRA[1]})
for i in ${!ARR[@]}
do
array[$i]=${dir}/${ARRA[0]}/${ARR[$i]}
done
joined_path=$(IFS=","; echo "${array[*]}")
echo "${ARRA[0]}	${joined_path}" >> ${file}_full
done < $file
