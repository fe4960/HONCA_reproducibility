#!/bin/sh

meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/fastq_paths_full1

while read line
do
IFS="	" ARR=($line)
for file in $( ls ${ARR[1]}/*R1_001.fastq.gz )
do
sh HCA_ON/scripts/1_SRAdownload/4read_len.sh $file
done 
done < $meta

