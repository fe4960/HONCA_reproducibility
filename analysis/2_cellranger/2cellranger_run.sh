#!/bin/sh
err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/2_cellranger/2cellranger
mkdir $err
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/fastq_paths_full
meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/fastq_paths_full1
#grep -P "GSM7553434" /dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/fastq_paths_full > /dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/fastq_paths_full2
#meta=/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/fastq_paths_full2


while read line
do
IFS="	" ARR=($line)
id=${ARR[0]}
dir=${ARR[1]}
echo "$id $dir"
#sbatch --ntasks-per-node=16 --mem=70GB --time=5-0 --error=${err}/${id}.err --output=${err}/${id}.out --account=ruic20_lab -p standard  /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/2_cellranger/2cellranger.sh $id $dir 
sbatch --ntasks-per-node=16 --mem=70GB --time=5-0 --error=${err}/${id}.err --output=${err}/${id}.out --account=ruic20_lab -p standard  /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/2_cellranger/2cellranger_r1_26.sh $id $dir 

done < $meta
