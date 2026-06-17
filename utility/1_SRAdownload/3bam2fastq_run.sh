#!/bin/sh
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/1_SRAdownload/3bam2fastq.sh
err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/1_SRAdownload/3bam2fastq
mkdir $err
#meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/SraRunTable.txt_ONONH"
meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/SraRunTable.txt_ONONH"

while read line
do
ARR=($line)
id=${ARR[3]}
#sample=${ARR[0]}
#bam=${ARR[5]}
sbatch -p standard --mem=40G --time=4-0 --account=ruic20_lab --error=${err}/${id}.err --output=${err}/${id}.out $comd $id 
done < $meta
