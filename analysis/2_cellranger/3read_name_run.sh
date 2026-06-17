#!/bin/sh
err=HCA_ON/scripts/2_cellranger/3read_name
mkdir $err
file=/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE267301_list
while read line
do
for i in R1 R2
do
sbatch --mem=20GB --time=0-10 -p standard --account=ruic20_lab --output=${err}/${line}.out --error=${err}/${line}.err ${err}.sh /dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE267301/${line} $i
done
done < $file
