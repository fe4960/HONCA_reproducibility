#!/bin/sh
err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/2_cellranger/2cellranger
mkdir $err
#meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/SRA_list_normal"

meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_AD_list_sam3"

#meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_AD_list_sam1"
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_AD
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/2_cellranger/2cellranger_GSE167494.sh
while read line
#for line in 
do
id=$line #${ARR[0]}
sbatch --ntasks-per-node=16 --mem=95GB --time=1-0 --error=${err}/${id}.err --output=${err}/${id}.out --account=ruic20_lab -p free-gpu --account=ruic20_lab_gpu  $comd ${id} $dir 

done < $meta




meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_list_sam"
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/2_cellranger/2cellranger_GSE167494.sh
while read line
#for line in 
do
id=$line #${ARR[0]}
##########3sbatch --ntasks-per-node=16 --mem=80GB --time=1-0 --error=${err}/${id}.err --output=${err}/${id}.out --account=ruic20_lab -p standard  $comd ${id} $dir 

done < $meta

#err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/2_cellranger/2cellranger
#mkdir $err
#meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/SRA_list_normal_sox9"
meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_SOX9_list_sam"
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_SOX9
while read line
#for line in SRR13782284 SRR13782283 SRR13782276 SRR13782277 
#for line in SRR13782282
do
id=$line #${ARR[0]}
############sbatch --ntasks-per-node=16 --mem=80GB --time=1-0 --error=${err}/${id}.err --output=${err}/${id}.out --account=ruic20_lab -p standard  $comd ${id} $dir 

done < $meta
