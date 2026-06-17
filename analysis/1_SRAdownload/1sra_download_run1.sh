#!/bin/sh
err=HCA_ON/scripts/1_SRAdownload/1sra_download

#err=HCA_ON/scripts/SRA_download/1sra_download
#output=/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE236566/

#output=/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE267301/

output=/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/PFC_SOX9
file=/dfs3b/ruic20_lab/junw42/HCA_ON/data/fastq/pub/GSE167494/SRA_list_normal_sox9
mkdir -p $output
mkdir -p $err
while read id
#for id in SRR29007591 SRR29007592 SRR29007593 SRR29007594 SRR29007595 SRR30562939 SRR30562940 SRR30562941 SRR30562942 SRR30562943 SRR30562944 SRR30562945
#for id in SRR25147535 SRR25147536 SRR25147537 SRR25147538 SRR25147539 SRR25147540 SRR25147541 SRR25147542 SRR25147543 SRR25147544 SRR25147545 SRR25147546 SRR25147547 SRR25147548 SRR25147549 SRR25147550 SRR25147551 SRR25147552 SRR25147553 SRR25147554 SRR25147555
do
sbatch -p standard --mem=5GB --time=4-0 --output=${err}/${id}.out --error=${err}/${id}.err --account=ruic20_lab ${err}.sh $id $output

	#sbatch -p standard --mem=5GB --time=4-0 --output=${err}/${id}.out --error=${err}/${id}.err --account=ruic20_lab HCA_ON/scripts/SRA_download/1sra_download.sh $id $output
done < $file
