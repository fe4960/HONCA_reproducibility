#!/bin/sh
#file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pri/cellranger_list"
error="/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/6_velocyto/1run_velocyto"
file="HCA_ON/scripts/6_velocyto/1run_velocyto_missing"
mkdir $error

#while read line
#do
#data="/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pri/${line}"
#sbatch -p standard --account=ruic20_lab --mem=60GB --time=2-0 --output=${error}/${line}.out --error=${error}/${line}.err  /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/6_velocyto/1run_velocyto.sh $data
#done < $file

file1=/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE236566/cellranger/cellranger_list_full

#while read line
#do
#arr=($line)
#echo "${arr[0]}\n"
#data=/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE236566/cellranger/${arr[0]}/${arr[0]}/outs
#sbatch -p standard --account=ruic20_lab --mem=60GB --time=2-0 --output=${error}/${arr[0]}.out --error=${error}/${arr[0]}.err  /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/6_velocyto/1run_velocyto.sh $data
#done < $file1

for line in GSM7553439
do
data=/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE236566/cellranger/${line}/${line}/outs
sbatch -p standard --account=ruic20_lab --mem=60GB --time=2-0 --output=${error}/${line}.out --error=${error}/${line}.err  /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/6_velocyto/1run_velocyto.sh $data
done 
