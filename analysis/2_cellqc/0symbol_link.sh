#!/bin/sh

#dir=/storage/chenlab/Data/HCA/human_ancestry_ret/RNA/cellranger
dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/"

mkdir $dir

file="/pub/junw42/proj/HCA_ON/scripts/1cellqc/ON_chen_crpath.txt_cellranger"



while read line
do

IFS="/" read -ra ARR <<< $line

#new_dir=${dir}/${ARR[-2]}/${ARR[-1]}
new_dir=${dir}/${ARR[-3]}

mkdir -p  ${new_dir}

for file in analysis.tar.gz  filtered_feature_bc_matrix.h5  metrics_summary.csv  possorted_genome_bam.bam  raw_feature_bc_matrix.h5  web_summary.html  cloupe.cloupe    filtered_feature_bc_matrix.tar.gz  molecule_info.h5  possorted_genome_bam.bam.bai  raw_feature_bc_matrix.tar.gz

do

ln -s ${line}/${file}  ${new_dir}/${file}

done

done < $file

#done


