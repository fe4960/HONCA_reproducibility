#!/bin/sh

#module load awscli/2.11.21

#indir="https://sra-pub-src-2.s3.amazonaws.com"
indir="https://sra-pub-src-1.s3.amazonaws.com"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/bam"

mkdir $outdir


id=$1 
sample=$2 
bam=$3 
dir=${outdir}/${id}
mkdir $dir
link=${indir}/${sample}/${bam}
#echo "$dir"
#echo "$link"
echo "curl -o  ${dir}/${id}.bam $link"
curl -o  ${dir}/${id}.bam $link



