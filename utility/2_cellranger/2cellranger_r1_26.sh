#!/bin/sh
cellranger=/dfs3b/ruic20_lab/junw42/software/cellranger-7.1.0/bin/cellranger
ref=/dfs3b/ruic20_lab/junw42/reference/cellranger/refdata-gex-GRCh38-2020-A
sample="bamtofastq"
id=$1
dir=$2

outdir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/cellranger/${id}
rm -r ${outdir}
mkdir -p $outdir
cd $outdir

$cellranger count --include-introns=true \
   --id=$id \
   --fastqs=$dir \
   --sample=$sample \
   --transcriptome=$ref \
   --localcores=16 \
   --localmem=70 \
   --r1-length=26

