#!/bin/sh
bam=/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/bam/${1}/${1}.bam
outdir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellranger/pub/GSE236566/fastq/${1}/
bam2fastq=/dfs3b/ruic20_lab/junw42/software/cellranger-7.1.0/lib/bin/bamtofastq
$bam2fastq --nthreads=10 $bam $outdir
