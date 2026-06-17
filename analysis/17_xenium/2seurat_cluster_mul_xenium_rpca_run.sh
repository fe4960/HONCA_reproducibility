#!/bin/sh
###cut=50
###f=10
#cut=20
#f=0
cut=20
f=5
#f=25
comd=HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca
mkdir $comd
#for l in oligo #fibro astro oligo
for l in astro fibro oligo #fibro #astro fibro #astro oligo
do
sh sb.sh -c ${comd}.sh -m 25 -t 0-2 -p fg -e $comd -j ${l}_${cut} -- $cut $l $f
done
