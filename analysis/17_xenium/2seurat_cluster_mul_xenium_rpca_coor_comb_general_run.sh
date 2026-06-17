#!/bin/sh

comd=HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general

mkdir $comd

label=fibro
cut=20
f=5
#clu=xenium_clusters_p5
clu="harmony_anno"
sh sb.sh -c ${comd}.sh -e $comd -j ${label}_${cut}_${f}_${clu} -p fg -m 30 -t 0-1 -- $label $cut $f $clu

label=Astro
cut=20
f=5
clu=xenium_clusters_p5

#sh sb.sh -c ${comd}.sh -e $comd -j ${label}_${cut}_${f}_${clu} -p fg -m 30 -t 0-1 -- $label $cut $f $clu


label=fibro
cut=20
f=3
clu=xenium_clusters_p5

#sh sb.sh -c ${comd}.sh -e $comd -j ${label}_${cut}_${f}_${clu} -p fg -m 30 -t 0-1 -- $label $cut $f $clu

label=Astro
cut=20
f=3
clu=xenium_clusters_p5

#sh sb.sh -c ${comd}.sh -e $comd -j ${label}_${cut}_${f}_${clu} -p fg -m 30 -t 0-1 -- $label $cut $f $clu

