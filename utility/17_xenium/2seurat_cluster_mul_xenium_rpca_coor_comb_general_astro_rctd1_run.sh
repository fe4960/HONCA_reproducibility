#!/bin/sh

comd=HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general_astro_rctd1

mkdir $comd

label=fibro
cut=20
f=5
clu=xenium_clusters_p5

#sh sb.sh -c ${comd}.sh -e $comd -j ${label}_${cut}_${f}_${clu} -p fg -m 30 -t 0-1 -- $label $cut $f $clu

label=astro
cut=20
f=5
clu="rctd1"
sc="subclass1"
#clu=xenium_clusters_p5

#sh sb.sh -c ${comd}.sh -e $comd -j ${label}_${cut}_${f}_${clu} -p fg -m 30 -t 0-1 -- $label $cut $f $clu $sc


label=fibro
cut=20
f=3
clu=xenium_clusters_p5

#sh sb.sh -c ${comd}.sh -e $comd -j ${label}_${cut}_${f}_${clu} -p fg -m 30 -t 0-1 -- $label $cut $f $clu

label=astro
cut=20
f=5
clu="rctd1" #xenium_clusters_p5
sc="subclass1"
idx="015"
#sh sb.sh -c ${comd}.sh -e $comd -j ${label}_${cut}_${f}_${clu}_${idx} -p fg -m 30 -t 0-1 -- $label $cut $f $clu $sc $idx

sc="subclass"
idx="025"
#sh sb.sh -c ${comd}.sh -e $comd -j ${label}_${cut}_${f}_${clu}_${idx} -p fg -m 30 -t 0-1 -- $label $cut $f $clu $sc $idx



label=oligo
cut=20
f=5
clu="rctd1"
sc="subclass"
idx="01"
#clu=xenium_clusters_p5

sh sb.sh -c ${comd}.sh -e $comd -j ${sc}_${label}_${cut}_${f}_${clu} -p s -m 30 -t 0-1   -- $label $cut $f $clu $sc $idx



label=fibro
cut=20
f=5
clu="rctd1"
sc="subclass1"
#idx="025"
#idx="n" #"n"
idx="015"
for label1 in seurat_v3 seurat_v3_recluster
do
#sh sb.sh -c ${comd}.sh -e $comd -j ${sc}_${label}_${cut}_${f}_${clu}_${label1} -p s  -m 30 -t 0-1 -- $label $cut $f $clu $sc $idx ${label1}
done

label=astro
cut=20
f=5
clu="rctd1"
#sc="subclass1"
sc="subclass3"
#idx="025"
idx="n"
#sh sb.sh -c ${comd}.sh -e $comd -j ${sc}_${label}_${cut}_${f}_${clu} -p fg -m 30 -t 0-1  -- $label $cut $f $clu $sc $idx


