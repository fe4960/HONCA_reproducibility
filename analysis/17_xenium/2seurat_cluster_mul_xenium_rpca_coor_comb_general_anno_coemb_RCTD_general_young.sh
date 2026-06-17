#!/bin/sh
source ~/.condainit
conda activate spatial1

#comd="HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general_anno_coemb_RCTD_general_young"
comd="HCA_ON/scripts/17_xenium/2seurat_cluster_mul_xenium_rpca_coor_comb_general_anno_coemb_RCTD_general_young1"

mkdir $comd

#label=oligo
label=astro
#label=fibro
cut=20
fc=5
#sc=subclass
sc=celltype2
#sc=celltype3
#sc="subclass1"
#sc="subclass1_seurat_final_6"
#sc=subclass3
id="015"
#id="n"
#id="020"
#id="020_ONHON_ON"
#id="01"
#id="025_default"
R --vanilla --args $label $cut $fc $sc $id < ${comd}.R > ${comd}/${label}_${cut}_${fc}_${sc}_${id}.out  2> ${comd}/${label}_${cut}_${fc}_${sc}_${id}.err
#Rscript --vanilla ${comd}.R $label $cut $fc $sc $id > ${comd}/${label}_${cut}_${fc}_${sc}_${id}.out  2> ${comd}/${label}_${cut}_${fc}_${sc}_${id}.err
