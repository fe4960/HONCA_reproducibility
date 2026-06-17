#!/bin/sh

main="/dfs3b/ruic20_lab/junw42"
err=${main}/HCA_ON/scripts/5_refine_major/6refine_scvi/
mkdir -p $err
s=7
#res=0.3
res=0.4
file=${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general_cellclass
#while read line
for line in  Astrocyte #Microglia  #Astrocyte Oligodendrocyte Microglia   #Oligodendrocyte_precursor_cell   #AC   #RPE RGC Oligodendrocyte_precursor_cell MG Endothelial_cell BC
do
#rml="rm_cluster_rs"
#rml="rm_cluster_rs_1"
#rml="rm_cluster_rs_0.4"

flavor="seurat_v3"
hvg=2000
epoch="none" 
gl="zinb" 
span=1
l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_sb
sbatch -p gpu --account=ruic20_lab_gpu --time=0-6 --mem=100GB --gres=gpu:1 --error=${err}/${l}.err --output=${err}/${l}.out  ${main}/HCA_ON/scripts/5_refine_major/6refine_scvi_general.sh $line $s $res $rml $flavor $hvg $epoch $gl $span

done 
#< $file 

