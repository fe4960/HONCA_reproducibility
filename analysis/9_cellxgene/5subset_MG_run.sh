#!/bin/sh

main="/dfs3b/ruic20_lab/junw42"
err=${main}/HCA_ON/scripts/5_refine_major/6refine_scvi/
mkdir -p $err
s=7
#res=0.2
res=1
#res=0.4
file=${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general_cellclass
#while read line
#for line in    Oligodendrocyte   #Endothelial_cell  #HC Cone #RPE #AC   #Oligodendrocyte  #Astrocyte #Microglia   #Oligodendrocyte_precursor_cell   #AC   #RPE RGC Oligodendrocyte_precursor_cell MG Endothelial_cell BC
for line in major #Fibroblast #Astrocyte
do
#rml="rm_cluster_rs_1_rm"
#rml="rm_cluster_rs_0.3"
#rml="rm_cluster_rs"
rml="none"
#rml="rm_cluster_rs_4_6_13"

#rml="rm_cluster_rs_46"
#rml="rm_cluster_rs_4"
#rml="rm_cluster_rs_1"
#rml="rm_cluster_rs_0.4"
#flavor="seurat"
#kp="kp_cl_06"
for kp in  kp_cl_16   #kp_cl_2 kp_cl_0  kp_cl_16 kp_cl_3 kp_cl_4 
do
#kp="kp_cl_2"
flavor="seurat_v3"
#flavor="seurat"
#hvg=2000
hvg=10000
#hvg=2000
epoch="none" 
gl="zinb" 
span=1
#l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_rm_rm
#l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_rm_sb
#l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_rm_nonsb
l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_${kp}_sb

#l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_rm_sb
sbatch -p gpu --account=ruic20_lab_gpu --time=0-2 --mem=150GB --gres=gpu:1 --error=${err}/${l}.err --output=${err}/${l}.out  ${main}/HCA_ON/scripts/9_cellxgene/5subset_MG.sh $line $s $res $rml $flavor $hvg $epoch $gl $span $kp
done

done 
#< $file 

