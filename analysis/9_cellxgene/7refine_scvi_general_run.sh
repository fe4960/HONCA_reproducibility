#!/bin/sh

main="/dfs3b/ruic20_lab/junw42"
err=${main}/HCA_ON/scripts/5_refine_major/6refine_scvi/
mkdir -p $err
s=7
#s=12345
#res=0.2
#res=0.6
res=0.4
#file=${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general_cellclass
#while read line
#for line in    Oligodendrocyte   #Endothelial_cell  #HC Cone #RPE #AC   #Oligodendrocyte  #Astrocyte #Microglia   #Oligodendrocyte_precursor_cell   #AC   #RPE RGC Oligodendrocyte_precursor_cell MG Endothelial_cell BC
for line in Astrocyte
#for line in Oligodendrocyte
do
#rml="rm_cl_1_8" #"none"
rml="none"
for kp in  none #kp_cl_0   #kp_cl_2 kp_cl_0  kp_cl_16 kp_cl_3 kp_cl_4 
do
#kp="kp_cl_2"
flavor="seurat_v3"
#flavor="seurat"
hvg=10000
#hvg=2000
epoch="none" 
gl="zinb" 
span=1
#l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_${kp}_ret_onh
#l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_${kp}_rm_non_spe
#l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_${kp}_${rml}_on
l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_${kp}_${rml}_rm_clu1_8

sbatch -p free-gpu --account=ruic20_lab_gpu --time=0-2 --mem=50GB --gres=gpu:1 --error=${err}/${l}.err --output=${err}/${l}.out  ${main}/HCA_ON/scripts/9_cellxgene/7refine_scvi_general.sh $line $s $res $rml $flavor $hvg $epoch $gl $span $kp
done

done 
#< $file 

