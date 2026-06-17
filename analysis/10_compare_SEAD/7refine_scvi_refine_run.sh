#!/bin/sh

main="/dfs3b/ruic20_lab/junw42"
err=${main}/HCA_ON/scripts/10_compare_SEAD/7refine_scvi_general
mkdir -p $err
s=7
res=0.6
for line in Oligodendrocyte #Astrocyte  Oligodendrocyte_precursor_cell
do
rml="none"
#rml="rm_cluster_rs_3"
#rml="rm_cluster_rs_3_4_5"
for kp in  none   #kp_cl_2 kp_cl_0  kp_cl_16 kp_cl_3 kp_cl_4 
do
#flavor="seurat_v3"
flavor="seurat"
hvg=10000
epoch="none" 
gl="zinb" 
span=1
l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_${kp}
sbatch -p gpu --account=ruic20_lab_gpu    --time=0-2 --mem=100GB --gres=gpu:1 --error=${err}/${l}.err --output=${err}/${l}.out  ${main}/HCA_ON/scripts/10_compare_SEAD/7refine_scvi_general.sh $line $s $res $rml $flavor $hvg $epoch $gl $span $kp
done

done 

