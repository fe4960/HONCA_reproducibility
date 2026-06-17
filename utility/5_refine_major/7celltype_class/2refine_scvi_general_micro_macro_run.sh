#!/bin/sh

main="/dfs3b/ruic20_lab/junw42"
err=${main}/HCA_ON/scripts/5_refine_major/6refine_scvi/
mkdir -p $err
s=7
res=0.6
file=${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general_cellclass
for line in major
do
rml="none"
for kp in  none   #kp_cl_2 kp_cl_0  kp_cl_16 kp_cl_3 kp_cl_4 
do
flavor="seurat_v3"
hvg=10000
epoch="none" 
gl="zinb" 
span=1
l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}_${kp}

sbatch -p gpu --account=ruic20_lab_gpu --time=0-2 --mem=30GB --gres=gpu:1 --error=${err}/${l}.err --output=${err}/${l}.out  ${main}/HCA_ON/scripts/5_refine_major/7celltype_class/2refine_scvi_general_micro_macro.sh  $line $s $res $rml $flavor $hvg $epoch $gl $span $kp
done

done 

