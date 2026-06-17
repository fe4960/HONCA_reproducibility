#!/bin/sh

main="/dfs3b/ruic20_lab/junw42"
err=${main}/HCA_ON/scripts/10_compare_SEAD/2refine_scvi_refine_run
mkdir -p $err
s=7
res=1
#for line in Astro Micro Oligo  
#for line in GSE267301 #Astrocyte   #Oligodendrocyte  
#for line in Astrocyte_GSE267301_merge Oligo_GSE267301_merge 
for line in GSE167494_more

#for line in GSE167494 
#for line in Astrocyte_PFC_ONONH Microglia_PFC_ONONH Oligodendrocyte_PFC_ONONH Oligodendrocyte_precursor_cell_PFC_ONONH    
do
rml="none"
for kp in none   #kp_cl_2 kp_cl_0  kp_cl_16 kp_cl_3 kp_cl_4 
do
flavor="seurat"
#flavor="seurat_v3"
hvg=10000
#hvg=10000
epoch="none" 
gl="zinb" 
span=1
l=${line}_hvg${num}_epoch${epoch}_sd${s}_res${res}_${flavor}

###sbatch -p gpu --account=ruic20_lab_gpu --time=0-2 --mem=100GB --gres=gpu:1 --error=${err}/${l}.err --output=${err}/${l}.out  ${main}/HCA_ON/scripts/10_compare_SEAD/2refine_scvi_general.sh  $line $s $res $rml $flavor $hvg $epoch $gl $span $kp

sbatch -p free-gpu --account=ruic20_lab_gpu --time=0-2 --mem=50GB --gres=gpu:1 --error=${err}/${l}.err --output=${err}/${l}.out  ${main}/HCA_ON/scripts/10_compare_SEAD/2refine_scvi_general_pub.sh  $line $s $res $rml $flavor $hvg $epoch $gl $span $kp


done

done 

