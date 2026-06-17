#!/bin/sh

main="/dfs3b/ruic20_lab/junw42"
err=${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general/
mkdir -p $err
file=${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general_cellclass
#while read line
for line in B_cell #Rod Cone HC #Microglia   #NK_T_cell
do
l=${line}_hvg2k_epochNone
sbatch -p gpu --account=ruic20_lab_gpu --time=0-6 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general.sh $line
done 

#< $file

#c="Astrocyte"
#l=Astrocyte_hvg2K_epochNone
#sbatch -p gpu --account=ruic20_lab_gpu --time=0-10 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general.sh $c


#c="Endothelial_cell"
#l=Endothelial_hvg20_epoch20
#sbatch -p gpu --account=ruic20_lab_gpu --time=0-10 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general.sh $c


#c="Mural_cell"
#l=Mural_hvg20_epoch20
#sbatch -p gpu --account=ruic20_lab_gpu --time=0-10 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general.sh $c


#c="NK/T_cell"
#l=NKT_hvg20_epoch20
#sbatch -p gpu --account=ruic20_lab_gpu --time=0-10 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general.sh $c

#c="Oligodendrocyte_precursor_cell"
#l=OPC_hvg20_epoch20
#sbatch -p gpu --account=ruic20_lab_gpu --time=0-10 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general.sh $c


#c="Schwann_cell"
#l=Schwann_hvg20_epoch20
#sbatch -p gpu --account=ruic20_lab_gpu --time=0-10 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general.sh $c
