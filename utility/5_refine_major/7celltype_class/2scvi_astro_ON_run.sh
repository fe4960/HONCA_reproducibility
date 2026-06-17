#!/bin/sh

main="/dfs3b/ruic20_lab/junw42"
err=${main}/HCA_ON/scripts/5_refine_major/6refine_scvi/
mkdir -p $err
#file=${main}/HCA_ON/scripts/5_refine_major/3scvi_cluster_general_cellclass
#while read c
#s=12345
s=7
#res=0.4
res=0.4
#res=0.3
#for c in 'Endothelial cell' 'Mural cell' 'NK/T cell' 'Oligodendrocyte precursor cell' 'Schwann cell'
#for c in Fibroblast #Astrocyte
#for c in Fibroblast
for c in Astrocyte
do
l=${c}_hvg2k_epoch20_sd${s}_res${res}_other
#sbatch -p gpu --account=ruic20_lab_gpu --time=0-12 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/7celltype_class/2fibroblast_otherPart.sh $c $s $res

sbatch -p gpu --account=ruic20_lab_gpu --time=0-2 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/7celltype_class/2scvi_astro_ON.sh $c $s $res


done 

#< $file
