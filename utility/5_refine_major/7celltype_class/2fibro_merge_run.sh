#!/bin/sh

main="/dfs3b/ruic20_lab/junw42"
err=${main}/HCA_ON/scripts/5_refine_major/6refine_scvi/
mkdir -p $err
s=7
res=1
for c in Fibroblast
do
l=${c}_hvg2k_epoch20_sd${s}_res${res}_other_final
#sbatch -p gpu --account=ruic20_lab_gpu --time=0-12 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/7celltype_class/2fibroblast_otherPart.sh $c $s $res

#sbatch -p gpu --account=ruic20_lab_gpu --time=0-2 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/7celltype_class/3astro_ON_retina.sh $c $s $res
sbatch -p gpu --account=ruic20_lab_gpu --time=0-2 --mem=150GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/7celltype_class/2fibro_merge.sh $c $s $res


done 

#< $file
