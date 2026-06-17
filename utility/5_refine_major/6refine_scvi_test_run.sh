#!/bin/sh

main="/dfs3b/ruic20_lab/junw42"
err=${main}/HCA_ON/scripts/5_refine_major/6refine_scvi/
mkdir -p $err
#s=12345
s=7
res=1
for c in Astrocyte
do
rml="rm_cluster_res_1_new"
#flavor="seurat"
flavor="seurat_v3"
num=2000
epoch="none"
gl="nb"
span=0.3
l=${c}_hvg${num}_epoch20_sd${s}_res${res}_${flavor}_${gl}
sbatch -p gpu --account=ruic20_lab_gpu --time=0-12 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/6refine_scvi_general.sh $c $s $res $rml $flavor $num $epoch $gl $span

#flavor="seurat"
#num=2000
#epoch=20
#l=${c}_hvg${num}_epoch20_sd${s}_res${res}_${flavor}
#sbatch -p gpu --account=ruic20_lab_gpu --time=0-12 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/6refine_scvi_general.sh $c $s $res $rml $flavor $num $epoch


#flavor="seurat"
#num=10000
#epoch="none"
#l=${c}_hvg${num}_epoch20_sd${s}_res${res}_${flavor}
#sbatch -p gpu --account=ruic20_lab_gpu --time=0-12 --mem=100GB --error=${err}/${l}.err --output=${err}/${l}.out --gres=gpu:1 ${main}/HCA_ON/scripts/5_refine_major/6refine_scvi_general.sh $c $s $res $rml $flavor $num $epoch




done 

