#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"


celltype=$1
seed=$2
res=$3
rmlist=$4
fl=$5
hvg=$6
max_epoch=$7
gl=$8
span=$9
kp=${10}
query=${main}/HCA_ON/data/9_Brain_h5ad/merged/${celltype}.h5ad
batch_key="sampleid"

bname="${celltype}_hvg_${hvg}_epoch_${max_epoch}_${fl}_rs_${res}"


label_key="none"
#outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
outdir=${main}/HCA_ON/data/9_Brain_h5ad/merged/
mkdir -p $outdir
label="majorclass"
sb=False
nlat=30


norm="n"
#obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_v3_rs_0.3_clean_rm4_6_13_sb_res_0.3.obs.gz
obs="none"
item="leiden" #"none"
ntrg=10
#mk="sanes_mk"
mk="none"
detail="n"
#indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
indir=$outdir
#dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
dirrm=$indir
rank="true"
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/10_compare_SEAD/2refine_scvi_general.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp
