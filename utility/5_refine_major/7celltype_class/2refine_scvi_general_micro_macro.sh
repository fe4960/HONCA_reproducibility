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
query=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/Micro_Macro_clean_final.h5ad

batch_key="sampleid"
bname="Micro_Macro_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_sb_${kp}"


label_key="none"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
mkdir -p $outdir
label="majorclass1"
#label="majorclass1"
sb=False
nlat=30


norm="n"
obs="none"
item="leiden" #"none"
ntrg=10
mk="none" #"sanes_mk"
detail="t"
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/

rank="true"
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp
