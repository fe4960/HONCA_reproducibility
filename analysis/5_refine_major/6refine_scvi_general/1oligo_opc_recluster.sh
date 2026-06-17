#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"


celltype="major"
seed=7
res=0.4
rmlist="none"

fl=seurat
hvg=10000
max_epoch="none"
gl="zinb"
span=1
kp=none
query=${main}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/oligo_opc_seurat_new.h5ad

batch_key="sampleid"
bname="oligo_opc_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_sb_${kp}"



label_key="none"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/
mkdir -p $outdir
label="majorclass1"
#label="majorclass1"
sb=True
nlat=30


norm="n"
obs="none"
item="leiden" #"none"
ntrg=10
mk="sanes_mk"
detail="t"
indir=${main}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/
dirrm=${main}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/

rank="true"
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp
