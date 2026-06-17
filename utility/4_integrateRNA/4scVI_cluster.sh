#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"

#query=${main}/HCA_ON/data/scvi/RNA/raw/ONONH_RNA_concated.h5ad
query=${main}/HCA_ON/data/4_scvi/RNA/raw/ONONH_RNA_full.h5ad

celltype="major"
batch_key="sampleid"
#bname="ONONH_hvg2k_epoch20_new"
#bname="ONONH_hvg2k_epochDefault_new"
#bname="ONONH_hvg2k_epoch400_new"
bname="ONONH_hvg2k_epoch20_new_correct"

label_key="none"
outdir=${main}/HCA_ON/data/4_scvi/RNA/scvi/
label="none"
#hvg=10000
hvg=2000
sb=False
nlat=30


norm="t"
obs="none"
item="none"

#max_epoch=400
max_epoch=20
rmlist="none"
seed=7
ntrg=10
mk="none"
indir=$outdir
res=1
detail="t"
fl="seurat"
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch

python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_hvg.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl
