#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"

#query=${main}/HCA_ON/data/scvi/RNA/raw/ONONH_RNA_concated.h5ad
query=${main}/HCA_ON/data/4_scvi/RNA/scvi/cellclass/Muller_glia.h5ad

celltype="major"
batch_key="sampleid"
bname="MG_hvg2k_epoch20"
label_key="none"
outdir=${main}/HCA_ON/data/4_scvi/RNA/scvi/cellclass/MG
mkdir $outdir
label="none"
#hvg=10000
hvg=2000
sb=False
nlat=30
norm="n"


python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm
