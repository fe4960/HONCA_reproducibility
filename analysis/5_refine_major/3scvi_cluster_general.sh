#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"

#query=${main}/HCA_ON/data/scvi/RNA/raw/ONONH_RNA_concated.h5ad
#query=${main}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_rawCount_rename_majorclass.h5ad
query=${main}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_rawCount_rename_majorclass_final.h5ad

celltype=$1
batch_key="sampleid"
#bname="ONONH_hvg2k_epochDefault_new"
#bname="ONONH_hvg2k_epoch400_new"

label_key="none"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
mkdir -p $outdir
label="majorclass1"

#label="majorclass"
#hvg=10000
hvg=2000
sb=False
nlat=30


norm="t"
obs="none"
item="none"
#max_epoch=20
max_epoch="none"
rmlist="none"

ntrg=10
mk="sanes_mk"  #"none" #sanes_mk"
detail="t"
seed=7
res=1
detail="t"
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}
fl="seurat_v3"
gl="zinb"
bname="${celltype}_hvg2k_epoch${max_epoch}_${fl}_${gl}"
span=1
dirrm=$indir
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm


#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster1.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span

#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch
