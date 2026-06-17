#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"

batch_key="sampleid"

rmlist="none" #"rm_cluster_rs" #"none" #rm_cluster_res_1


celltype=$1
seed=7

outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/

#query=${outdir}/${celltype}_subclass_new.h5ad

query=${outdir}/${celltype}_subclass.h5ad
#query=${outdir}/major_hvg2000_epochnone_seurat_rs_1_clean_ON_ONH_scvi_trg.h5ad

#fl="seurat"
fl="seurat_v3"
#hvg=2000
hvg=10000
max_epoch="none"
gl="zinb"
span=1
res=1


label_key="none"
label="majorclass1"
sb=True
nlat=30


norm="n"
obs="none"

item="none"
ntrg=10
mk="none"
detail="n"
rank="none"
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_span_${span}_subclass"

python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank

