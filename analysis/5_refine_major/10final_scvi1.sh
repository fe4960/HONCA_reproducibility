#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"

batch_key="sampleid"

rmlist="none" #"rm_cluster_rs" #"none" #rm_cluster_res_1


celltype="major"
seed=7

outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/

query=${outdir}/ON_ONH_clean_new.h5ad


#query=${outdir}/ON_ONH_clean.h5ad
#query=${outdir}/major_hvg2000_epochnone_seurat_rs_1_clean_ON_ONH_scvi_trg.h5ad

#fl="seurat"
fl="seurat_v3"
hvg=2000
#hvg=10000
max_epoch="none"
#max_epoch=20
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
detail="t"
rank="none"
kp="none"
#indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_ON_ONH"
bname="${celltype}_hvg${hvg}_nonsb_epoch${max_epoch}_${fl}_rs_${res}_clean_ON_ONH_new"
#bname="${celltype}_hvg${hvg}_sb_epoch${max_epoch}_${fl}_rs_${res}_clean_ON_ONH_new"

python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp

