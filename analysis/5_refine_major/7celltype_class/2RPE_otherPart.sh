#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"

#query=${main}/HCA_ON/data/scvi/RNA/raw/ONONH_RNA_concated.h5ad

celltype="major"
seed=7
res=0.3


query="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/RPE/clean/RPE_rpe_osf_tmcb_ononh.h5ad"
batch_key="sampleid"

bname="${celltype}_seed_${seed}_res_${res}_clean_other"


label_key="none"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
mkdir -p $outdir
label="majorclass"
hvg=2000
sb=False
nlat=30


norm="t"
obs="none"
item="none"
max_epoch="none"
#rmlist="rm_cluster_res_1"
#rmlist="rm_cluster_res_rm_res0.4_clu4"
#rmlist="rm_cluster_res_rm_res0.4_clu4_clu5"
rmlist="none"
ntrg=10
mk="sanes_mk"
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}
detail="n"
fl="seurat_v3"
gl="zinb"
span=1
dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
#rmcluster=""
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm
