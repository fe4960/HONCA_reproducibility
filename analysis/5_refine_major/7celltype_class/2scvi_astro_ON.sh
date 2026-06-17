#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"

#query=${main}/HCA_ON/data/scvi/RNA/raw/ONONH_RNA_concated.h5ad

celltype=$1
seed=$2
res=$3


query="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte_ON/Astro_ON.h5ad"

batch_key="sampleid"
#bname="${celltype}_hvg10k_epoch20_clean"

#bname="${celltype}_hvg2k_epoch20_clean_sd${seed}_clean_clu4"
#bname="${celltype}_seed_${seed}_res_${res}_clean_other"
#bname="${celltype}_seed_${seed}_res_${res}_clean_other_clean"
#bname="${celltype}_seed_${seed}_res_${res}_clean_ononh_retina"

#bname="${celltype}_hvg2k_epoch20_clean_sd${seed}"
#bname="ONONH_hvg2k_epochDefault_new"
#bname="ONONH_hvg2k_epoch400_new"

label_key="none"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
mkdir -p $outdir
label="majorclass"
#hvg=10000
hvg=2000
sb=False
nlat=30


norm="n"
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
rank="none"
#rank="true"
bname="${celltype}_hvg_${hvg}_fl_${fl}_res_${res}_sd_${seed}_ononh_retina_mg_sb"

#rmcluster=""
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $indir $rank
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $indir $rank
