#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"

#query=${main}/HCA_ON/data/scvi/RNA/raw/ONONH_RNA_concated.h5ad

celltype="major" #$1
seed=7 #$2
res=1 #$3


#query="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_rpe_osf_tmcb_ononh_clean_final.h5ad"

query="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_all_rawcount.h5ad"

batch_key="sampleid"

label_key="none"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
mkdir -p $outdir
label="majorclass"
hvg=10000
#hvg=2000
sb=False
nlat=30


norm="t"
obs="none"
item="none"
max_epoch="none"
rmlist="rm_cluster_ret"
ntrg=10
mk="sanes_mk"
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}
detail="n"
fl="seurat_v3"
gl="zinb"
span=1
rank="none"
#rank="true"
#bname="${celltype}_hvg_${hvg}_fl_${fl}_res_${res}_sd_${seed}_ononh_retina_mg_sb"
#bname="${celltype}_hvg_${hvg}_fl_${fl}_res_${res}_sd_${seed}_clean_other"
bname="${celltype}_hvg_${hvg}_fl_${fl}_res_${res}_sd_${seed}_rmRet"
kp="none"
#rmcluster=""
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw_check.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $indir $rank $kp
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $indir $rank
