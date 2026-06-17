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
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_LC6.h5ad
#query=${main}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/Oligo_GSE267301_merge_rawcount_rm_res0.4_clu3.h5ad
query="/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/${celltype}_PFC_ONONH.h5ad"
#query=${main}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge.h5ad
batch_key="sampleid"
#bname="Oligo_GSE267301_merge_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_rm3_${kp}"
bname="${celltype}_PFC_ONONH_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_rm_${rmlist}_kp_${kp}"
label_key="none"
outdir=${main}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/
mkdir -p $outdir
label="majorclass"
sb=False
nlat=30


norm="t"
#obs=${outdir}/Oligo_GSE267301_merge_hvg_10000_epoch_none_seurat_v3_rs_0.4_res_0.4.obs.gz
obs="none" #${outdir}/Oligo_GSE267301_merge_hvg10000_epochnone_seurat_v3_rs_0.4_rm3_none_scvi_trg.obs.gz
obs_kp="none"
item="leiden" #"none"
ntrg=10
mk="none"
detail="n"
indir=$outdir #${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
dirrm=$outdir #{main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
#dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/

rank="true"
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp ${obs_kp}
