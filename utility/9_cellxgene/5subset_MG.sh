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
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg2000_epochnone_seurat_v3_rs_0.4_clean_rm_scvi_trg.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg20000_epochnone_seurat_v3_rs_0.4_clean_rm4_15_16_scvi_trg.h5ad

#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_res_1_clean_rm_res_1.0.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_res_1_clean_res_1.0.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg_2000_fl_seurat_v3_res_0.4_sd_7_ononh_retina_nonsb_scvi_trg.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/${celltype}_hvg2k_epochnone_seurat_v3_zinb_scvi_trg.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_res_1_clean_res_1.0.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_v3_rs_1_span_1_subclass_scvi_trg.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_LC6.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg_10000_fl_seurat_v3_res_1_sd_7_rmRet_scvi_trg.h5ad
query=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_all_rawcount.h5ad

#query=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_all_normcount.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_v3_rs_0.6_clean_rm4_6_13_sb_scvi_trg.h5ad

batch_key="sampleid"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_sb"

#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm4_15"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm4_15_sb"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm4_6_13_sb"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm4_sb_${kp}"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm4_6_13_sb_${kp}"
bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_sb_${kp}_MG"

#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm4_6_sb"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm4_6_nonsb"

#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm4_15_16"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm4_15_16_rm3"

#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm_rm"


label_key="none"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
mkdir -p $outdir
label="majorclass"
#label="majorclass1"
sb=False
nlat=30

norm="t"
#norm="n"
#obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_v3_rs_1_span_1_subclass.obs.gz
#obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_v3_res_0.3_clean_rm4_sb_res_0.3.obs.gz
#obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_v3_rs_0.3_clean_rm4_6_13_sb_res_0.3.obs.gz
#obs="none"
obs="HCA_ON/data/5_refine_major/scvi/major/clean/major_hvg_10000_fl_seurat_v3_res_1_sd_7_rmRet.obs.gz"
#obs="${outdir}/${celltype}_hvg_2000_fl_seurat_v3_res_0.3_sd_7_ononh_retina_nonsb_res_0.3.obs.gz"

#obs="${outdir}/${celltype}_hvg_2000_fl_seurat_v3_res_0.3_sd_7_ononh_retina_nonsb.obs.gz"
#obs="${outdir}/${celltype}_res_1_clean_res_1.0.obs.gz" #"none"
#obs="${outdir}/${celltype}_res_1_clean_rm_rm_res_1.0.obs.gz" #"none"
#obs="${outdir}/${celltype}_res_0.5_clean_rm_res_0.5.obs.gz" #"none"

item="leiden" #"none"
ntrg=10
#mk="sanes_mk"
mk="none"
detail="t"
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
#dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/

rank="true"
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp
