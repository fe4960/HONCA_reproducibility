#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"


celltype=Oligodendrocyte_precursor_cell
seed=7
res=0.4
rmlist="rm_cluster_rs1_1"

#rmlist="rm_cluster_rs_0.4"
#fl=seurat_v3
fl=seurat
hvg=10000
max_epoch="none"
gl="zinb"
span=1
kp=none
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_v3_rs_0.4_clean_rm4_15_scvi_trg.h5ad
query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_seurat.h5ad
batch_key="sampleid"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm1_sb_${kp}"
#bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_res0.3rm1_sb_${kp}"
bname="${celltype}_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_rm1_sb_${kp}"



label_key="none"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
mkdir -p $outdir
label="majorclass1"
#label="majorclass1"
sb=True
nlat=30


norm="n"
##obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_res_0.4_clean_sb_res_0.4.obs.gz
#obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_res_0.3_clean_sb_res_0.3.obs.gz
#obs="none"
obs="${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_rs_1_clean_sb_none.obs.gz"
obs_kp="none"
item="leiden" #"none"
ntrg=10
mk="sanes_mk"
detail="t"
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/

rank="true"
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp $obs_kp
