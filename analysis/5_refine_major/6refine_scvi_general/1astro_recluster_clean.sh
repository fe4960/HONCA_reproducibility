#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"


celltype=Astrocyte
seed=7
#seed=12345
res=0.6
#rmlist=rm_cluster_0_15_16
rmlist=rm_cluster_5_9

#rmlist=rm_cluster_0_12_15_16
fl="seurat"
hvg=10000
max_epoch="none"
gl="zinb"
span=1
kp=none
batch_key="sampleid"

###query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_new2_ret_onh.h5ad
######bname="${celltype}_subclass_new2_ret_onh_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_sb"

######query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_new2_on.h5ad
#######bname="${celltype}_subclass_new2_on_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_nCount_covar_clean_sb"
######bname="${celltype}_subclass_new2_on_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_rm1_8_clean_sb"


#####query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_new_rm_non_spec.h5ad
######bname="${celltype}_subclass_new_rm_non_spec_onh_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_sb"

########query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_new.h5ad
#######bname="${celltype}_subclass_new_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_covar_nCount_clean_sb"

#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_new2_rm_clu1_8.h5ad
#query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_rs_0.4_clean_sb_none_scvi_trg.h5ad
query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad

#bname="${celltype}_subclass_new2_rm_clu1_8_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_seed_${seed}_clean_sb"
#bname="${celltype}_hvg10000_epochnone_seurat_rs_0.4_clean_sb_rm0_12_15_16" #######
#bname=${celltype}_hvg${hvg}_epoch${max_epoch}_seurat_rs_${res}_clean_sb_rm0_12_15_16
#bname=${celltype}_hvg${hvg}_epoch${max_epoch}_seurat_rs_${res}_clean_sb_rm0_15_16
bname=${celltype}_hvg${hvg}_epoch${max_epoch}_seurat_rs_${res}_clean_sb_rm0_15_16_seed12345_rs0.8_rm5_9_seed_${seed}

label_key="none"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
mkdir -p $outdir
label="majorclass1"
sb=False
nlat=30


norm="n"
#obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_v3_rs_0.3_clean_rm4_6_13_sb_res_0.3.obs.gz
#obs="none"
#obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_rs_1_clean_sb_none.obs.gz
obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_rs_0.8_clean_sb_rm0_15_16.obs.gz
#obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_new2_on_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_sb.obs.gz
item="leiden" #"none"
ntrg=10
mk="sanes_mk"
detail="t"
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/

rank="true"
kp_obs="none"
#covar="total_counts"
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp $kp_obs
#echo " python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw_nCount.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp $covar  "
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw_nCount.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp $covar
