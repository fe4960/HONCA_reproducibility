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

query=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_new2_rm_clu1_8.h5ad
bname="${celltype}_subclass_new2_rm_clu1_8_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_seed_${seed}_clean_sb"


label_key="none"
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/
mkdir -p $outdir
label="majorclass1"
sb=False
nlat=30


norm="n"
#obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_v3_rs_0.3_clean_rm4_6_13_sb_res_0.3.obs.gz
obs="none"
#obs=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_subclass_new2_on_hvg${hvg}_epoch${max_epoch}_${fl}_rs_${res}_clean_sb.obs.gz
item="leiden" #"none"
ntrg=10
mk="sanes_mk"
detail="t"
indir=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/
dirrm=${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/

rank="true"
#covar="total_counts"
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank
python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp
#echo " python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw_nCount.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp $covar  "
#python /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/4scVI_cluster_raw_nCount.py $query $celltype ${batch_key} $bname ${label_key} $outdir $label $hvg $sb $nlat $norm $obs $item $max_epoch $rmlist $seed $ntrg $mk $indir $res $detail $fl $gl $span $dirrm $rank $kp $covar
