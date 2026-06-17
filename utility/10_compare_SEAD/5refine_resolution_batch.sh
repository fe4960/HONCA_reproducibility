#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"

batch_key="sampleid"

#err="${main}/HCA_ON/scripts/5_refine_major/5refine_resolution_sb"

err="${main}/HCA_ON/scripts/10_compare_SEAD/5refine_resolution_batch"

mkdir -p $err

#for cell in Astrocyte_GSE267301_merge Oligo_GSE267301_merge
#for cell in Oligo_GSE267301_merge
for cell in Astrocyte Microglia Oligodendrocyte Oligodendrocyte_precursor_cell 
do	
#outdir="HCA_ON/data/scvi/RNA/raw/${cell}" 
outdir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/${cell}_PFC_ONONH/
#outdir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean
#indir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/

indir=$outdir

for res in 0.3  0.4 0.6 
do

for kp in none #kp_cl_3 

do
#h5ad=${outdir}/${cell}_hvg10000_epochnone_seurat_v3_rs_0.4_rm3_4_5_none_scvi_trg.h5ad
#h5ad=${outdir}/${cell}_hvg_10000_epoch_none_seurat_v3_rs_1_scvi_trg.h5ad
h5ad=${outdir}/${cell}_PFC_ONONH_hvg_10000_epoch_none_seurat_v3_rs_1_sb_scvi_trg.h5ad
#bname=${cell}_hvg_10000_epoch_none_seurat_v3_rs_${res}
#bname=${cell}_hvg10000_epochnone_seurat_v3_rs_${res}_rm3_4_5_none
bname=${cell}_PFC_ONONH_hvg_10000_epoch_none_seurat_v3_rs_${res}_sb
rmlist="none" #"rm_cluster_rs" #"none" #rm_cluster_res_1
mk="none" #save="n"
save="n"
sbatch --mem=30GB -p free --account=ruic20_lab --time=0-2 --output=${err}/${cell}_${res}_${kp}.out --error=${err}/${cell}_${res}_${kp}.err  HCA_ON/scripts/5_refine_major/5refine_resolution.sh $h5ad $res ${batch_key} $bname $outdir $rmlist $mk $save $indir
done
done
done
