#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"

batch_key="sampleid"

#err="${main}/HCA_ON/scripts/5_refine_major/5refine_resolution"
err="${main}/HCA_ON/scripts/5_refine_major/5refine_resolution_sb"

mkdir -p $err

for cell in Endothelial_cell 
do	
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean
indir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/

for res in 0.8 1.2 1.3 1.4 1.5  
do

for kp in none #kp_cl_3 #kp_cl_1   #kp_cl_1 kp_cl_3 kp_cl_4 kp_cl_06

do

h5ad=${outdir}/${cell}_hvg10000_epochnone_seurat_v3_rs_1_subclass_sb_seurat_rmRPE_scvi_trg.h5ad
#bname=${cell}_hvg10000_epochnone_seurat_rs_${res}_clean_rm_sb_rm_4_6_none_scvi_trg.h5ad
bname=${cell}_hvg10000_epochnone_seurat_v3_rs_${res}_subclass_sb_seurat_rmRPE
rmlist="none" #"rm_cluster_rs" #"none" #rm_cluster_res_1
mk="sanes_mk" #save="n"
save="n"
sbatch --mem=10GB -p standard --account=ruic20_lab --time=0-3 --output=${err}/${cell}_${res}_${kp}.out --error=${err}/${cell}_${res}_${kp}.err  HCA_ON/scripts/5_refine_major/5refine_resolution.sh $h5ad $res ${batch_key} $bname $outdir $rmlist $mk $save $indir
done
done
done
