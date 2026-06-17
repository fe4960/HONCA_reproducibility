#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"

batch_key="sampleid"

#err="${main}/HCA_ON/scripts/5_refine_major/5refine_resolution"
err="${main}/HCA_ON/scripts/5_refine_major/5refine_resolution_sb"

mkdir -p $err

for cell in Mural_cell
do	
outdir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean
indir=${main}/HCA_ON/data/5_refine_major/scvi/${cell}/

for res in 0.6 0.7 0.8 0.9 1 #0.2 0.3 0.4 0.5 
do

for kp in none #kp_cl_3 #kp_cl_1   #kp_cl_1 kp_cl_3 kp_cl_4 kp_cl_06

do

#h5ad=${outdir}/${cell}_hvg10000_epochnone_seurat_rs_0.6_rmnone_sb_scvi_trg.h5ad
h5ad=${outdir}/${cell}_hvg10000_epochnone_seurat_rs_0.6_rmrm_clu1_sb_scvi_trg.h5ad
bname=${cell}_hvg10000_epochnone_seurat_rs_${res}_rmrm_clu1_sb

#bname=${cell}_hvg10000_epochnone_seurat_rs_${res}_rmnone_sb
rmlist="none" #"rm_cluster_rs" #"none" #rm_cluster_res_1
mk="sanes_mk" #save="n"
save="n"
sbatch --mem=30GB -p free-gpu --account=ruic20_lab --time=0-3 --output=${err}/${cell}_${res}_${kp}.out --error=${err}/${cell}_${res}_${kp}.err  HCA_ON/scripts/5_refine_major/5refine_resolution.sh $h5ad $res ${batch_key} $bname $outdir $rmlist $mk $save $indir
done
done
done
