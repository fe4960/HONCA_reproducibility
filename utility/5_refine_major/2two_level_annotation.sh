#!/bin/sh
source ~/.condainit
conda activate scvi
source $(jlbashwrapperpubutils -e)
main="/dfs3b/ruic20_lab/junw42"
#input="${main}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_hvg2k_epoch20_new_scvi.h5ad"
#input="${main}/HCA_ON/data/4_scvi/RNA/scvi/cellclass/Muller_glia.h5ad"
#input="${main}/HCA_ON/data/4_scvi/RNA/scvi/cellclass/MG/MG_hvg2k_epoch20_scvi.h5ad"
#input=${main}/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_hvg_10000_fl_seurat_v3_res_0.4_sd_7_clean_final_sb_scvi_trg.h5ad
input=${main}/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/
outdir="${main}/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/twoLevel_pan_eye_Fibro"
mkdir -p ${outdir}
#outdir="${main}/HCA_ON/data/4_scvi/RNA/scvi/cellclass/twoLevel"
#yaml="${main}/HCA_ON/scripts/4_integrateRNA/9two_level.yaml"
#yaml="${main}/HCA_ON/scripts/5_refine_major/2two_level.yaml"
yaml="${main}/HCA_ON/scripts/5_refine_major/2two_level.yaml"

mkdir $outdir

scrnah5adscvi2leidenbyreswkfl  -d "$outdir"  -c $yaml -- $input
