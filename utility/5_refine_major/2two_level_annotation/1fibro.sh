#!/bin/sh
source ~/.condainit
conda activate scvi
source $(jlbashwrapperpubutils -e)
main="/dfs3b/ruic20_lab/junw42"
#input=${main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new3.h5ad
#input=${main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.4_clean_sb_rm0_12_15_16_scvi_trg.h5ad
#input=${main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7_scvi_trg.h5ad
input=${main}/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_hvg10000_epochnone_seurat_v3_rs_1_clean_sb_none_seed_7_rmRPE_scvi_trg.h5ad

#input=${main}/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_hvg10000_epochnone_seurat_rs_1_clean_sb_none_seed_7_rmRPE_scvi_trg.h5ad
#input=${main}/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/
#outdir="${main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/twoLevel_pan_eye_Astro_seed12345_new"
#outdir="${main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/twoLevel_pan_eye_Astro_seed7_new"
#outdir="${main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/twoLevel_pan_eye_Astro_seed12345_new_update"
#outdir="${main}/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/twoLevel_Fibro_seed7"
outdir="${main}/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/twoLevel_Fibro_seed7_seurat_v3"

#outdir="${main}/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/twoLevel_Fibro_seed12345"
#outdir="${main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/twoLevel_pan_eye_Astro"
mkdir -p ${outdir}
#outdir="${main}/HCA_ON/data/4_scvi/RNA/scvi/cellclass/twoLevel"
#yaml="${main}/HCA_ON/scripts/4_integrateRNA/9two_level.yaml"
#yaml="${main}/HCA_ON/scripts/5_refine_major/2two_level.yaml"
#yaml="${main}/HCA_ON/scripts/5_refine_major/2two_level.yaml"
#yaml="${main}/HCA_ON/scripts/5_refine_major/2two_level_annotation/1astro.yaml"
#yaml="${main}/HCA_ON/scripts/5_refine_major/2two_level_annotation/1fibro.yaml"
yaml="${main}/HCA_ON/scripts/5_refine_major/2two_level_annotation/1fibro.yaml"

mkdir $outdir

scrnah5adscvi2leidenbyreswkfl  -d "$outdir"  -c $yaml -- $input
