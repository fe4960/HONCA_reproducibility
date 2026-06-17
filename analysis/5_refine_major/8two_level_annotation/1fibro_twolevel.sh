#!/bin/sh
source ~/.condainit
conda activate scvi
source $(jlbashwrapperpubutils -e)
main="/dfs3b/ruic20_lab/junw42"
celltype="Fibroblast"
input="${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/${celltype}_hvg10000_epochnone_seurat_v3_rs_0.4_clean_rm4_sb_scvi_trg.h5ad"
outdir="${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/twoLevel"
yaml="${main}/HCA_ON/data/5_refine_major/scvi/${celltype}/clean/twoLevel/two_level.yaml"
#yaml="${main}/HCA_ON/scripts/4_integrateRNA/9two_level2.yaml"
mkdir $outdir

scrnah5adscvi2leidenbyreswkfl  -d "$outdir"  -c $yaml -t 5 -- $input
