#!/bin/sh
source ~/.condainit
conda activate scvi
source $(jlbashwrapperpubutils -e)
main="/dfs3b/ruic20_lab/junw42"
input="${main}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_hvg2k_epoch20_new_scvi.h5ad"
outdir="${main}/HCA_ON/data/4_scvi/RNA/scvi/twoLevel"
yaml="${main}/HCA_ON/scripts/4_integrateRNA/9two_level1.yaml"
#yaml="${main}/HCA_ON/scripts/4_integrateRNA/9two_level2.yaml"
mkdir $outdir

scrnah5adscvi2leidenbyreswkfl  -d "$outdir"  -c $yaml -t 2 -- $input
