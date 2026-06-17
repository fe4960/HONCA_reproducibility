#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"
outdir=${main}/HCA_ON/data/scvi/RNA/scvi/
mkdir -p $outdir
yaml=${main}/HCA_ON/scripts/4_integrateRNA/3scvi_config.yaml
indir=${main}/HCA_ON/data/scvi/RNA/raw/
bname="ONONH_HGV2k_scvi"
scrnah5adfiles2scviwkfl -d "$outdir" -b "$bname" -e scvi -t 1 -c $yaml -- "$indir"/*.h5ad
