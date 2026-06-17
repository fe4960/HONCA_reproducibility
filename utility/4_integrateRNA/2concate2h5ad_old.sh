#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"
outdir=${main}/HCA_ON/data/scvi/RNA/raw/
bname="ONONH_RNA_concated"
mkdir -p $outdir
fname=${main}/HCA_ON/data/sample_list/RNA_sample_list_meta.gz
scrnascanpyconcat2h5ad.sh -t 12 -d "$outdir" -b "$bname" -- $fname
