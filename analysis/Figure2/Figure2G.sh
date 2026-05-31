#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug

# 1. Run SCENIC+ pipeline
atac=ON_ATAC.h5ad
rna=ON_RNA.h5ad
bname=$(basename "$atac" .h5ad)
outdir=atacrnah5ad2scenicplusregulonwkfl
configfile=Figure2G_config.yaml
atacrnah5ad2scenicplusregulonwkfl -c "$configfile" -e scenicplus_v1.0a2 -d "$outdir" -b "$bname" -t 40 -r "$rna" -- "$atac"

# 2. generate heatmap
infile=atacrnah5ad2scenicplusregulonwkfl/ataccistopic2runscenicpluspipeline/ON_ATAC/scplusmdata.h5mu
bname=ON_ATAC
featurefile=atacrnah5ad2scenicplusregulonwkfl/scenicplusreguloncorrelation2subsetfeature/ON_ATAC.txt.gz
outdir=scenicplush5mu2eregulonheatmap
scenicplush5mu2eregulonheatmap -e scenicplus_v1.0a2 -d "$outdir" -b "$bname" -l majorclass -W 20 -H 5 -n "$featurefile" -- "$infile"
