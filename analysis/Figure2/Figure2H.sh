#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug

# 1. Use the SCENIC+ results in Figure2G.sh
atac=ON_ATAC.h5ad
rna=ON_RNA.h5ad
bname=$(basename "$atac" .h5ad)
outdir=atacrnah5ad2scenicplusregulonwkfl
configfile=Figure2G_config.yaml
atacrnah5ad2scenicplusregulonwkfl -c "$configfile" -e scenicplus_v1.0a2 -d "$outdir" -b "$bname" -t 40 -r "$rna" -- "$atac"

# 2. Plot RSS
infile=atacrnah5ad2scenicplusregulonwkfl/scenicplush5mu2eregulonrss/ON_ATAC.txt.gz
bname=$(basename "$infile" .txt.gz)
outdir=scenicplusrsstxt2plotrss
scenicplusrsstxt2plotrss -e scenicplus_v1.0a2 -d "$outdir" -b "$bname" -n 10 -m 4 -W 3.5 -H 4 -- "$infile"
