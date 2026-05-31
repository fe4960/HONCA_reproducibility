#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug

# 1. DAR
archrprjdir=ON_ATAC
bname=$(basename "$archrprjdir")
outdir=archrproject2getmarkerfeatures
archrproject2getmarkerfeatures -e r-archr_v1.0.3 -d "$outdir" -b "${bname}" -g majorclass -n PeakMatrix -m wilcoxon -z T -- "$archrprjdir"

# 2. Heatmap
infile=archrproject2getmarkerfeatures/ON_ATAC.rds
bname=$(basename "$infile" .rds)
outdir=archrmarkerfeature2heatmap
archrmarkerfeature2heatmap -e r-archr_v1.0.3 -d "$outdir" -b "$bname" -c '"FDR <= 0.01 & Log2FC >= 0.5"' -W 6 -H 6 -g F -- "$infile"
