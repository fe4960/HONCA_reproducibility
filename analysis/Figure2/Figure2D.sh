#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug

infile=Figure2D_ON_ATAC_peakset.rds
bname=$(basename "$infile" _peakset.rds)
outdir=archrpeaksetrds2peaktypebarplot
archrpeaksetrds2peaktypebarplot -e r-archr_v1.0.3 -d "$outdir" -b "$bname" -W 3.5 -H 3.5 -y "'Number of OCRs'" -- "$infile"
