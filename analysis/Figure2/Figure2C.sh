#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug

orderfile=Figure2C_marker/majorclass_order.txt
markerfile=Figure2C_marker/majorclass_onecol.txt
outdir=scrnah5ad2seuratv5dotplot

# ATAC gene score
infile=ON_ATAC.h5ad
bname=$(basename "$infile" .h5ad)
color=blue
scrnah5ad2seuratv5dotplot -e r-seurat_v5.4.0 -n T -s T -d "$outdir" -b "${bname}" -m "$markerfile" -c "$orderfile" -g majorclass -H 3.5 -W 5.5 --clean -r "$color" -- "$infile"

# RNA
infile=ON_ATAC.h5ad
bname=$(basename "$infile" .h5ad)
color=red
scrnah5ad2seuratv5dotplot -e r-seurat_v5.4.0 -n T -d "$outdir" -b "${bname}" -m "$markerfile" -c "$orderfile" -g majorclass -H 3.5 -W 5.5 --clean -r "$color" -- "$infile"
