#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug

# 1. Integrate ATAC and RNA
atac=ON_ATAC
rna=ON_RNA.h5ad
bname=$(basename "$atac")
outdir=archratacrnasupervised2addgeneintegrationmatrix
archratacrnasupervised2addgeneintegrationmatrix -e r-archr_v1.0.3 -d "$outdir" -b "$bname" -r "$rna" -l majorclass -L majorclass -c Harmony -n 10000 -a cca -t 20 -- "$atac"

# 2. Add peak2gene links
archrdir=archratacrnasupervised2addgeneintegrationmatrix/ON_ATAC
bname=$(basename "$archrdir")
outdir=archrproject2addpeak2genelinks
archrproject2addpeak2genelinks -e r-archr_v1.0.3 -d "$outdir" -b "$bname" -m GeneIntegrationMatrix -c Harmony -- "$archrdir"

# 3. plot heatmap
infile=archrproject2addpeak2genelinks/ON_ATAC
bname=$(basename "$infile")
outdir=archrproject2plotpeak2geneheatmap
archrproject2plotpeak2geneheatmap -e r-archr_v1.0.3 -d "$outdir" -b "$bname" -g majorclass -c 0.5 -f 0.01 -k 15 -W 10 -H 10 -G F -s T -- "$f"
