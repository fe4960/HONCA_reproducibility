#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

source trapdebug

# 1. scGLUE preprocessing
atac=ON_ATAC.h5ad
rna=ON_RNA.h5ad
gtf=cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz
outdir=gluernaatac2preproc
nhvg=10000
bname=ON_coembed
gluernaatac2preproc -e scglue -d "$outdir" -b "$bname" -r "$rna" -a "$atac" -g "$gtf" --rnabatchkey sampleid --rnahvgntop "$nhvg" --rnahvgflavor seurat -c 0.2 --gtfby gene_name -R -A -G


# 2. scGLUE co-embedding
indir=gluernaatac2preproc
outdir=glueh5adgraphatacrnatype2coembed
bname=ON_coembed
rnarep=X_scVI
atacrep=Harmony
rna=$indir/${bname}_RNA.h5ad
atac=$indir/${bname}_ATAC.h5ad
graphguidance=$indir/${bname}_guidance.graphml.gz
glueh5adgraphatacrnatype2coembed -e scglue -d "$outdir" -b "$name" -r "$rna" -a "$atac" -g "$graphguidance" -l majorclass -L majorclass -R sampleid -A sampleid -m "$rnarep" -n "$atacrep"


# 3. Plot co-embedding umap
infile=glueh5adgraphatacrnatype2coembed/ON_coembed_10000_X_scVI_Harmony.h5ad
bname=$(basename "$infile" .h5ad)
outdir=scrnah5ad2seuratumapby
scrnah5ad2seuratumapby -e r-seurat_v4.4.0 -d "$outdir" -b "$bname" -s scglue_domain__ -g scglue_celltype__ -g majorclass -L 1 -t F -H 3.5 -W 3.5 -a T -c T -- "$infile"
