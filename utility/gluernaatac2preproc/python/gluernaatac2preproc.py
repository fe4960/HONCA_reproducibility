#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import scanpy as sc
import networkx as nx
import scglue
import seaborn as sns
import episcanpy as epi

def preproc_rna(file):
	if file.endswith('.h5ad'):
		rna=sc.read_h5ad(file)
	else:
		print(f"Error: format is not supported for {file}", file=sys.stderr)
		sys.exit(-1)

	print('==> rna')
	print(rna)
	print(rna.X)
	print(rna.X.data)

	rna.layers['counts']=rna.X.copy()
	sc.pp.normalize_total(rna)
	sc.pp.log1p(rna)

	# 1. HVG
	## seurat_v3 uses raw counts, while others use normalized counts
	if rnahvgflavor=='seurat_v3':
		sc.pp.highly_variable_genes(
			rna,
			flavor=rnahvgflavor,
			n_top_genes=rnahvgntop,
			subset=True,
			layer='counts',
			batch_key=rnabatchkey,
		)
	else:
		sc.pp.highly_variable_genes(
			rna,
			flavor=rnahvgflavor,
			n_top_genes=rnahvgntop,
			subset=True,
			batch_key=rnabatchkey,
		)

	# 2. UMAP
	sc.pp.scale(rna)
	sc.tl.pca(rna, n_comps=100, svd_solver='auto', random_state=seed)
	sc.pp.neighbors(rna, random_state=seed, metric='cosine')
	sc.tl.umap(rna, random_state=seed)

	# 3. plot UMAP and save .h5ad
	sc.set_figure_params(dpi_save=500, figsize=(5, 5))
	for splitby in rna.obs.columns:
		ncolor=len(rna.obs[splitby].value_counts())
		if ncolor<100:
			sc.pl.umap(rna, color=splitby, frameon=False, show=False, save=f"{bname}_RNA_umap_{splitby}_wolabel.png")
			sc.pl.umap(rna, color=splitby, frameon=False, show=False, save=f"{bname}_RNA_umap_{splitby}_ondata.png",
				legend_loc='on data',
				legend_fontsize='xx-small',
				legend_fontweight='normal',
				legend_fontoutline=1,
				)
		else:
			palette=sns.husl_palette(ncolor)
			sc.pl.umap(rna, color=splitby, palette=palette, frameon=False, show=False, save=f"{bname}_RNA_umap_{splitby}_wolabel.png")
			sc.pl.umap(rna, color=splitby, palette=palette, frameon=False, show=False, save=f"{bname}_RNA_umap_{splitby}_onlabel.png") # duplicate

	# 4. Genomic regions for gene symbols
	scglue.data.get_gene_annotation(rna, gtf=gtffile, gtf_by=gtfby)
	## bug fix symbols without regions from GTF
	feature=rna.var.loc[~rna.var[['chromStart', 'chromEnd']].isna().any(axis=1)].index
	rna=rna[:, feature].copy()

	sc.write(filename=f"{outdir}/{bname}_RNA.h5ad", adata=rna)
	return rna

def preproc_atac(file):
	if file.endswith('.h5ad'):
		atac=sc.read_h5ad(file)
	else:
		print(f"Error: format is not supported for {file}", file=sys.stderr)
		sys.exit(-1)

	print('==> atac')
	print(atac)
	print(atac.X)
	print(atac.X.data)

	# bug fix: do not store counts layer to save memory
	# atac.layers['counts']=atac.X.copy()

	# 1. HVG

	# epi.pp.binarize(atac)
	## Do not binarize peak values but use raw counts instead
	## Ref: 1. https://github.com/gao-lab/GLUE/issues/116
	## Ref: 2. https://www.nature.com/articles/s41592-023-02112-6

	epi.pp.filter_features(atac, min_cells=10)
	epi.pp.cal_var(atac)
	var_annot=atac.var.sort_values(ascending=False, by='variability_score')
	nb_features=int(atac.shape[1]*atacpeakfrac)
	min_score=var_annot['variability_score'][nb_features]
	atac.var['highly_variable']=False
	atac.var.loc[atac.var['variability_score']>min_score, 'highly_variable']=True

	print(f"{nb_features=}, {min_score=}", file=sys.stderr)
	print(atac.var, file=sys.stderr)
	print(atac.var['highly_variable'].value_counts(), file=sys.stderr)

	# 2. LSI + UMAP
	scglue.data.lsi(atac, n_components=100, n_iter=15, use_highly_variable=True)
	sc.pp.neighbors(atac, use_rep='X_lsi', random_state=seed, metric='cosine')
	sc.tl.umap(atac, random_state=seed)

	# 3. plot UMAP and save .h5ad
	for splitby in atac.obs.columns:
		ncolor=len(atac.obs[splitby].value_counts())
		if ncolor<100:
			sc.pl.umap(atac, color=splitby, frameon=False, show=False, save=f"{bname}_ATAC_umap_{splitby}_wolabel.png")
			sc.pl.umap(atac, color=splitby, frameon=False, show=False, save=f"{bname}_ATAC_umap_{splitby}_ondata.png",
				legend_loc='on data',
				legend_fontsize='xx-small',
				legend_fontweight='normal',
				legend_fontoutline=1,
				)
		else:
			palette=sns.husl_palette(ncolor)
			sc.pl.umap(atac, color=splitby, palette=palette, frameon=False, show=False, save=f"{bname}_ATAC_umap_{splitby}_wolabel.png")
			sc.pl.umap(atac, color=splitby, palette=palette, frameon=False, show=False, save=f"{bname}_ATAC_umap_{splitby}_onlabel.png") # duplicate

	# 4. Split peak feature
	atacvarsplit=atac.var_names.str.split(r"[:-]")
	atac.var['chrom']=atacvarsplit.map(lambda x: x[0])
	atac.var['chromStart']=atacvarsplit.map(lambda x: x[1]).astype(int)
	atac.var['chromEnd']=atacvarsplit.map(lambda x: x[2]).astype(int)

	# 5. Subset the HVF to speedup downstream calculation
	atac_hvf=atac[:, atac.var['highly_variable']].copy()

	# 6. atac -> bname_full.h5ad, atac_hvf -> bname.h5ad
	sc.write(filename=f"{bname}_ATAC_full.h5ad", adata=atac)
	sc.write(filename=f"{bname}_ATAC.h5ad", adata=atac_hvf)
	return atac_hvf

def preproc_graphguidance(rna, atac):
	# 1. Graph guidance
	guidance=scglue.genomics.rna_anchored_guidance_graph(rna, atac)
	scglue.graph.check_graph(guidance, [rna, atac])

	# 2. Save objects
	nx.write_graphml(guidance, f"{bname}_guidance.graphml.gz")

if __name__ == "__main__":
	# 1. RNA
	if refreshrna:
		rna=preproc_rna(rnafile)
	else:
		rna=sc.read_h5ad(rnafile)

	# 2. ATAC
	if refreshatac:
		atac=preproc_atac(atacfile)
	else:
		atac=sc.read_h5ad(atacfile)

	# 3. Graph guidance
	if refreshgraph:
		preproc_graphguidance(rna, atac)
