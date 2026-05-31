#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import scanpy as sc
import networkx as nx
import scglue
from multiprocess import Pool
import itertools
import seaborn as sns
import pandas as pd
import anndata as ad

# Run co-embedding among low-representations from rna and atac
def cmd(*args):
	rnarep, atacrep=args
	runname=f"{bname}_{rnarep}_{atacrep}"

	# 1. RNA
	if rnafile.endswith('.h5ad'):
		rna=sc.read_h5ad(rnafile)
		rna.obs['scglue_domain__']='RNA'
		rna.obs['scglue_celltype__']=rna.obs[rnalabel]
		rna.obs['scglue_batch__']=rna.obs[rnabatchkey]
	else:
		print(f"Error: format is not supported for {rnafile}", file=sys.stderr)
		sys.exit(-1)
	print('==> rna')
	print(rna)
	
	# 2. ATAC
	if atacfile.endswith('.h5ad'):
		atac=sc.read_h5ad(atacfile)
		atac.obs['scglue_domain__']='ATAC'
		atac.obs['scglue_celltype__']=atac.obs[ataclabel]
		atac.obs['scglue_batch__']=atac.obs[atacbatchkey]
	else:
		print(f"Error: format is not supported for {atacfile}", file=sys.stderr)
		sys.exit(-1)
	print('==> atac')
	print(atac)
	
	# 3. Graph guidance
	guidance=nx.read_graphml(guidefile)
	print('==> guidance')
	print(guidance)
	guidance_hvf=guidance.subgraph(
		itertools.chain(
			rna.var.query('highly_variable').index,
			atac.var.query('highly_variable').index,
			),
		).copy()
	print('==> guidance_hvf')
	print(guidance_hvf)

	if rnarep not in rna.obsm:
		print(f"Warning: {rnarep} is not found in {rnafile}. Skipping {runname}!")
		return

	if atacrep not in atac.obsm:
		print(f"Warning: {atacrep} is not found in {atacfile}. Skipping {runname}!")
		return

	scglue.models.configure_dataset(
		rna,
		'NB',
		use_highly_variable=True,
		use_layer='counts',
		use_rep=rnarep,
		use_batch='scglue_batch__',
		use_cell_type='scglue_celltype__',
	)

	scglue.models.configure_dataset(
		atac,
		'NB',
		use_highly_variable=True,
		use_rep=atacrep,
		use_batch='scglue_batch__',
		use_cell_type='scglue_celltype__',
	)

	# Training
	glue=scglue.models.fit_SCGLUE(
		{'RNA': rna, 'ATAC': atac},
		guidance_hvf,
		init_kws={
			'latent_dim': latent_dim,
			'h_depth': h_depth,
			'random_seed': seed,
			},
		fit_kws={
			'data_batch_size': data_batch_size,
			'graph_batch_size': graph_batch_size,
			'align_burnin': align_burnin,
			'max_epochs': max_epochs,
			'patience': patience,
			'reduce_lr_patience': reduce_lr_patience,
			'directory': f"glue_{runname}",
			},
	)
	print(f"==> glue: {runname}")
	print(glue)

	glue.save(f"{runname}.dill")

	dx=scglue.models.integration_consistency(
		glue,
		{'RNA': rna, 'ATAC': atac},
		guidance_hvf,
	)
	dx.to_csv(f"{runname}_consistency.txt.gz", sep='\t', index=False)

	rna.obsm['X_glue']=glue.encode_data('RNA', rna)
	atac.obsm['X_glue']=glue.encode_data('ATAC', atac)

	combined=ad.concat([rna, atac])
	sc.pp.neighbors(combined, use_rep='X_glue', metric='cosine', random_state=seed)
	sc.tl.umap(combined, random_state=seed)

	sc.write(filename=f"{runname}.h5ad", adata=combined)
	combined.obs.insert(loc=0, column='barcode', value=combined.obs.index)
	combined.obs.to_csv(f"{runname}_obs.txt.gz", sep='\t', index=False)
	combined.var.insert(loc=0, column='symbol', value=combined.var.index)
	combined.var.to_csv(f"{runname}_var.txt.gz", sep='\t', index=False)

	feature_embeddings=glue.encode_graph(guidance_hvf)
	feature_embeddings=pd.DataFrame(feature_embeddings, index=glue.vertices)
	feature_embeddings.insert(loc=0, column='scglue_feature', value=feature_embeddings.index)
	feature_embeddings.to_csv(f"{runname}_feature_embed.txt.gz", sep='\t', index=False)

	# rna.varm['X_glue']=feature_embeddings.reindex(rna.var_names).to_numpy()
	# atac.varm['X_glue']=feature_embeddings.reindex(atac.var_names).to_numpy()

	# RNA
	sc.write(filename=f"{runname}_RNA.h5ad", adata=rna)
	rna.obs.insert(loc=0, column='barcode', value=rna.obs.index)
	rna.obs.to_csv(f"{runname}_RNA_obs.txt.gz", sep='\t', index=False)
	rna.var.insert(loc=0, column='symbol', value=rna.var.index)
	rna.var.to_csv(f"{runname}_RNA_var.txt.gz", sep='\t', index=False)
	# ATAC
	sc.write(filename=f"{runname}_ATAC.h5ad", adata=atac)
	atac.obs.insert(loc=0, column='barcode', value=atac.obs.index)
	atac.obs.to_csv(f"{runname}_ATAC_obs.txt.gz", sep='\t', index=False)
	atac.var.insert(loc=0, column='symbol', value=atac.var.index)
	atac.var.to_csv(f"{runname}_ATAC_var.txt.gz", sep='\t', index=False)
	# HVG Graph guidance
	nx.write_graphml(guidance_hvf, f"{runname}_guidance_hvf.graphml.gz")

	# Plot UMAP
	for splitby in ['scglue_domain__', 'scglue_celltype__', 'scglue_batch__']:
		ncolor=len(combined.obs[splitby].value_counts())
		if ncolor<100:
			sc.pl.umap(combined, color=splitby, frameon=False, show=False, save=f"{runname}_ATAC_umap_{splitby}_wolabel.png")
			sc.pl.umap(combined, color=splitby, frameon=False, show=False, save=f"{runname}_ATAC_umap_{splitby}_ondata.png",
				legend_loc='on data',
				legend_fontsize='xx-small',
				legend_fontweight='normal',
				legend_fontoutline=1,
				)
		else:
			palette=sns.husl_palette(ncolor)
			sc.pl.umap(combined, color=splitby, palette=palette, frameon=False, show=False, save=f"{runname}_ATAC_umap_{splitby}_wolabel.png")
			sc.pl.umap(combined, color=splitby, palette=palette, frameon=False, show=False, save=f"{runname}_ATAC_umap_{splitby}_onlabel.png") # duplicate

if __name__ == "__main__":
	# Run the job in parallel
	with Pool(processes=numthreads) as p:
		p.starmap(cmd, itertools.product(rnareps, atacreps))
