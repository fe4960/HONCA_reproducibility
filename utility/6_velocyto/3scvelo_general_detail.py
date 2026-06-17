import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import sys
#scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

indir=sys.argv[1]
h5ad=sys.argv[2]
dir1=f"{indir}/veloc"
bname=sys.argv[3]
hvg=int(sys.argv[4])
import os
if not os.path.exists(dir1):
    os.makedirs(dir1)

#adata=sc.read(f"{indir}/{h5ad}.h5ad")


#adata.X=adata.layers["counts"]

#ldata = sc.read(f"{indir}/{h5ad}_sp_loom.h5ad")

#adata = scv.utils.merge(adata, ldata)

adata=sc.read(f'{dir1}/{bname}_velocity.h5ad')

scv.pl.velocity_graph(adata, threshold=.1, save=f'{dir1}/{bname}_velocity_graph.png', color="subclass")

df = adata.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]
df.to_csv(f'{dir1}/{bname}_velocity_gene.csv')

scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95], save=f'{dir1}/{bname}_cell_cycle.png')

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95],save=f'{dir1}/{bname}_velocity_length_confidence.png')

#df = adata.obs.groupby('subclass')[keys].mean().T
#df.style.background_gradient(cmap='coolwarm', axis=1, save=f'{dir1}/{bname}_velocity_length_confidence_table.png')

x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell=70)
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax, save=f'{dir1}/{bname}_cell_transitions.png')

# this is needed due to a current bug - bugfix is coming soon.
adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='subclass')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5,save=f'{dir1}/{bname}_paga_transition.png')

adata.write(f'{dir1}/{bname}_velocity.h5ad')
