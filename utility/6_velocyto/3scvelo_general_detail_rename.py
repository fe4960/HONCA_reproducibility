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

dt={"COP_LAMA2+":"NFOL_UTRN+"}
adata.obs["subclass1"]=adata.obs["subclass"].replace(dt)
scv.tl.paga(adata, groups='subclass1')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5,save=f'{dir1}/{bname}_paga_transition1.png')

adata.write(f'{dir1}/{bname}_velocity1.h5ad')
