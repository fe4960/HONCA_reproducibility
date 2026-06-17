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

adata=sc.read(f'{dir1}/{bname}_velocity.h5ad')

#scv.tl.velocity(adata, mode="dynamical", use_latent_time=True)
#scv.tl.velocity_graph(adata)
#scv.pl.velocity_embedding_stream(adata, basis='umap', save=f'{dir1}/{bname}_velocity_latent.png', color="subclass")
#adata.write(f'{dir1}/{bname}_velocity_latent.h5ad')

adata.obs["root_cells0"]=adata.obs["root_cells"]

adata.obs["end_points0"]=adata.obs["end_points"]


adata.obs["root_cells"]=1-adata.obs["root_cells0"]
adata.obs["end_points"]=1-adata.obs["end_points0"]

adata.obs["latent_time0"]=adata.obs["latent_time"]

adata.obs["latent_time"]=1-adata.obs["latent_time0"]
scv.tl.velocity(adata, mode="dynamical", use_latent_time=True)
scv.tl.velocity_graph(adata,tkey="latent_time")
scv.pl.velocity_embedding_stream(adata, basis='umap', save=f'{dir1}/{bname}_velocity_reservse.png', color="subclass")
adata.write(f'{dir1}/{bname}_velocity_reservse.h5ad')

