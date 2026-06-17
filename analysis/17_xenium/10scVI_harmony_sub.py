import matplotlib.pyplot as plt
#import seaborn as sns
#import scvelo as scv
import numpy as np
import pandas as pd
import anndata
import os
from os.path import exists
#import plotly.express as px
import scanpy as sc
import scvi
import sys

sys.setrecursionlimit(10000)

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)



ref=sys.argv[1]
query=sys.argv[2]
syst=sys.argv[3]
outdir=sys.argv[4]
label=sys.argv[5]
cl=sys.argv[6]
sc.settings.figdir = f"{outdir}/figures"

adata_ref = sc.read(ref)

adata_ref.obs["system"]="snRNA"
adata_ref.obs["ident"]="ref"


adata_query=sc.read(query)
adata_query.obs[cl]="unknown"
adata_query=adata_query[(adata_query.obs["lr_celltype"]==label)&(adata_query.obs["system"]=="xenium")].copy()
adata_query.X=adata_query.layers["counts"]

gene=adata_query.var.index.intersection(adata_ref.var.index)

adata_query=adata_query[:,gene].copy()
adata_ref=adata_ref[:,gene].copy()

adata=anndata.concat([adata_query,adata_ref])


adata.layers["counts"] = adata.X.copy()


sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)



scvi.settings.seed = 7
scvi.model.SCVI.setup_anndata(adata, batch_key = "sampleid", layer="counts")

vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="zinb")

vae.train( accelerator="gpu")

adata.obsm["X_scVI"] = vae.get_latent_representation()


sc.external.pp.harmony_integrate(adata, key="system", basis='X_scVI', adjusted_basis='X_scVI_harmony')

sc.pp.neighbors(adata, random_state=42, use_rep='X_scVI_harmony')
sc.tl.umap(adata, random_state=42)

sc.pl.umap(adata, color="system", title="System", save=f"_{label}_pca_sys_harmony.png" , frameon=False )

sc.pl.umap(adata, color="majorclass", title="majorclass", save=f"_{label}_pca_ct_harmony.png", frameon=False )
sc.pl.umap(adata, color="sampleid", title="sampleid", save=f"_{label}_pca_ct_sampleid.png", frameon=False )

adata.write(f"{outdir}/sn_st_scvi_harmony.h5ad")


