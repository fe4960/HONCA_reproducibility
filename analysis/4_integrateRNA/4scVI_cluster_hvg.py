import matplotlib.pyplot as plt
#import seaborn as sns
#import scvelo as scv
import numpy as np
import pandas as pd
import anndata
from os.path import exists
#import plotly.express as px
import scanpy as sc
import scvi
import sys

sys.setrecursionlimit(10000)

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

query=sys.argv[1]
celltype=sys.argv[2]
batch_key=sys.argv[3]
bname=sys.argv[4]
label_key=sys.argv[5]
outdir=sys.argv[6]
label=sys.argv[7]
hvg=int(sys.argv[8])
sb=sys.argv[9]
nlat=int(sys.argv[10])
norm=sys.argv[11]
obs=sys.argv[12]
item=sys.argv[13]
max_epoch=sys.argv[14]
rm=sys.argv[15]
scvi.settings.seed = int(sys.argv[16])
ntrg=int(sys.argv[17])
mk=sys.argv[18]
indir=sys.argv[19]
res=float(sys.argv[20])
detail=sys.argv[21]
flavor=sys.argv[22]
dataset1=sc.read(query)

sc.settings.figdir = f"{outdir}/figures"

if label != "none":
    dataset1=dataset1[dataset1.obs[label] != "unassigned"]

if celltype == "major":
	adata_query=dataset1
elif celltype == "NN":
	adata_query=dataset1[dataset1.obs[label].isin(["MG","RPE","Microglia","Astrocyte"])].copy()
else:
	adata_query=dataset1[dataset1.obs[label] == celltype].copy()

if rm != "none":
    lis=pd.read_csv(f'{indir}/{rm}',header=None)
    rmcluster=lis[0].values.astype(str)
    adata_query=adata_query[-adata_query.obs["leiden"].isin(rmcluster)]

df=adata_query.obs[batch_key].value_counts()
ssample=df.loc[df==1].index.values

adata_query=adata_query[-adata_query.obs[batch_key].isin(ssample)]

if obs != "none":
    df=pd.read_csv(obs,header=0,index_col=0)
    idx=adata_query.obs.index.intersection(df.index)
    key=item.split(",")
    for k in key:
        adata_query.obs.loc[idx,k]=df.loc[idx,k]

adata_query.obs[batch_key] = adata_query.obs[batch_key].astype("str").astype("category")

if label_key != "none":
    adata_query.obs[label_key] = adata_query.obs[label_key].astype("str").astype("category")

#adata_query.X=adata_query.layers["counts"].copy()
if norm == "t":
    adata_query.layers["counts"] = adata_query.X.copy()
    sc.pp.normalize_total(adata_query,target_sum=1e4)
    sc.pp.log1p(adata_query)
#adata_query.raw = adata_query
if flavor == "seurat_v3":
    sc.pp.highly_variable_genes(
		adata_query,
		flavor="seurat_v3",
		n_top_genes=hvg,
		batch_key = batch_key,
		subset=False,
                layer="counts",
		span=1
	)
else:
#    try: 
    sc.pp.highly_variable_genes(adata_query,
#		flavor="seurat_v3",
                n_top_genes=hvg,
		batch_key = batch_key,
		subset=False,
		span=1 )

adata_query.var["highly_variable"].to_csv(f'{outdir}/{bname}_hvg.csv')
#except:
#	print("An exception occurred! Fix the batch or run without batch now.")
#	sc.pp.highly_variable_genes(
#		adata_query, flavor="seurat_v3", n_top_genes=2000, subset=True
#		adata_query, n_top_genes=hvg, subset=False, span=1
#	)

#scvi.settings.seed = 7
if label_key == "none":
	scvi.model.SCVI.setup_anndata(adata_query, batch_key = batch_key, layer="counts")
#	scvi.data.setup_anndata(adata_query, batch_key = batch_key)
else:
	scvi.model.SCVI.setup_anndata(adata_query, batch_key = batch_key, layer="counts" , labels_key = label_key)

vae = scvi.model.SCVI(adata_query, n_layers=2, n_latent=nlat, gene_likelihood="zinb")

if max_epoch != "none":
    vae.train(max_epochs=int(max_epoch), accelerator="gpu")
else :
    vae.train( accelerator="gpu")

adata_query.obsm["X_scVI"] = vae.get_latent_representation()
#sc.pp.neighbors(adata_query, use_rep = "X_scVI")
#sc.tl.leiden(adata_query, resolution=res)
#sc.tl.umap(adata_query)
adata_query.write(f'{outdir}/{bname}_scvi.h5ad')

