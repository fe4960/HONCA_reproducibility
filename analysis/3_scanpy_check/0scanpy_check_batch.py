import matplotlib.pyplot as plt
#import seaborn as sns
#import scvelo as scv
import numpy as np
import pandas as pd
import anndata
from os.path import exists
#import plotly.express as px
import scanpy as sc
#import scvi
import sys
sys.setrecursionlimit(10000)


query_list=sys.argv[1]
celltype=sys.argv[2]
queryname=sys.argv[3]
in_dir=sys.argv[4]
batch_key="sample"
sample=[]
############
with open(query_list, "r") as sl:
	for line in sl:
		line1=line.strip()
		sample.append(line1)

print(f'done reading list')
adata_list=[]


for sam in sample:
	print(f'read {sam}')
	path=f'{in_dir}/{sam}.h5ad'
	adata=sc.read(path)
	adata.obs["sample"]=sam
	print(f'mody index')
	adata.obs.index=sam+"_"+adata.obs.index
	adata_list.append(adata)
	print(f'done reading {sam}')

dataset1=anndata.concat(adata_list)


if celltype == "major":
	adata_query=dataset1
else:
	adata_query=dataset1[(dataset1.obs[label] == celltype )]



labels_key = ""
adata_query.layers["counts"] = adata_query.X.copy()
sc.pp.normalize_total(adata_query,target_sum=1e4)
sc.pp.log1p(adata_query)
adata_query.raw = adata_query

sc.pp.highly_variable_genes(
#		adata_query, flavor="seurat_v3", n_top_genes=2000, subset=True
		adata_query, n_top_genes=2000, subset=True,span=1, batch_key=batch_key)

sc.tl.pca(adata_query)

sc.pp.neighbors(adata_query)
sc.tl.umap(adata_query)
sc.tl.leiden(adata_query, resolution=1)

#print(adata_query.obs.leiden) 
author_cell_type = adata_query.obs.leiden.astype(str).astype("category")
sc.pl.embedding(adata_query, basis="X_umap", color=["sample","leiden"],ncols=1,frameon=False,save=f'_{queryname}_{celltype}_scanpy_merge.png', palette="tab20")
#adata_query.obs["sample1"]=adata_query.obs["sample"]
for sample in adata_query.obs["sample"].keys():
#   adata_query.obs.at[adata_query.obs["sample"]!=sample,"sample1"]="other"
    adata=adata_query[adata_query.obs["sample"]==sample]
    sc.pl.embedding(adata, basis="X_umap", color=["sample"],ncols=1,frameon=False,save=f'_{queryname}_{sample}_scanpy_merge.png', palette="tab20")
