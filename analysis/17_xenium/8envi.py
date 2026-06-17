import os

os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]="-1" # Change to -1 if you want to use CPU!

import warnings
warnings.filterwarnings('ignore')

import scenvi

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd
import scanpy as sc
import colorcet
import umap.umap_ as umap
import sys
import anndata


#sc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.h5ad"
h5ad=sys.argv[1]
output=sys.argv[2]
target_tissue=sys.argv[3]
label=sys.argv[4]
cl=sys.argv[5]
query_list=sys.argv[6]
sc_data = sc.read_h5ad(h5ad)


sample=[]

with open(query_list, "r") as sl:
    for line in sl:
        line1=line.strip()
        sample.append(line1)

print(f'done reading list')
adata_list=[]
i=0
for sam in sample:
    i=i+1
    print(f'read {sam}')
    path=f'{sam}'
#    if os.path.exists(path):
    adata=sc.read(path)
#    else :
#        path=f'{sam}'
#        adata=sc.read(path)
#    adata.obs["sample"]=sam
#    print(f'mody index')
#    if sam not in array:
    adata.obs.index=str(i)+"_"+adata.obs.index
    xs=f"{output}/{target_tissue}/xenium_spatial_{label}_{i}.csv"
    coords = pd.read_csv(xs, index_col=0)
    adata.obsm["spatial"] = coords.values
    adata.obs["batch"]=f"{label}_{i}"
#    adata.obs["system"]=syst
#    adata.obs[cl]="unknown"
#    adata.obs["ident"]="query"
    adata_list.append(adata)
    print(f'done reading {sam}')

st_data=anndata.concat(adata_list)


########
st_data.layers['count'] = st_data.X.copy()

sc.pp.normalize_total(st_data)
sc.pp.log1p(st_data)

st_data.raw=st_data

st_data.layers['log'] = st_data.X.copy() #sc.pp.log1p(sc_data, copy=True).X

st_data.X=st_data.layers['count'].copy()

#st_data.write(f"{output}/{target_tissue}/xenium_spatial_{label}_concat.h5ad")

fit = umap.UMAP(
    n_neighbors = 100,
    min_dist = 0.8,
    n_components = 2,
)


sc_data.layers['count'] = sc_data.X.copy()

sc.pp.normalize_total(sc_data)
sc.pp.log1p(sc_data)

sc_data.raw=sc_data

sc_data.layers['log'] = sc_data.X.copy() #sc.pp.log1p(sc_data, copy=True).X

sc_data.X=sc_data.layers['count'].copy()

gene=st_data.var.index.intersection(sc_data.var.index)



sc.pp.highly_variable_genes(sc_data, layer = 'log', n_top_genes = 2048, batch_key="sampleid")

sc_data.var.loc[gene,'highly_variable']=True


# log1p handles sparse correctly
#X_log = sc_data[:, sc_data.var['highly_variable']].X
#X_log = np.log1p(X_log.A) if not isinstance(X_log, np.ndarray) else np.log1p(X_log)

sc_data1=sc_data[:, sc_data.var['highly_variable']].copy()
sc_data1

#sc_data.obsm['UMAP_exp'] =fit.fit_transform(sc.pp.log1p(sc_data1, copy=True).X)
#sc_data.obsm['UMAP_exp'] =fit.fit_transform(sc_data1.layers['log'])
sc_data.obsm['UMAP_exp'] = sc_data.obsm['X_umap']

# run UMAP
#sc_data.obsm['UMAP_exp'] = fit.fit_transform(X_log)

#sc_data.obsm['UMAP_exp'] = fit.fit_transform(np.log(sc_data[:, sc_data.var['highly_variable']].X + 1))

#sc_data.obsm['UMAP_exp'] = fit.fit_transform(np.log(sc_data[:, sc_data.var['highly_variable']].X + 1))

fig = plt.figure(figsize = (10,10))
sns.scatterplot(x = sc_data.obsm['UMAP_exp'][:, 0], y = sc_data.obsm['UMAP_exp'][:, 1],  hue = sc_data.obs[cl], s = 16,
                palette = "tab20", legend = True)
plt.tight_layout()
plt.axis('off')
plt.title('scRNA-seq Data')
plt.show()
plt.savefig(f"{output}/{target_tissue}/xenium_spatial_{label}_ref_scRNA.png")


st_data.write(f"{output}/{target_tissue}/xenium_spatial_{label}_st_data.h5ad")
sc_data.write(f"{output}/{target_tissue}/xenium_spatial_{label}_sc_data.h5ad")

#########

