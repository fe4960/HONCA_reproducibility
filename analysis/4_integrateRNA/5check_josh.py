import pandas as pd
import scanpy as sc
import numpy as np
adata_query=sc.read_h5ad("/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONONH_hvg2k_epoch20_new_scvi.h5ad")

pnas=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/3_ref/joshSanes_meta/SCP_meta.csv",header=0)


samplelist=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONH_ON_public_meta.csv'
meta=pd.read_csv(samplelist,header=0,sep=",")
dt=meta[["sample","sampleid"]].set_index("sample").to_dict()["sampleid"]
pnas["bc"]=pnas["NAME"].replace(dt, regex=True)


adata_query.obs["celltype"]="chen"
pnas=pnas.set_index("bc")
idx=pnas.index.intersection(adata_query.obs.index)

adata_query.obs.loc[idx,"celltype"]=pnas.loc[idx,"cell_type__ontology_label"]


sc.pl.embedding(adata_query, basis="X_umap", color=["celltype"],legend_loc="on data", ncols=1,frameon=False,save=f'_{bname}_joshCellType_onData.png', palette="tab20")

df=adata_query.obs.groupby(["leiden","celltype"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat.to_csv(f'{bname}_leiden_josh_scvi_conf_mat_table.csv',header=True, index=True)
row_max=conf_mat.max(axis=1)
dominate_name=row_max[row_max > 0.7].index

df=adata_query[adata_query.obs["celltype"]!="chen"].obs.groupby(["leiden","celltype"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat.to_csv(f'{bname}_leiden_joshOnly_scvi_conf_mat_table.csv',header=True, index=True)
row_max=conf_mat.max(axis=1)

conf_mat.idxmax(axis=1).to_csv(f"{bname}_leiden_joshOnly_scvi_conf_mat_table_max.csv")


adata_query1=adata_query[adata_query.obs["leiden"].isin(["17","7","36","30","42","22","23","26","28","29","33","35","43","3","31"])]
#adata_query1=adata_query[adata_query.obs["leiden"].isin(["17","7","36","30","42","22","23","26","28","29","33","35","43","3","31"])].write(f"{outdir}/{bname}_unknown_for_retina.h5ad")

adata_query1.X=adata_query1.layers["counts"]
adata_query1.write(f"{outdir}/{bname}_unknown_for_retina_raw.h5ad")
