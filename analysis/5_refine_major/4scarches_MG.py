import scanpy as sc
import numpy as np
import pandas as pd

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/cellclass/MG/"
bname="MG_hvg2k_epoch20_scarches"
adata_query=sc.read_h5ad("/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/query_query_latent.h5ad")
ct=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/4_integrateRNA/7celltype",sep="\t",header=None)
dt=ct.set_index(0).to_dict()[1]
adata_query.obs["majorclass"]=adata_query.obs["scANVI_predictions"].replace(dt)


adata_query1=sc.read_h5ad(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/cellclass/MG/MG_hvg2k_epoch20_scvi.h5ad")

#idx=adata_query.obs.index.intersection(adata_query1.obs.index)
#adata_query1.obs.loc[idx,"retina_class"]=adata_query.obs.loc[idx,"majorclass"]
#df=adata_query1.obs.groupby(["leiden","retina_class"]).size().unstack(fill_value=0)
#conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
#conf_mat.to_csv(f'{outdir}/{bname}_leiden_retina_class_conf_mat_table.csv',header=True, index=True)
#row_max=conf_mat.max(axis=1)
#dominate_name=row_max[row_max > 0.7].index
#conf_mat.idxmax(axis=1)
#conf_mat.idxmax(axis=1).to_csv(f"{outdir}/{bname}_leiden_retina_class_conf_mat_table_max.csv")


idx=adata_query[adata_query.obs["scANVI_prediction_max_probability"]>0.9].obs.index.intersection(adata_query1.obs.index)
adata_query1.obs.loc[idx,"retina_class"]=adata_query.obs.loc[idx,"majorclass"]
df=adata_query1.obs.groupby(["leiden","retina_class"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat.to_csv(f'{outdir}/{bname}_leiden_retina_class_conf_mat_table_scANVI_h0.9.csv',header=True, index=True)
row_max=conf_mat.max(axis=1)
dominate_name=row_max[row_max > 0.7].index
conf_mat.idxmax(axis=1)
conf_mat.idxmax(axis=1).to_csv(f"{outdir}/{bname}_leiden_retina_class_conf_mat_table_max_scANVI_h0.9.csv")
