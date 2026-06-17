import scanpy as sc
import pandas as pd
from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

dir="/dfs3b/ruic20_lab/junw42/"

sc.settings.figdir = f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/figures"
adata_query=sc.read_h5ad(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_hvg2k_epoch20_new_scvi.h5ad")
ct=pd.read_csv(f"{dir}/HCA_ON/scripts/4_integrateRNA/8leiden_celltype",sep="\t",header=None)
dt=ct.set_index(0).to_dict()[1]
bname="ONONH_hvg2k_epoch20_new"
#adata_query1=adata_query[adata_query.obs["leiden"].isin(["14","21","41","45"])]
adata_query.obs["majorclass"]=adata_query.obs["leiden"].astype(int).replace(dt,regex=True)

#######
adata_query1=sc.read_h5ad(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/cellclass/MG/MG_hvg2k_epoch20_scvi.h5ad")
idx=adata_query1[adata_query1.obs["leiden"].isin(["4","9"])].obs.index
adata_query.obs.loc[idx,"majorclass"]="Astrocyte"
sc.pl.embedding(adata_query, basis="X_umap", color=["majorclass"],legend_loc="on data", ncols=1,frameon=False,save=f'{bname}_celltype_beforeClean_onData.png', palette="tab20")
sc.pl.embedding(adata_query, basis="X_umap", color=["majorclass"], ncols=1,frameon=False,save=f'{bname}_celltype_beforeClean_wolabel.png', palette="tab20")


########
adata_query1=adata_query[-adata_query.obs["leiden"].isin(["14","21","41","45","48"])]
sc.pl.embedding(adata_query1, basis="X_umap", color=["majorclass"],legend_loc="on data", ncols=1,frameon=False,save=f'{bname}_celltype_v1_onData.png', palette="tab20")
sc.pl.embedding(adata_query1, basis="X_umap", color=["majorclass"], ncols=1,frameon=False,save=f'{bname}_celltype_v1_wolabel.png', palette="tab20")
del adata_query1.raw
adata_query1.write(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass.h5ad")

adata_query2=adata_query1.copy()

del adata_query2.raw
del adata_query2.layers["counts"]
adata_query2.write(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_normCount.h5ad")


adata_query2.obs.to_csv(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_normCount.obs.gz",sep="\t")
adata_query2.var.to_csv(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_normCount.var.gz",sep="\t")

adata_query1.X=adata_query1.layers["counts"]
del adata_query1.layers["counts"]
del adata_query1.raw

adata_query1.write(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_rawCount.h5ad")

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt



adata_query1=sc.read_h5ad(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_rawCount.h5ad")

df=pd.DataFrame(adata_query1.obs["majorclass"].value_counts())
df["major class"]=df.index
df["cell number"]=df["majorclass"]
df = df.sort_values(by='cell number', ascending=False)
fig=plt.figure(figsize=(22, 10))
sns.barplot(x='cell number', y="major class", data=df, order=df["major class"],palette="tab20")
plt.show()
fig.savefig(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/figures/ONONH_majorclass.png", dpi=300)



