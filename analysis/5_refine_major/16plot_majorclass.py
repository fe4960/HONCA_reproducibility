import pandas as pd
import scanpy as sc
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

rcParams["figure.figsize"] = (5,5)

dir="/dfs3b/ruic20_lab/junw42"
f1="hvg2000_nonsb_epochnone_seurat_v3_rs"
bname=f"_major_{f1}_1_clean_ON_ONH_new"
adata=sc.read(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/major_ON_ONH_new_final.h5ad")
df=pd.DataFrame(adata.obs["majorclass1"].value_counts())
df["major class"]=df.index
df["cell number"]=df["majorclass1"]
df["major class"]=df["major class"].replace("Oligodendrocyte_precursor_cell","OPC")
df["major class"]=df["major class"].replace("Unknown?","Immune_neuron_unk")



df = df.sort_values(by='cell number', ascending=False)
fig=plt.figure(figsize=(10, 8))
sns.barplot(x='cell number', y="major class", data=df, order=df["major class"],palette="tab20")
plt.show()
fig.savefig(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass.png", dpi=300)

df.to_csv(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass.csv" )

sc.pl.embedding(adata, basis="X_umap", color=["majorclass1"], ncols=1,frameon=False, save=f"_major_{f1}_1_clean_ON_ONH_new_majorclass_onData.png", palette="tab20", legend_loc="on data")
#sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,palette="tab20", save=f"_major_{f1}_1_clean_ON_ONH_new_subclass.png")

sc.settings.figdir = f"{dir}/HCA_ON/data/5_refine_major/scvi/major/figures"

plot=["race",  "age_year", "sampleid"]

for pl in plot:
    sc.pl.embedding(adata, basis="X_umap", color=[pl],ncols=1,frameon=False,save=f'_{bname}_{pl}.png', palette="tab20")


dt={"Male": "lightblue", "Female": "blue"}
sc.pl.embedding(adata, basis="X_umap", color=["gender"],ncols=1,frameon=False,save=f'_{bname}_gender.png', palette=dt)


#adata.obs.set_levels(["source"], ascending=False, inplace=True)
dt={"Chen": "lightblue", "Sanes": "blue"}
sc.pl.embedding(adata, basis="X_umap", color=["source"],ncols=1,frameon=False,save=f'_{bname}_source.png', palette=dt)

dt={"ON": "lightblue", "ONH": "blue"}
sc.pl.embedding(adata, basis="X_umap", color=["tissue"],ncols=1,frameon=False,save=f'_{bname}_tissue.png', palette=dt)



#adata.obs.groupby(['nCount_RNA','tissue','sampleid','majorclass1']).count()

adata.obs[["sampleid","tissue"]].drop_duplicates().groupby(["tissue"]).count()
adata.obs[["age_year","donor","tissue"]].drop_duplicates().groupby(["tissue"]).mean()
adata.obs[["age_year","donor","tissue"]].drop_duplicates().groupby(["tissue"]).std()
adata.obs[["age_year","donor"]].drop_duplicates().std()
adata.obs[["source","donor","tissue"]].drop_duplicates().groupby(["tissue","source"]).count()


df=adata.obs[["nCount_RNA","majorclass1","sampleid","tissue","source"]].groupby(["tissue","source","sampleid","majorclass1"]).count()

df.to_csv(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass_proportion_sample.csv")

obs=pd.read_csv(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/major_ON_ONH_new_final.obs.gz",header=0,sep="\t",index_col=0)
obs["majorclass1"]=obs["majorclass1"].replace("Oligodendrocyte_precursor_cell","OPC")
obs["majorclass1"]=obs["majorclass1"].replace("Unknown?","Immune_neuron_unk")

s=["Chen","Sanes"]
t=["ON","ONH"]
for s1 in s:
    for t1 in t:
        df=obs.loc[(obs["source"]==s1)&(obs["tissue"]==t1),["nCount_RNA","majorclass1","sampleid","tissue","source"]].groupby(["tissue","source","sampleid","majorclass1"]).count()
        df.to_csv(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass_proportion_sample_{t1}_{s1}.csv")

df.to_csv(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass_proportion_sample.csv")

df=adata.obs[["nCount_RNA","majorclass1","sampleid","tissue","source"]].groupby(["tissue","source","sampleid","majorclass1"]).count()

df.to_csv(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass_proportion_sample.csv")



import numpy as np

df=adata.obs.groupby(["majorclass","source"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat.to_csv(f'{dir}/HCA_ON/data/5_refine_major/scvi/major/{bname}_source_conf_mat_table.csv',header=True, index=True)

import pandas as pd

pnas=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/3_ref/joshSanes_meta/SCP_meta.csv",header=0)
samplelist=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONH_ON_public_meta.csv'
meta=pd.read_csv(samplelist,header=0,sep=",")
dt=meta[["sample","sampleid"]].set_index("sample").to_dict()["sampleid"]
pnas["bc"]=pnas["NAME"].replace(dt, regex=True)


#adata_query.obs["majorclass_sanes"]="chen"
pnas=pnas.set_index("bc")
idx=pnas.index.intersection(adata.obs.index)

adata_query2=adata[idx]

adata_query2.obs["majorclass_sanes"]=pnas.loc[idx,"cell_type__ontology_label"]

df = adata_query2.obs.groupby(["majorclass_sanes", "majorclass1"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)

import matplotlib.pyplot as plt

fig=plt.figure(figsize=(15, 15))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("majorclass_sanes")
plt.ylabel("majorclass")
fig.savefig(f'{dir}/HCA_ON/data/5_refine_major/scvi/major/{bname}_sanes_current_conf.png')



