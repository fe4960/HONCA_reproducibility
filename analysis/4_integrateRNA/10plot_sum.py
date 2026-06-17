import scanpy as sc
dir="/dfs3b/ruic20_lab/junw42/"
adata_query1=sc.read_h5ad(f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_normCount.h5ad")

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)
sc.settings.figdir = f"{dir}/HCA_ON/data/4_scvi/RNA/scvi/figures"
bname="ONONH_majorclass"
adata_query1.obs["age1"]=adata_query1.obs["age"]
adata_query1.obs["age1"]=adata_query1.obs["age1"].cat.add_categories(["0"])
adata_query1.obs.loc[adata_query1.obs["age"]=="1day","age1"]="0"
adata_query1.obs["age1"]=adata_query1.obs["age1"].cat.remove_unused_categories()
t=["tissue","race","age1","sampleid","gender"]

adata_query1.obs["age1"]=adata_query1.obs["age1"].astype("int").copy()
for i in t:
    sc.pl.embedding(adata_query1, basis="X_umap", color=[i], ncols=1,frameon=False,save=f'{bname}_{i}_wolabel.png', palette="tab20")


markers={"B cell": ["BANK1", "BLK", "BCL11A"], "Schwann cell": ["SCN7A","INSC","L1CAM","NGFR"], "T cell": ["BCL11B", "THEMIS", "IL7R", "CAMK4", "CD2"], "Astrocyte": ["AC092957.1", "DCLK1", "PAX8", "GFAP", "ALDH1L1"], "Endothelial cell": ["PTPRB", "LDB2", "VWF", "FLT1", "ANO2", "MECOM","PECAM1"], "Fibroblast": ["CLMP","BICC1","SCARA5","FHL2","DCN"], "Macrophage": ["CD163", "F13A1","MRC1", "STAB1"], "Melanocyte": ["ABCB5","PAX3","SYNPR","PMEL","GALNTL6"], "Microglia": ["MPZ","MLIP","GALNT17","PRX","CD74","CX3CR1"], "NK cell": ["KLRD1","CD247","TXK","KLRF1"], "Oligodendrocyte": ["ST18","CERCAM","TMEM144","MOBP","SH3GL3","MOG"], "Oligodendrocyte precursor cell": ["CSMD1","TNR","OPCML","SNTG1","CSMD3","OLIG2","OLIG1"], "Mular pericyte": ["RGS5","ABCC9","TRPC4"], "RPE": ["CCBE1","SLC38A11","ATP6V1C2","OTX2","BEST1","RPE65"], "Mular vascular associated smooth muscle cell": ["RGS6","MYH11","ACTA2","UNC13C"] , "Rod": "PDE6A", "Cone": "ARR3", "BC": "VSX2", "AC": "PAX6", "RGC": "RBPMS", "HC": "ONECUT1", "MG": "RLBP1"}

markers1 = dict(sorted(markers.items()))

sc.pl.dotplot(adata_query1, markers1, groupby='majorclass',  dendrogram=False, save=f"{bname}_NSForest_marker.png", use_raw=False)




known_markers={"B cell": ["DOCK2", "PTPRC", "EBF1","BCL11A"], "T cell": ["DOCK2","PTPRC","CD2","THEMIS"], "Astrocyte": ["ALDH1L1", "GFAP"], "Endothelial cell": ["FLT1","PECAM1"], "Fibroblast": ["COL1A2","DCN"], "Macrophage": ["DOCK2","PTPRC","MRC1","CD163"], "Microglia": ["DOCK2", "PTPRC","P2RY12", "ITGAM", "CD74", "CX3CR1"], "NK cell": ["DOCK2", "PTPRC", "KLRD1","KLRF1"], "Oligodendrocyte": ["MOG","MOBP"], "Oligodendrocyte precursor cell": ["OLIG1","OLIG2"], "Mular pericyte": ["TAGLN","ACTA2","TRPC4","RGS5","LAMA2"], "Mular vacular associated smooth muscle cell": ["TAGLN","ACTA2","UNC13C","RGS6","CASQ2","FGF7","SLC7A2","CLSTN2"],"Schwann cell": ["SCN7A","NCAM1","L1CAM","NGFR"], "Schwann cell myelinating": ["MPZ","PRX","MBP"], "RPE": ["MLANA","TYR","BEST1","RPE65"], "Melanocyte": ["PMEL","PAX3"], "Rod": "PDE6A", "Cone": "ARR3", "BC": "VSX2", "AC": "PAX6", "RGC": "RBPMS", "HC": "ONECUT1", "MG": "RLBP1" }

known_markers1 = dict(sorted(known_markers.items()))

sc.pl.dotplot(adata_query1, known_markers1, groupby='majorclass',  dendrogram=False, save=f"{bname}_known_marker.png", use_raw=False)

indi_markers={"AC": "PAX6", "Astrocyte": "GFAP", "B cell": "EBF1", "BC": "VSX2", "Cone": "ARR3", "Endothelial cell": "PECAM1", "Fibroblast": ["COL1A2"], "HC": "ONECUT1", "MG": "RLBP1", "Macrophage": "CD163", "Melanocyte": "PAX3", "Microglia": "CD74", "Mural cell": "ACTA2", "NK/T cell": "THEMIS", "Oligodendrocyte" : "MOBP", "Oligodendrocyte precursor cell": "CSMD1", "RGC": "RBPMS", "RPE": "BEST1", "Rod" : "PDE6A", "Schwann cell": "SCN7A"}

indi_markers1 = dict(sorted(indi_markers.items()))

sc.pl.dotplot(adata_query1, indi_markers1, groupby='majorclass',  dendrogram=False, save=f"{bname}_indi_marker.png", use_raw=False)


adata_query1.obs["source"]="Chen"
adata_query1.obs.loc[adata_query1.obs["sampleid"].str.contains("GSM",regex=False),"source"]="Sanes"

i="source"

sc.pl.embedding(adata_query1, basis="X_umap", color=[i], ncols=1,frameon=False,save=f'{bname}_{i}_wolabel.png', palette="tab20")


import numpy as np

df=adata_query1.obs.groupby(["majorclass","source"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat.to_csv(f'{dir}/HCA_ON/data/4_scvi/RNA/scvi/{bname}_source_conf_mat_table.csv',header=True, index=True)

import pandas as pd

pnas=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/3_ref/joshSanes_meta/SCP_meta.csv",header=0)
samplelist=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONH_ON_public_meta.csv'
meta=pd.read_csv(samplelist,header=0,sep=",")
dt=meta[["sample","sampleid"]].set_index("sample").to_dict()["sampleid"]
pnas["bc"]=pnas["NAME"].replace(dt, regex=True)


#adata_query.obs["majorclass_sanes"]="chen"
pnas=pnas.set_index("bc")
idx=pnas.index.intersection(adata_query1.obs.index)

adata_query2=adata_query1[idx]

adata_query2.obs["majorclass_sanes"]=pnas.loc[idx,"cell_type__ontology_label"]

df = adata_query2.obs.groupby(["majorclass_sanes", "majorclass"]).size().unstack(fill_value=0)
norm_df = df / df.sum(axis=0)

import matplotlib.pyplot as plt

fig=plt.figure(figsize=(8, 8))
_ = plt.pcolor(norm_df)
_ = plt.xticks(np.arange(0.5, len(df.columns), 1), df.columns, rotation=90)
_ = plt.yticks(np.arange(0.5, len(df.index), 1), df.index)
plt.xlabel("majorclass_sanes")
plt.ylabel("majorclass")
fig.savefig(f'{bname}_{i}_sanes_current_conf.png')


adata_query1.obs.to_csv(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/4_scvi/RNA/scvi/ONONH_majorclass_normCount_detailed.obs.gz")

#row_max=conf_mat.max(axis=1)

#conf_mat.idxmax(axis=1).to_csv(f"{bname}_leiden_joshOnly_scvi_conf_mat_table_max.csv")

#sc.pl.dotplot(adata_query1.cat.reorder_categories(["B cell", "Schwann cell", "T cell", "astrocyte", "endothelia cell", "fibroblast", "macrophage","melanocyte","microglia", "natural killer cell", "oligodendrocyte","oligodendrocyte precursor cell", "pericyte", "pigmented epithelial cell", "vascular associated smooth muscle cell", "Rod", "Cone", "BC", "AC", "RGC", "HC", "MG"], inplace=True), markers, groupby='leiden', dendrogram=False, save=f"{bname}_NSForest_marker.png", use_raw=False)
#sc.pl.dotplot(adata_query1, markers, groupby='majorclass', categories_order=["B cell", "Schwann cell ", "NK/T cell", "Astrocyte", "Endothelial cell", "Fibroblast", "Macrophage","Melanocyte","Microglia", "Oligodendrocyte","Oligodendrocyte precursor cell", "Mural cell", "RPE",  "Rod", "Cone", "BC", "AC", "RGC", "HC", "MG"], dendrogram=False, save=f"{bname}_NSForest_marker.png", use_raw=False)
#sc.pl.dotplot(adata_query1, known_markers,  groupby='leiden', categories_order=["B cell", "Schwann cell ", "NK/T cell", "Astrocyte", "Endothelial cell", "Fibroblast", "Macrophage","Melanocyte","Microglia", "Oligodendrocyte","Oligodendrocyte precursor cell", "Mural cell", "RPE",  "Rod", "Cone", "BC", "AC", "RGC", "HC", "MG"]  ,  dendrogram=False, save=f"{bname}_known_marker.png", use_raw=False)
