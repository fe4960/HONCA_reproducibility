import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

rcParams["figure.figsize"] = (5,5)
celltype="Oligodendrocyte_precursor_cell"
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"
mk="sanes_mk"
#adata=sc.read_h5ad(f"{outdir}/{celltype}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad")
adata1=sc.read_h5ad(f"{outdir}/{celltype}_hvg10000_epochnone_seurat_rs_0.3_clean_sb_none_scvi_trg.h5ad")
adata=sc.read_h5ad(f"{outdir}/{celltype}_subclass.h5ad")
bname=f"{celltype}_subclass_seurat_cycling_rmRPE"
#bname=f"{celltype}_subclass_seurat_cycling"
adata.obs["subclass"]="Oligodendrocyte_precursor_cell"
sc.settings.figdir = f"{outdir}/figures"

obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_hvg10000_epochnone_seurat_rs_0.8_clean_sb_none.obs.gz",header=0, index_col=0)
idx=obs.loc[obs["leiden"]==10].index

adata=adata[adata1.obs.index]

del adata.obsm["X_scVI"]
adata.obsm["X_scVI"]=adata1.obsm["X_scVI"]
del adata.obsm["X_umap"]
adata.obsm["X_umap"]=adata1.obsm["X_umap"]
del adata.uns["umap"]
adata.uns["umap"]=adata1.uns["umap"]
del adata.obs["leiden"]
adata.obs["leiden"]=adata1.obs["leiden"]


#obs=pd.read_csv(f"{outdir}/{celltype}_res_1_clean_res_1.0.obs.gz",header=0, index_col=0)
#adata.obs.loc[obs.index,"leiden1"]=obs["leiden"]
#obs=pd.read_csv(f"{outdir}/{celltype}_res_1_clean_res_1.0.obs.gz",header=0, index_col=0)
#idx=obs[obs["leiden"]==13].index
#adata=adata[adata.obs["leiden"]!="6"]
#adata.obs["subclass"]="Macrophage"
#idx=obs[obs["leiden"]==9].index
#adata.obs.loc[idx,"subclass"]="Macrophage_VCAN+"
#idx=obs[obs["leiden"]==10].index
#adata.obs.loc[idx,"subclass"]="Macrophage_RORA+"
#idx=obs[obs["leiden"]==11].index
#adata.obs.loc[idx,"subclass"]="Macrophage_RRM2+"
#idx=obs[obs["leiden"]==12].index
#adata.obs.loc[idx,"subclass"]="Macrophage_JUN+"

adata=adata[(adata.obs["sampleid"]!="MMD_23_17976_ONH_RNA")&(adata.obs["sampleid"]!="MMD_23_20181_ONH_RNA")]

celltype1="OPC"

dt={"0": f"{celltype1}", "2": f"{celltype1}_GLIS3+", "1": f"MFOL", "3": f"COP_NFOL", "4": f"NFOL_UTRN+"}
adata.obs["subclass"]=adata.obs["leiden"].replace(dt)
adata.obs["subclass"]=adata.obs["subclass"].cat.add_categories(["OPC_Cycling"])
#adata.obs.loc[adata.obs["subclass"].isin(["0"]),"subclass"]=f"{celltype}_ON_1"
#adata.obs.loc[adata.obs["subclass"].isin(["2"]),"subclass"]=f"{celltype}_ON_2"
#adata.obs.loc[adata.obs["subclass"].isin(["3"]),"subclass"]=f"{celltype}_ON_3"
#adata.obs.loc[adata.obs["subclass"].isin(["4"]),"subclass"]=f"{celltype}_ON_4"
#adata.obs.loc[adata.obs["subclass"].isin(["1"]),"subclass"]=f"{celltype}_ONH"
adata.obs.loc[idx,"subclass"]="OPC_Cycling"

#adata.obs.loc[idx,"subclass"]="Astrocyte_retina"

sc.pl.embedding(adata, basis="X_umap", color=["subclass"],legend_loc="on data", ncols=1,frameon=False,save=f'{bname}_onData.png', palette="tab20")
sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,save=f'{bname}.png', palette="tab20")

markers={"B cell": ["BANK1", "BLK", "BCL11A"], "Schwann cell": ["SCN7A","INSC","L1CAM","NGFR"], "T cell": ["BCL11B", "THEMIS", "IL7R", "CAMK4", "CD2"], "Astrocyte": ["AC092957.1", "DCLK1", "PAX8", "GFAP", "ALDH1L1"], "Endothelial cell": ["PTPRB", "LDB2", "VWF", "FLT1", "ANO2", "MECOM","PECAM1"], "Fibroblast": ["CLMP","BICC1","SCARA5","FHL2","DCN"], "Macrophage": ["CD163", "F13A1","MRC1", "STAB1"], "Melanocyte": ["ABCB5","PAX3","SYNPR","PMEL","GALNTL6"], "Microglia": ["MPZ","MLIP","GALNT17","PRX","CD74","CX3CR1"], "NK cell": ["KLRD1","CD247","TXK","KLRF1"], "Oligodendrocyte": ["ST18","CERCAM","TMEM144","MOBP","SH3GL3","MOG"], "Oligodendrocyte precursor cell": ["CSMD1","TNR","OPCML","SNTG1","CSMD3","OLIG2","OLIG1"], "Mular pericyte": ["RGS5","ABCC9","TRPC4"], "RPE": ["CCBE1","SLC38A11","ATP6V1C2","OTX2","BEST1","RPE65"], "Mular vascular associated smooth muscle cell": ["RGS6","MYH11","ACTA2","UNC13C"] , "Rod": "PDE6A", "Cone": "ARR3", "BC": "VSX2", "AC": "PAX6", "RGC": "RBPMS", "HC": "ONECUT1", "MG": "RLBP1"}

markers1 = dict(sorted(markers.items()))

known_markers={"B cell": ["DOCK2", "PTPRC", "EBF1","BCL11A"], "T cell": ["DOCK2","PTPRC","CD2","THEMIS"], "Astrocyte": ["ALDH1L1", "GFAP"], "Endothelial cell": ["FLT1","PECAM1"], "Fibroblast": ["COL1A2","DCN"], "Macrophage": ["DOCK2","PTPRC","MRC1","CD163"], "Microglia": ["DOCK2", "PTPRC","P2RY12", "ITGAM", "CD74", "CX3CR1"], "NK cell": ["DOCK2", "PTPRC", "KLRD1","KLRF1"], "Oligodendrocyte": ["MOG","MOBP"], "Oligodendrocyte precursor cell": ["OLIG1","OLIG2"], "Mular pericyte": ["TAGLN","ACTA2","TRPC4","RGS5","LAMA2"], "Mular vacular associated smooth muscle cell": ["TAGLN","ACTA2","UNC13C","RGS6","CASQ2","FGF7","SLC7A2","CLSTN2"],"Schwann cell": ["SCN7A","NCAM1","L1CAM","NGFR"], "Schwann cell myelinating": ["MPZ","PRX","MBP"], "RPE": ["MLANA","TYR","BEST1","RPE65"], "Melanocyte": ["PMEL","PAX3"], "Rod": "PDE6A", "Cone": "ARR3", "BC": "VSX2", "AC": "PAX6", "RGC": "RBPMS", "HC": "ONECUT1", "MG": "RLBP1" }

known_markers1 = dict(sorted(known_markers.items()))

sc.pl.dotplot(adata, markers1, groupby='subclass', dendrogram=False, save=f"{bname}_marker.png", use_raw=False)
sc.pl.dotplot(adata, known_markers1, groupby='subclass', dendrogram=False, save=f"{bname}_known_marker.png", use_raw=False)

sc.tl.dendrogram(adata, groupby='subclass', use_rep="X_scVI", use_raw=False, n_pcs=30)

adata.obs[["nCount_RNA","subclass","sampleid"]].groupby(["subclass","sampleid"]).count().to_csv(f"{outdir}/{bname}_sample_number.txt",sep="\t")

lis=pd.read_csv(f'{indir}/{mk}',header=None)
dt=lis[0].values #.set_index(0).to_dict()[1]

sc.pl.dotplot(adata, dt, groupby='subclass', dendrogram=True, save=f"{bname}_marker_sanes_ct.png", use_raw=False)

mk={"other":["GAP43","LAMA2","EBF1","FOXP1"],"OPC":["PDGFRA","CSPG4","OLIG2"],"COP":["GPR17","VCAN","NEU4","BMP4","FYN","SOX6","PTPRZ1","PDCD4"],"NFOL":["TCF7L2","CEMIP2"],"MFOL":["MBP","OPALIN","CTPS1"],"MOL":["MBP","PLP1","MOG","GRM3","KLK6","PTGDS","ANXA5","HOPX"]}


sc.pl.dotplot(adata, dt, groupby='subclass', dendrogram=True, save=f"{bname}_marker_other_ct.png", use_raw=False)


sc.pl.violin(adata,keys=["pct_counts_mt","n_genes_by_counts", "total_counts"], groupby="subclass", rotation=75, save=f"{bname}_QC.png", multi_panel=True)

plot=["race", "gender", "age_year"]

for pl in plot:
    sc.pl.embedding(adata, basis="X_umap", color=[pl],ncols=1,frameon=False,save=f'_{bname}_{pl}.png', palette="tab20")


#adata.obs.set_levels(["source"], ascending=False, inplace=True)
dt={"Chen": "lightblue", "Sanes": "blue"}
sc.pl.embedding(adata, basis="X_umap", color=["source"],ncols=1,frameon=False,save=f'_{bname}_source.png', palette=dt)

dt={"ON": "lightblue", "ONH": "blue"}
sc.pl.embedding(adata, basis="X_umap", color=["tissue"],ncols=1,frameon=False,save=f'_{bname}_tissue.png', palette=dt)

rcParams["figure.figsize"] = (14,5)

#total=adata.obs[["nCount_RNA","subclass"]].groupby(["subclass"]).count()
adata.obs[["nCount_RNA","source","subclass"]].groupby(["subclass","source"]).count().to_csv(f"{outdir}/{bname}_source.csv")
df=pd.read_csv(f"{outdir}/{bname}_source.csv")
total=df.groupby(["subclass"])["nCount_RNA"].transform("sum")
df["nCount_RNA"]=df["nCount_RNA"]/total

plt.barh( df[df["source"]!="Chen"]["subclass"], df[df["source"]!="Chen"]["nCount_RNA"], color='blue' )
plt.barh(df[df["source"]=="Chen"]["subclass"], df[df["source"]=="Chen"]["nCount_RNA"],  left=df[df["source"]!="Chen"]["nCount_RNA"], color='lightblue' )
plt.legend(["Sanes", "Chen"])
plt.xlabel("Cell proportion")
plt.savefig(f"{outdir}/figures/{bname}_source.png")

#top_bar = mpatches.Patch(color='darkblue', label='Sanes')
#bottom_bar = mpatches.Patch(color='lightblue', label='Chen')
#plt.legend(handles=[top_bar, bottom_bar])

#sns_plot=sns.barplot(data=df, y="subclass", x="nCount_RNA", hue="source", dodge=True)
#plt.xlabel("Cell proportion")
#plt.savefig(f"{outdir}/{bname}_source.png")

#bar1 = sns.barplot(x="nCount_RNA",  y="subclass", data=total, color='darkblue')
#bar2 = sns.barplot(x="nCount_RNA", y="subclass", data=df[df["source"]!="Chen"], color='lightblue')

# add legend
#top_bar = mpatches.Patch(color='darkblue', label='Chen')
#bottom_bar = mpatches.Patch(color='lightblue', label='Sanes')
# show the graph
#plt.show()



df=adata.obs[["nCount_RNA","tissue","subclass"]].groupby(["subclass","tissue"]).count().to_csv(f"{outdir}/{bname}_tissue.csv")
df=pd.read_csv(f"{outdir}/{bname}_tissue.csv")
total=df.groupby(["subclass"])["nCount_RNA"].transform("sum")
df["nCount_RNA"]=df["nCount_RNA"]/total

plt.barh( df[df["tissue"]!="ON"]["subclass"], df[df["tissue"]!="ON"]["nCount_RNA"], color='blue' )
plt.barh(df[df["tissue"]=="ON"]["subclass"], df[df["tissue"]=="ON"]["nCount_RNA"],  left=df[df["tissue"]!="ON"]["nCount_RNA"], color='lightblue' )
plt.legend(["ONH", "ON"])
plt.xlabel("Cell proportion")
plt.savefig(f"{outdir}/figures/{bname}_tissue.png")



sc.tl.rank_genes_groups(adata, "subclass", method="wilcoxon", use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=10,  save=f'{bname}_TRG10.png', dendrogram=False)
adata.write(f"{outdir}/{bname}.h5ad")

