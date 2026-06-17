import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

rcParams["figure.figsize"] = (5,5)
rcParams.update({'font.size': 14})

celltype="Fibroblast"
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"
mk="sanes_mk"
adata=sc.read_h5ad(f"{outdir}/{celltype}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad")
#bname=f"{celltype}_subclass_clean_new"
bname=f"{celltype}_subclass_clean_new_rename"

#bname=f"{celltype}_subclass"
adata.obs["subclass"]="Fibroblast"
sc.settings.figdir = f"{outdir}/figures"

adata1=sc.read_h5ad(f"{outdir}/{celltype}_hvg10000_epochnone_seurat_rs_1_clean_sb_none_seed_7_scvi_trg.h5ad")


#adata1=sc.read_h5ad(f"{outdir}/{celltype}_hvg10000_epochnone_seurat_v3_rs_0.6_clean_rm4_6_13_sb_scvi_trg.h5ad")
adata=adata[adata1.obs.index]

del adata.obsm["X_scVI"]
adata.obsm["X_scVI"]=adata1.obsm["X_scVI"]
del adata.obsm["X_umap"]
adata.obsm["X_umap"]=adata1.obsm["X_umap"]
del adata.uns["umap"]
adata.uns["umap"]=adata1.uns["umap"]
obs=pd.read_csv(f"{outdir}/{celltype}_hvg10000_epochnone_seurat_rs_0.8_clean_sb_none_seed_7.obs.gz",header=0,index_col=0)
adata.obs["leiden"]=obs.loc[adata.obs.index,"leiden"].astype("str")
obs1=pd.read_csv(f"{outdir}/{celltype}_hvg10000_epochnone_seurat_rs_0.7_clean_sb_none_seed_7.obs.gz",header=0,index_col=0)
#lc=adata1[adata1.obs["leiden"]=="5"].obs.index
ara_1=obs1.loc[obs1["leiden"]==7].index

celltype1="Fibro"
#dt={"1": f"{celltype1}_arachnoid", "6": f"{celltype1}_arachnoid",   "0": f"{celltype1}_pia", "3": f"{celltype1}_sclera", "2": f"{celltype1}_dura", "4": f"{celltype1}_RPEchoroid", "5" : f"{celltype1}_x"}

dt={"0":"Fibro_arachnoid_PLCB4+" ,"12":"Fibro_arachnoid_PLCB4+" ,"13": "Fibro_arachnoid_PLCB4+","11":"Fibro_dura","2":"Fibro_dura","6":"Fibro_RPEchoroid_SMOC2+","8":"Fibro_perivascular","1":"Fibro_pia_ABCA10+","10":"Fibro_pia_ABCA10+","7": "Fibro_HMCN1+","4":"Fibro_RPEchoroid_BMP5+","3":"Fibro_sclera_NOX4+","5":"Fibro_sclera_KCNMA1+","9":"Fibro_x"}


#dt={"0":"Fibro_arachnoid_PLCB4+" ,"12":"Fibro_arachnoid_PLCB4+" ,"13": "Fibro_arachnoid_PLCB4+","11":"Fibro_dura","2":"Fibro_dura","6":"Fibro_RPEchoroid_SMOC2+","8":"Fibro_perivascular","1":"Fibro_pia_ABCA10+","10":"Fibro_pia_ABCA10+","7": "Fibro_pia_HMCN1+_LC?","4":"Fibro_RPEchoroid_BMP5+","3":"Fibro_sclera_NOX4+","5":"Fibro_sclera_KCNMA1+","9":"Fibro_x"}
adata.obs["subclass"]=adata.obs["leiden"].replace(dt)
#adata.obs.loc[adata.obs["subclass"].isin(["0"]),"subclass"]=f"{celltype}_ON_1"
#adata.obs.loc[adata.obs["subclass"].isin(["2"]),"subclass"]=f"{celltype}_ON_2"
#adata.obs.loc[adata.obs["subclass"].isin(["3"]),"subclass"]=f"{celltype}_ON_3"
#adata.obs.loc[adata.obs["subclass"].isin(["4"]),"subclass"]=f"{celltype}_ON_4"
#adata.obs.loc[adata.obs["subclass"].isin(["1"]),"subclass"]=f"{celltype}_ONH"

adata.obs.loc[ara_1,"subclass"]="Fibro_arachnoid_STXBP6+"
##########adata.obs["subclass"]=adata.obs["subclass"].cat.add_categories(["Fibro_arachnoid_TRPM3_hi","Fibro_arachnoid_STK39_hi","Fibro_pia_ABCA10_hi",f"{celltype1}_pia_CASC15_hi",f"{celltype1}_dura_1",f"{celltype1}_dura_2",f"{celltype1}_lamina_cribrosa"])

#adata.obs.loc[ara_1,"subclass"]=f"{celltype1}_arachnoid_TRPM3_hi"
#adata.obs.loc[ara_2,"subclass"]=f"{celltype1}_arachnoid_STK39_hi"

#adata.obs.loc[pia_1,"subclass"]=f"{celltype1}_pia_ABCA10_hi"
#adata.obs.loc[pia_2,"subclass"]=f"{celltype1}_pia_CASC15_hi"


#adata.obs.loc[dura_1,"subclass"]=f"{celltype1}_dura_1"
#adata.obs.loc[dura_2,"subclass"]=f"{celltype1}_dura_2"


#adata.obs.loc[lc,"subclass"]=f"{celltype1}_lamina_cribrosa"

#adata=adata[-adata.obs["subclass"].isin([f"{celltype1}_arachnoid", f"{celltype1}_dura"])]
#adata.obs["cluster_sub"]=data.obs.subclass.cat.remove_unused_categories()

sc.pl.embedding(adata, basis="X_umap", color=["subclass"],legend_loc="on data", ncols=1,frameon=False,save=f'{bname}_onData.png', palette="tab20")
sc.pl.embedding(adata, basis="X_umap", color=["subclass"],  ncols=1,frameon=False,save=f'{bname}.png', palette="tab20")

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
sc.settings.figdir = f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/figures"

mk=["ACTA2"]#,"GFAP"]
sc.pl.dotplot(adata,mk,groupby="subclass", save="Fibroblast_subclass_LC6_ACTA2.png")

mk=["FN1","TNXB","TNR","VTN","THBS1","THBS2","THBS4","COL1A1","ACTA2","GFAP"]
sc.pl.dotplot(adata,mk,groupby="subclass",save="_Fibroblast_subclass_clean_new_LC_mk.png")
