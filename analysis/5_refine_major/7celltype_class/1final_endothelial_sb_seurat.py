import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

rcParams["figure.figsize"] = (5,5)
celltype="Endothelial_cell"
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"
mk="sanes_mk"
adata=sc.read_h5ad(f"{outdir}/{celltype}_hvg2000_epochnone_seurat_rs_0.4_clean_rm_scvi_trg.h5ad")

sc.settings.figdir=f"{indir}/clean/figures"
#adata1=sc.read_h5ad(f"{outdir}/{celltype}_hvg10000_epochnone_seurat_v3_rs_1_span_1_subclass_scvi_trg.h5ad")

adata1=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_hvg10000_epochnone_seurat_rs_0.6_clean_rm_sb_rm_4_6_none_scvi_trg.h5ad")

idx=adata.obs.index.intersection(adata1.obs.index)
adata=adata[idx]


c=["X_scVI","X_umap"]

for c1 in c:
    del adata.obsm[c1]
    adata.obsm[c1]=adata1[idx].obsm[c1]

c=['hvg', 'leiden', 'leiden_colors', 'neighbors', 'umap']

for c1 in c:
    del adata.uns[c1]
    adata.uns[c1]=adata1[idx].uns[c1]


bname=f"{celltype}_subclass_sb_seurat_rmRPE"

######
obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_hvg10000_epochnone_seurat_rs_0.5_clean_rm_sb_rm_4_6_none_scvi_trg.h5ad.obs.gz", header=0, index_col=0)

adata.obs["leiden1"]=obs.loc[adata.obs.index,"leiden"].astype(str)


obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_hvg10000_epochnone_seurat_rs_1_clean_rm_sb_rm_4_6_none_scvi_trg.h5ad.obs.gz", header=0, index_col=0)

idx=obs.loc[obs["leiden"].isin([3,9])].index
adata.obs.loc[idx,"leiden1"]="leiden"+"_"+obs.loc[idx,"leiden"].astype("str")

adata.obs.loc[adata.obs["leiden1"]=="3","leiden1"]="leiden_3"
#####

adata.obs["subclass"]="Endo"
adata=adata[adata.obs["leiden1"]!="5"].copy()
dt={"0": "Capillary", "1": "Vein_Venule", "2": "Artery_Arteriole","leiden_9": "Venule_postcapillary", "4": "Capillary_fenestrated", "leiden_3": "Venule_POSTN+_PLVAP+"}
adata.obs["subclass"]=adata.obs["leiden1"].replace(dt)
#adata.obs.loc[idx1,"subclass"]="Venule_POSTN+"


adata=adata[(adata.obs["sampleid"]!="MMD_23_17976_ONH_RNA")&(adata.obs["sampleid"]!="MMD_23_20181_ONH_RNA")]

sc.pl.embedding(adata, basis="X_umap", color=["subclass"],legend_loc="on data", ncols=1,frameon=False,save=f'{bname}_onData.png', palette="tab20")

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

