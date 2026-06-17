import matplotlib.pyplot as plt
#import seaborn as sns
#import scvelo as scv
import numpy as np
import pandas as pd
import anndata
from os.path import exists
#import plotly.express as px
import scanpy as sc
import scvi
import sys

sys.setrecursionlimit(10000)

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

query=sys.argv[1]
celltype=sys.argv[2]
batch_key=sys.argv[3]
bname=sys.argv[4]
label_key=sys.argv[5]
outdir=sys.argv[6]
label=sys.argv[7]
hvg=int(sys.argv[8])
sb=sys.argv[9]
nlat=int(sys.argv[10])
norm=sys.argv[11]
obs=sys.argv[12]
item=sys.argv[13]
max_epoch=sys.argv[14]
rm=sys.argv[15]
scvi.settings.seed = int(sys.argv[16])
ntrg=int(sys.argv[17])
mk=sys.argv[18]
indir=sys.argv[19]
res=float(sys.argv[20])
detail=sys.argv[21]
flavor=sys.argv[22]
gl=sys.argv[23]
span=float(sys.argv[24])
dirrm=sys.argv[25]
rank=sys.argv[26]
kp=sys.argv[27]
dataset1=sc.read(query)

sc.settings.figdir = f"{outdir}/figures"

if label != "none":
    dataset1=dataset1[dataset1.obs[label] != "unassigned"]

if celltype == "major":
	adata_query=dataset1
elif celltype == "NN":
	adata_query=dataset1[dataset1.obs[label].isin(["MG","RPE","Microglia","Astrocyte"])].copy()
else:
	adata_query=dataset1[dataset1.obs[label] == celltype].copy()

if obs != "none":
    df=pd.read_csv(obs,header=0,index_col=0)
    idx=adata_query.obs.index.intersection(df.index)
    adata_query=adata_query[idx]
    key=item.split(",")
    for k in key:
        adata_query.obs[k]=df[k]
        adata_query.obs.loc[idx,k]=df.loc[idx,k]
    adata_query.obs["leiden"]=adata_query.obs["leiden"].astype(str)



if rm != "none":
    lis=pd.read_csv(f'{dirrm}/{rm}',header=None)
    rmcluster=lis[0].values.astype(str)
    adata_query=adata_query[-adata_query.obs["majorclass"].isin(rmcluster)]

if kp != "none":
    lis=pd.read_csv(f'{dirrm}/{kp}',header=None)
    kpcluster=lis[0].values.astype(str)
    print(kpcluster)
    adata_query=adata_query[adata_query.obs["leiden"].isin(kpcluster)]



df=adata_query.obs[batch_key].value_counts()
ssample=df.loc[df==1].index.values

adata_query=adata_query[-adata_query.obs[batch_key].isin(ssample)]


adata_query.obs[batch_key] = adata_query.obs[batch_key].astype("str").astype("category")

if label_key != "none":
    adata_query.obs[label_key] = adata_query.obs[label_key].astype("str").astype("category")

adata_query.raw=adata_query

#adata_query.X=adata_query.layers["counts"].copy()
if norm == "t":
    adata_query.layers["counts"] = adata_query.X.copy()
    sc.pp.normalize_total(adata_query,target_sum=1e4)
    sc.pp.log1p(adata_query)
#adata_query.raw = adata_query
if flavor == "seurat_v3":
    print(f"{flavor}   seurat_v3")
    sc.pp.highly_variable_genes(
		adata_query,
		flavor="seurat_v3",
		n_top_genes=hvg,
		batch_key = batch_key,
		subset=True,
                layer="counts",
		span=span
	)
elif flavor == "seurat" :

    print(f"{flavor}   seurat")
    sc.pp.highly_variable_genes(adata_query,
                n_top_genes=hvg,
		batch_key = batch_key,
		subset=True,
		span=span )
#quit()
#except:
#	print("An exception occurred! Fix the batch or run without batch now.")
#	sc.pp.highly_variable_genes(
#		adata_query, flavor="seurat_v3", n_top_genes=2000, subset=True
#		adata_query, n_top_genes=hvg, subset=False, span=1
#	)

#scvi.settings.seed = 7
if label_key == "none":
	scvi.model.SCVI.setup_anndata(adata_query, batch_key = batch_key, layer="counts")
#	scvi.data.setup_anndata(adata_query, batch_key = batch_key)
else:
	scvi.model.SCVI.setup_anndata(adata_query, batch_key = batch_key, layer="counts" , labels_key = label_key)

vae = scvi.model.SCVI(adata_query, n_layers=2, n_latent=nlat, gene_likelihood=gl)

if max_epoch != "none":
    vae.train(max_epochs=int(max_epoch), accelerator="gpu")
else :
    vae.train( accelerator="gpu")


adata_query.obsm["X_scVI"] = vae.get_latent_representation()

#vae.get_latent_representation().to_csv(f'{outdir}/{bname}_X_scVI.csv')

#quit()

sc.pp.neighbors(adata_query, use_rep = "X_scVI")
sc.tl.leiden(adata_query, resolution=res)
sc.tl.umap(adata_query)
#adata_query.write(f'{outdir}/{bname}_scvi.h5ad')
if detail == "t":
    adata_query.obs["source"]="Chen"
    adata_query.obs.loc[adata_query.obs["sampleid"].str.contains("GSM",regex=False),"source"]="Sanes"
    adata_query.obs["age_year"]=adata_query.obs["age"]
    adata_query.obs["age_year"]=adata_query.obs["age_year"].cat.add_categories(["0"])
    adata_query.obs.loc[adata_query.obs["age"]=="1day","age_year"]="0"
    adata_query.obs["age_year"]=adata_query.obs["age_year"].cat.remove_unused_categories()
    adata_query.obs["age_year"]=adata_query.obs["age_year"].astype("int").copy()



    plot=[batch_key, "leiden","tissue","race", "gender", "source", "age_year"]

    for pl in plot:
        sc.pl.embedding(adata_query, basis="X_umap", color=[pl],ncols=1,frameon=False,save=f'_{bname}_{pl}.png', palette="tab20")




adata_query.var["mt"] = adata_query.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata_query, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

sc.pl.embedding(adata_query, basis="X_umap", color=["leiden"],legend_loc="on data", ncols=1,frameon=False,save=f'_{bname}_leiden_onData.png', palette="tab20")


sc.pl.violin(adata_query,keys=["MALAT1", "pct_counts_mt","n_genes_by_counts", "total_counts"], groupby="leiden", rotation=75, save=f"_{bname}_QC.png", multi_panel=True)

#markers={"B cell": ["BANK1", "BLK", "BCL11A"], "Schwann cell": ["SCN7A","INSC","L1CAM","NGFR"], "T cell": ["BCL11B", "THEMIS", "IL7R", "CAMK4", "CD2"], "astrocyte": ["AC092957.1", "DCLK1", "PAX8", "GFAP", "ALDH1L1"], "endothelium": ["PTPRB", "LDB2", "VWF", "FLT1", "ANO2", "MECOM","PECAM1"], "fibroblast": ["CLMP","BICC1","SCARA5","FHL2","DCN"], "macrophage": ["CD163", "F13A1","MRC1", "STAB1"], "melanocyte": ["ABCB5","PAX3","SYNPR","PMEL","GALNTL6"], "microglial cell": ["MPZ","MLIP","GALNT17","PRX"], "natural killer cell": ["KLRD1","CD247","TXK","KLRF1"], "oligodendrocyte": ["ST18","CERCAM","TMEM144","MOBP","SH3GL3","MOG"], "oligodendrocyte precursor cell": ["CSMD1","TNR","OPCML","SNTG1","CSMD3","OLIG2","OLIG1"], "pericyte": ["RGS5","ABCC9","TRPC4"], "pigmented epithelial cell": ["CCBE1","SLC38A11","ATP6V1C2","OTX2","BEST1","RPE65"], "vascular associated smooth muscle cell": ["RGS6","MYH11","ACTA2","UNC13C"] }

#known_markers={"B cell": ["DOCK2", "PTPRC", "EBF1","BCL11A"], "T cell": ["DOCK2","PTPRC","CD2","THEMIS"], "astrocyte": ["ALDH1L1", "GFAP"], "endothelium": ["FLT1","PECAM1"], "fibroblast": ["COL1A2","DCN"], "macrophage": ["DOCK2","PTPRC","MRC1","CD163"], "microglia cell": ["DOCK2", "PTPRC","P2RY12", "ITGAM"], "natural killer cell": ["DOCK2", "PTPRC", "KLRD1","KLRF1"], "oligodendrocyte": ["MOG","MOBP"], "oligodendrocyte precursor cell": ["OLIG1","OLIG2"], "pericyte": ["TAGLN","ACTA2","TRPC4","RGS5","LAMA2"], "vacular associated smooth muscle cell": ["TAGLN","ACTA2","UNC13C","RGS6","CASQ2","FGF7","SLC7A2","CLSTN2"],"Schwann cell": ["SCN7A","NCAM1","L1CAM","NGFR"], "myelinating Schwann cell": ["MPZ","PRX","MBP"], "pigmented epithelial cell": ["MLANA","TYR","BEST1","RPE65"], "melanocyte": ["PMEL","PAX3"] }
markers={"B cell": ["BANK1", "BLK", "BCL11A"], "Schwann cell": ["SCN7A","INSC","L1CAM","NGFR"], "T cell": ["BCL11B", "THEMIS", "IL7R", "CAMK4", "CD2"], "Astrocyte": ["AC092957.1", "DCLK1", "PAX8", "GFAP", "ALDH1L1"], "Endothelial cell": ["PTPRB", "LDB2", "VWF", "FLT1", "ANO2", "MECOM","PECAM1"], "Fibroblast": ["CLMP","BICC1","SCARA5","FHL2","DCN"], "Macrophage": ["CD163", "F13A1","MRC1", "STAB1"], "Melanocyte": ["ABCB5","PAX3","SYNPR","PMEL","GALNTL6"], "Microglia": ["MPZ","MLIP","GALNT17","PRX","CD74","CX3CR1"], "NK cell": ["KLRD1","CD247","TXK","KLRF1"], "Oligodendrocyte": ["ST18","CERCAM","TMEM144","MOBP","SH3GL3","MOG"], "Oligodendrocyte precursor cell": ["CSMD1","TNR","OPCML","SNTG1","CSMD3","OLIG2","OLIG1"], "Mular pericyte": ["RGS5","ABCC9","TRPC4"], "RPE": ["CCBE1","SLC38A11","ATP6V1C2","OTX2","BEST1","RPE65"], "Mular vascular associated smooth muscle cell": ["RGS6","MYH11","ACTA2","UNC13C"] , "Rod": "PDE6A", "Cone": "ARR3", "BC": "VSX2", "AC": "PAX6", "RGC": "RBPMS", "HC": "ONECUT1", "MG": "RLBP1"}

markers1 = dict(sorted(markers.items()))

known_markers={"B cell": ["DOCK2", "PTPRC", "EBF1","BCL11A"], "T cell": ["DOCK2","PTPRC","CD2","THEMIS"], "Astrocyte": ["ALDH1L1", "GFAP"], "Endothelial cell": ["FLT1","PECAM1"], "Fibroblast": ["COL1A2","DCN"], "Macrophage": ["DOCK2","PTPRC","MRC1","CD163"], "Microglia": ["DOCK2", "PTPRC","P2RY12", "ITGAM", "CD74", "CX3CR1"], "NK cell": ["DOCK2", "PTPRC", "KLRD1","KLRF1"], "Oligodendrocyte": ["MOG","MOBP"], "Oligodendrocyte precursor cell": ["OLIG1","OLIG2"], "Mular pericyte": ["TAGLN","ACTA2","TRPC4","RGS5","LAMA2"], "Mular vacular associated smooth muscle cell": ["TAGLN","ACTA2","UNC13C","RGS6","CASQ2","FGF7","SLC7A2","CLSTN2"],"Schwann cell": ["SCN7A","NCAM1","L1CAM","NGFR"], "Schwann cell myelinating": ["MPZ","PRX","MBP"], "RPE": ["MLANA","TYR","BEST1","RPE65"], "Melanocyte": ["PMEL","PAX3"], "Rod": "PDE6A", "Cone": "ARR3", "BC": "VSX2", "AC": "PAX6", "RGC": "RBPMS", "HC": "ONECUT1", "MG": "RLBP1" }

known_markers1 = dict(sorted(known_markers.items()))

sc.pl.dotplot(adata_query, markers1, groupby='leiden', dendrogram=False, save=f"_{bname}_marker.png", use_raw=True)
sc.pl.dotplot(adata_query, known_markers1, groupby='leiden', dendrogram=False, save=f"_{bname}_known_marker.png", use_raw=True)



#adata_query.obs[["reactionid","leiden","sampleid"]].groupby(["leiden","sampleid"]).count().to_csv(f"{outdir}/{bname}_sample_number.txt",sep="\t")
adata_query.obs[["nCount_RNA","leiden","sampleid"]].groupby(["leiden","sampleid"]).count().to_csv(f"{outdir}/{bname}_sample_number.txt",sep="\t")

adata_query.write(f'{outdir}/{bname}_scvi_trg.h5ad')


if rank != "none":
    sc.tl.rank_genes_groups(adata_query, "leiden", method="wilcoxon", use_raw=True)
    sc.pl.rank_genes_groups(adata_query, n_genes=ntrg,sharey=False, save=f'_{bname}_TRG.png')
    sc.pl.rank_genes_groups_dotplot(adata_query, n_genes=ntrg,  save=f'_{bname}_TRG.png', dendrogram=False)





df=adata_query.obs.groupby(["leiden","sampleid"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat.to_csv(f'{outdir}/{bname}_conf_mat_table.csv',header=True, index=True)
conf_mat.idxmax(axis=1).to_csv(f"{outdir}/{bname}_conf_mat_table_max.csv")

row_max=conf_mat.max(axis=1)

sc.tl.dendrogram(adata_query, groupby='leiden',use_rep="X_scVI", use_raw=True, n_pcs=30 )

if mk != "none":
    lis=pd.read_csv(f'{indir}/{mk}',header=None)
    dt=lis[0].values #.set_index(0).to_dict()[1]
    sc.pl.dotplot(adata_query, dt, groupby='leiden', dendrogram=True, save=f"{bname}_marker_sanes_ct.png", use_raw=True)


dominate_name=row_max[row_max > 0.7].index
print(dominate_name)
row_max[row_max > 0.7].to_csv(f'{outdir}/{bname}_dominate_sample.csv',header=True, index=True)



#dbt_name=row_max[row_max < 0.7].index
#print(dbt_name)
#CategoricalIndex(['13', '42', '47'], categories=['0', '1', '2', '3', '4', '5', '6', '7', ...], ordered=False, dtype='category', name='leiden')
