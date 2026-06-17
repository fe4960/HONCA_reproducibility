import scanpy as sc
import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd
sys.setrecursionlimit(10000)

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

h5ad=sys.argv[1]
res=float(sys.argv[2])
batch_key=sys.argv[3]
bname=sys.argv[4]
outdir=sys.argv[5]
rm=sys.argv[6]
marker=sys.argv[7]
save=sys.argv[8]
indir=sys.argv[9]
adata=sc.read(h5ad)


if rm != "none":
    lis=pd.read_csv(f'{outdir}/{rm}',header=None)
    rmcluster=lis[0].values.astype(str)
    adata=adata[-adata.obs["leiden"].isin(rmcluster)]


sc.settings.figdir = f"{outdir}/figures"

sc.tl.leiden(adata, resolution=res)
sc.tl.umap(adata)

#plot=[batch_key, "tissue","race", "gender", "source", "age_year"]

#for pl in plot:
#    sc.pl.embedding(adata, basis="X_umap", color=[pl],ncols=1,frameon=False,save=f'_{bname}_{pl}.png', palette="tab20")

sc.pl.embedding(adata, basis="X_umap", color=["leiden"],legend_loc="on data", ncols=1,frameon=False,save=f'_{bname}_leiden_onData.png', palette="tab20")

sc.pl.violin(adata,keys=["pct_counts_mt","n_genes_by_counts", "total_counts"], groupby="leiden", rotation=75, save=f"_{bname}_QC.png", multi_panel=True)

markers={"B cell": ["BANK1", "BLK", "BCL11A"], "Schwann cell": ["SCN7A","INSC","L1CAM","NGFR"], "T cell": ["BCL11B", "THEMIS", "IL7R", "CAMK4", "CD2"], "Astrocyte": ["AC092957.1", "DCLK1", "PAX8", "GFAP", "ALDH1L1"], "Endothelial cell": ["PTPRB", "LDB2", "VWF", "FLT1", "ANO2", "MECOM","PECAM1"], "Fibroblast": ["CLMP","BICC1","SCARA5","FHL2","DCN"], "Macrophage": ["CD163", "F13A1","MRC1", "STAB1"], "Melanocyte": ["ABCB5","PAX3","SYNPR","PMEL","GALNTL6"], "Microglia": ["MPZ","MLIP","GALNT17","PRX","CD74","CX3CR1"], "NK cell": ["KLRD1","CD247","TXK","KLRF1"], "Oligodendrocyte": ["ST18","CERCAM","TMEM144","MOBP","SH3GL3","MOG"], "Oligodendrocyte precursor cell": ["CSMD1","TNR","OPCML","SNTG1","CSMD3","OLIG2","OLIG1"], "Mular pericyte": ["RGS5","ABCC9","TRPC4"], "RPE": ["CCBE1","SLC38A11","ATP6V1C2","OTX2","BEST1","RPE65"], "Mular vascular associated smooth muscle cell": ["RGS6","MYH11","ACTA2","UNC13C"] , "Rod": "PDE6A", "Cone": "ARR3", "BC": "VSX2", "AC": "PAX6", "RGC": "RBPMS", "HC": "ONECUT1", "MG": "RLBP1"}

markers1 = dict(sorted(markers.items()))

known_markers={"B cell": ["DOCK2", "PTPRC", "EBF1","BCL11A"], "T cell": ["DOCK2","PTPRC","CD2","THEMIS"], "Astrocyte": ["ALDH1L1", "GFAP"], "Endothelial cell": ["FLT1","PECAM1"], "Fibroblast": ["COL1A2","DCN"], "Macrophage": ["DOCK2","PTPRC","MRC1","CD163"], "Microglia": ["DOCK2", "PTPRC","P2RY12", "ITGAM", "CD74", "CX3CR1"], "NK cell": ["DOCK2", "PTPRC", "KLRD1","KLRF1"], "Oligodendrocyte": ["MOG","MOBP"], "Oligodendrocyte precursor cell": ["OLIG1","OLIG2"], "Mular pericyte": ["TAGLN","ACTA2","TRPC4","RGS5","LAMA2"], "Mular vacular associated smooth muscle cell": ["TAGLN","ACTA2","UNC13C","RGS6","CASQ2","FGF7","SLC7A2","CLSTN2"],"Schwann cell": ["SCN7A","NCAM1","L1CAM","NGFR"], "Schwann cell myelinating": ["MPZ","PRX","MBP"], "RPE": ["MLANA","TYR","BEST1","RPE65"], "Melanocyte": ["PMEL","PAX3"], "Rod": "PDE6A", "Cone": "ARR3", "BC": "VSX2", "AC": "PAX6", "RGC": "RBPMS", "HC": "ONECUT1", "MG": "RLBP1" }

known_markers1 = dict(sorted(known_markers.items()))

#sc.pl.dotplot(adata, markers1, groupby='leiden', dendrogram=False, save=f"{bname}_marker.png", use_raw=True)
#sc.pl.dotplot(adata, known_markers1, groupby='leiden', dendrogram=False, save=f"{bname}_known_marker.png", use_raw=True)

sc.pl.dotplot(adata, markers1, groupby='leiden', dendrogram=False, save=f"{bname}_marker.png", use_raw=True)
sc.pl.dotplot(adata, known_markers1, groupby='leiden', dendrogram=False, save=f"{bname}_known_marker.png", use_raw=True)

sc.tl.dendrogram(adata, groupby='leiden',use_rep="X_scVI", use_raw=True, n_pcs=30 )

#sc.tl.dendrogram(adata, groupby='leiden',use_rep="X_scVI", use_raw=True, n_pcs=30 )
sc.pl.dendrogram(adata, groupby='leiden', save=f"{bname}_dendrogram.png")



#sc.pl.dotplot(adata, markers, groupby='bulk_labels', dendrogram=True)

if marker != "none":
    lis=pd.read_csv(f'{indir}/{marker}',header=None)
    dt=lis[0].values #.set_index(0).to_dict()[1]
    sc.pl.dotplot(adata, dt, groupby='leiden', dendrogram=True, save=f"{bname}_{marker}.png", use_raw=True)    
#    sc.pl.dotplot(adata, dt, groupby='leiden', dendrogram=True, save=f"{bname}_{marker}.png", use_raw=True)    


adata.obs[["nCount_RNA","leiden","sampleid"]].groupby(["leiden","sampleid"]).count().to_csv(f"{outdir}/{bname}_sample_number.txt",sep="\t")

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)

#sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)

sc.pl.rank_genes_groups(adata, n_genes=20,sharey=False, save=f'_{bname}_TRG20.png')

sc.pl.rank_genes_groups_dotplot(adata, n_genes=10,  save=f'{bname}_TRG10.png', dendrogram=False)


top_genes = pd.DataFrame({
    group: adata.uns['rank_genes_groups']['names'][group][:20]
    for group in adata.uns['rank_genes_groups']['names'].dtype.names
})
print(top_genes)

top_genes.to_csv(f"{outdir}/{bname}_top20.txt",sep="\t")


df=adata.obs.groupby(["leiden","sampleid"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:, np.newaxis]
conf_mat.to_csv(f'{outdir}/{bname}_sample_conf_mat_table.csv',header=True, index=True)
row_max=conf_mat.max(axis=1)
if save == "t":
    adata.write(f'{outdir}/{bname}.h5ad')
dominate_name=row_max[row_max > 0.7].index
print(dominate_name)
adata.obs.to_csv(f'{outdir}/{bname}.obs.gz')

#to_csv(f'{outdir}/{bname}_dominate_sample.csv',header=True, index=True)
