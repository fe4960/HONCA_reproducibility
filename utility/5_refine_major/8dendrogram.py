import scanpy as sc
import pandas as pd

celltype="Astrocyte"

query=f"HCA_ON/data/5_refine_major/scvi/{celltype}/clean/{celltype}_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad"

indir=f"HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"HCA_ON/data/5_refine_major/scvi/{celltype}/clean/"

sc.settings.figdir = f"{outdir}/figures"
seed=7
res=0.7
rmlist="none" #"rm_cluster_rs_0.3"

fl="seurat"
hvg="10000"
max_epoch="none"
gl="zinb"
span="1"
kp="none"
batch_key="sampleid"
obs=pd.read_csv(f"HCA_ON/data/5_refine_major/scvi/{celltype}/clean/{celltype}_hvg{hvg}_epoch{max_epoch}_{fl}_rs_{res}_clean_sb_{kp}.obs.gz",index_col=0,header=0)
bname=f"{celltype}_hvg{hvg}_epoch{max_epoch}_{fl}_rs_{res}_clean_sb_{kp}"
adata=sc.read_h5ad(query)
adata.obs["leiden"]=obs.loc[adata.obs.index,"leiden"].astype("str").astype("category")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sampleid")
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.dendrogram(adata, groupby="leiden", n_pcs=30, use_rep="X_pca")
mk="sanes_mk"
lis=pd.read_csv(f'{indir}/{mk}',header=None)
dt=lis[0].values #.set_index(0).to_dict()[1]
adata.write(f"{outdir}/{bname}_dendrogram.h5ad")
sc.pl.dotplot(adata, dt, groupby='leiden', dendrogram=True, save=f"{bname}_marker_sanes_ct.png", use_raw=False)


