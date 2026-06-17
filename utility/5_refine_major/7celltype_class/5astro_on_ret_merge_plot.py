import scanpy as sc
from matplotlib import rcParams
import pandas as pd
import matplotlib.pyplot as plt
rcParams["figure.figsize"] = (5,5)
dir1="/dfs3b/ruic20_lab/junw42"
sc.settings.figdir = f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/figures"
obs=pd.read_csv(f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.3_clean_sb_rm0_12_15_16_none_rm4_9_rm2_seed_7.obs.gz",header=0,index_col=0)
adata=sc.read(f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7_scvi_trg.h5ad")
adata=adata[obs.index]
adata.obs["leiden1"]=obs.loc[adata.obs.index,"leiden"].astype("str")
sc.pl.umap(adata,color="leiden1",save="Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7_leiden1_onData.png", legend_loc="on data")
bname="Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7"
outdir=f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
df=adata.obs[["nCount_RNA","tissue","leiden1"]].groupby(["leiden1","tissue"]).count().to_csv(f"{outdir}/{bname}_tissue.csv")
df=pd.read_csv(f"{outdir}/{bname}_tissue.csv")
total=df.groupby(["leiden1"])["nCount_RNA"].transform("sum")
df["nCount_RNA"]=df["nCount_RNA"]/total

plt.barh( df[df["tissue"]!="ON"]["leiden1"], df[df["tissue"]!="ON"]["nCount_RNA"], color='blue' )
plt.barh(df[df["tissue"]=="ON"]["leiden1"], df[df["tissue"]=="ON"]["nCount_RNA"],  left=df[df["tissue"]!="ON"]["nCount_RNA"], color='lightblue' )
plt.legend(["ONH", "ON"])
plt.xlabel("Cell proportion")
plt.savefig(f"{outdir}/figures/{bname}_tissue.png")



sc.tl.rank_genes_groups(adata, "leiden1", method="wilcoxon", use_raw=False)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=10,  save=f'{bname}_TRG10_leiden1.png', dendrogram=False)
adata.write(f"{dir1}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg10000_epochnone_seurat_rs_0.6_clean_sb_rm0_12_15_16_merge_seed_7_scvi_trg_update.h5ad")

celltype="Astrocyte"

