import scanpy as sc
import pandas as pd
import sys
celltype=sys.argv[1]
fn=sys.argv[2]
#ct=["Oligodendrocyte","Oligodendrocyte_precursor_cell"]
#ct=["Oligodendrocyte"]
#ct=["Astrocyte"]
#ct=["Fibroblast"]
#for celltype in ct: 
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"
#    bname=f"{celltype}_subclass_new4_clean"
#    bname=f"{celltype}_subclass_seurat"
bname=f"{celltype}_{fn}"

#    bname=f"{celltype}_subclass_clean_new"
adata=sc.read_h5ad(f"{outdir}/{bname}.h5ad")
st=adata.obs["subclass"].unique()
for s in st:
    df=pd.DataFrame({"gene":adata.uns['rank_genes_groups']["names"][s], "logfoldchanges": adata.uns['rank_genes_groups']["logfoldchanges"][s], "pvals":adata.uns['rank_genes_groups']["pvals"][s], "pvals_adj":adata.uns['rank_genes_groups']["pvals_adj"][s], "scores":adata.uns['rank_genes_groups']["scores"][s]})
    out=f"{outdir}/{bname}.{s}.trg.csv"
    df.to_csv(out,index=False)

