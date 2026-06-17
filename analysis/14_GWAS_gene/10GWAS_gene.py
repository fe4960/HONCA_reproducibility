import pandas as pd
import scanpy as sc
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt

plt.rcParams.update({
    "xtick.labelsize": 15,   # column labels
    "ytick.labelsize": 15    # row labels
})

fn=["gwas-association-downloaded_2026-04-04-pubmedId_29785010.tsv",  "gwas-association-downloaded_2026-04-04-pubmedId_38242088.tsv", "gwas-association-downloaded_2026-04-04-pubmedId_27064256.tsv",  "gwas-association-downloaded_2026-04-04-pubmedId_33627673.tsv", "gwas-association-downloaded_2026-04-04-pubmedId_28073927.tsv",  "gwas-association-downloaded_2026-04-04-pubmedId_37386247.tsv"]

main="/dfs3b/ruic20_lab/junw42/reference/GWAS/glaucoma_sign_GWAS/"

df_all=pd.DataFrame()

for f in fn:
    file=f"{main}/{f}"
    df=pd.read_csv(file, header=0, index_col=0, sep="\t")
    df["MAPPED_GENE"]=df["MAPPED_GENE"].str.split(r",\s+|\s+-\s+")
#    df["MAPPED_GENE"]=df["MAPPED_GENE"].str.split("\s+-\s+")
    df=df.explode("MAPPED_GENE")
    df = df[df["MAPPED_GENE"].notna()]
    df["MAPPED_GENE"] = df["MAPPED_GENE"].str.strip()
    df_all=pd.concat([df_all, df], axis=0)

df_all=df_all.reset_index(drop=True)
df_all.to_csv(f"{main}/combined_GWAS_association", sep="\t")

df_all.loc[df_all["DISEASE/TRAIT"].str.contains("open"),"MAPPED_GENE"].value_counts()

#poag=df_all.loc[df_all["DISEASE/TRAIT"].str.contains("open"),"MAPPED_GENE"].drop_duplicates()

##########poag=df_all.loc[df_all["DISEASE/TRAIT"].str.contains("open"),"MAPPED_GENE"].value_counts().reset_index()

poag=df_all["MAPPED_GENE"].value_counts().reset_index()


poag=poag.loc[poag["count"]>=2]
#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_raw_normcount_only_2k_celltype_rmUnk_simple.h5ad")


dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final"

dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/"

path=f"{dir}/ONONH_all_raw_normcount_only.h5ad"
adata=sc.read_h5ad(path)

adata=adata[(~adata.obs["majorclass1"].isin(["Schwann_cell", "Pigmented_cell"]))&(adata.obs["celltype"]!="Unknown?")]

#poag_gene=adata.var_names.interaction(poag)

#poag_gene
sc.settings.figdir=f"{main}/figures"
#sc.pl.violin(adata,poag_gene, groupby="celltype", save="_GWAS_gene_celltype.png")
#sc.pl.violin(adata,poag_gene, groupby="majorclass", save="_GWAS_gene_majorclass.png")


df=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/gene_exp_celltype_cpm.txt.gz", sep="\t", header=0, index_col=0)

df.drop(columns=["Pigmented_cell_1", "Pigmented_cell_2", "Schwann_cell_myelinating", "Schwann_cell_nonmyelinating"], inplace=True)

df_filtered = df[df.max(axis=1) >= 500]


poag_gene=adata.var_names[(adata.var_names.isin(poag["MAPPED_GENE"]))&(adata.var_names.isin(df_filtered.index))]

##########X = df_filtered.loc[poag_gene, :].values 

###########Z = sch.linkage(X, method='ward')
##########gene_order = sch.leaves_list(Z)

############genes_ordered = df_filtered.index[gene_order]




sc.tl.dendrogram(adata, groupby="celltype")
sc.pl.dotplot(adata, poag_gene, groupby="celltype", save="_GWAS_gene_celltype.png", dendrogram=True)
#sc.pl.dotplot(adata, genes_ordered, groupby="celltype", save="_GWAS_gene_celltype_cpm500_ct2.pdf", swap_axes=True, dendrogram=True)

df=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/gene_exp_majorclass_cpm.txt.gz", sep="\t", header=0, index_col=0)

df.drop(columns=[ "Pigmented_cell", "Schwann_cell"], inplace=True)


df_filtered = df[df.max(axis=1) >= 500]

#poag_gene=adata.var.loc[(adata.var_names.isin(poag["MAPPED_GENE"]))&(adata.var_names.isin(df_filtered.index)),].index

poag_gene=adata.var_names[(adata.var_names.isin(poag["MAPPED_GENE"]))&(adata.var_names.isin(df_filtered.index))]


#############X = df_filtered.loc[poag_gene, :].values

############Z = sch.linkage(X, method='ward')
#############gene_order = sch.leaves_list(Z)

############genes_ordered = df_filtered.index[gene_order]

#sc.pl.dotplot(adata,poag_gene, groupby="majorclass", save="_GWAS_gene_majorclass.png")
sc.tl.dendrogram(adata, groupby="majorclass1")
sc.pl.dotplot(adata, poag_gene , groupby="majorclass1", save="_GWAS_gene_majorclass_cpm500_ct2.pdf", swap_axes=False, dendrogram=True)



