import scanpy as sc

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/figures/"

h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.h5ad"

data=sc.read(h5ad)

sc.pl.umap(data,color="majorclass", save="_ONONH_all_rawcount_subset.png")
