import scanpy as sc
#adata
adata_query=sc.read_h5ad("/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/scvi/_ONONH_hvg2k_epoch20_scvi.h5ad")
#adata_query=sc.read_h5ad("/dfs3b/ruic20_lab/junw42/ONONH_test.h5ad")
#adata_query.raw = adata_query
#del adata_query.var['gene_symbol']
sc.tl.rank_genes_groups(adata_query, "leiden", method="wilcoxon",use_raw=False)
#adata_query.write("/dfs3b/ruic20_lab/junw42/ONONH_test_trg.h5ad")
sc.pl.rank_genes_groups(adata_query, n_genes=25  , save=f'_ONONH_hvg2k_epoch20_scvi_TRG.png', sharey=False)


adata_query.var["mt"] = adata_query.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata_query, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)


sc.pl.violin(data,keys=["pct_counts_mt","n_genes_by_counts", "total_counts"], groupby="leiden", rotation=75, save="_ONONH_hvg2k_epoch20_scvi_vln_QC.png", multi_panel=True)

