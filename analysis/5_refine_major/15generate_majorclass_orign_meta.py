import scanpy as sc
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/major_hvg2000_nonsb_epochnone_seurat_v3_rs_1_clean_ON_ONH_new_scvi_trg.h5ad")
df=adata.obs[['reactionid', 'sampleid', 'donor', 'race', 'gender', 'age', 'tissue','sampleid_legacy','age_year']]
df1=df.drop_duplicates()
df1.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_origin",index=False)
