import scanpy as sc
main="/dfs3b/ruic20_lab/junw42"
adata=sc.read(f"{main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_seed_7_res_0.4_clean_ononh_retina_scvi_trg.h5ad")
adata=adata[adata.obs["leiden"]!="10"]
adata.write(f"{main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_seed_7_res_0.4_clean_ononh_retina_scvi_trg_rm_clu10.h5ad")
