import scanpy as sc
main="/dfs3b/ruic20_lab/junw42"
#adata=sc.read(f"{main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg_2000_fl_seurat_v3_res_0.4_sd_7_ononh_retina_nonsb_scvi_trg.h5ad")
#adata=adata[adata.obs["subclass"]!="Astro_retina_HRCA"]
#adata=adata[adata.obs["leiden"]!=6]
#dt={"1":"Astro_ONH","3":"Astro_retina", "0":"Astro_ON"}
#adata.obs["subclass_original"]=adata.obs["subclass"].copy()
#adata.obs["subclass"]=adata.obs["leiden"].replace(dt)
#adata.obs.loc[adata.obs["leiden"].isin(["0","2","4","5","6","7","8","9"]),"subclass"]="Astro_ON"
#adata.obs["subclass"]=adata.obs["subclass"].cat.remove_unused_categories()
#adata.write(f"{main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg_2000_fl_seurat_v3_res_0.4_sd_7_ononh_retina_nonsb_subclass.h5ad")

adata=sc.read(f"{main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg_2000_fl_seurat_v3_res_0.4_sd_7_ononh_retina_nonsb_subclass.h5ad")

adata1=sc.read(f"{main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new.h5ad")

adata1=adata1[adata.obs.index]
adata1.obs["subclass"]=adata.obs.loc[adata1.obs.index,"subclass"]
adata1.write(f"{main}/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new_update.h5ad")
