import scanpy as sc
from matplotlib import rcParams
import pandas as pd
rcParams["figure.figsize"] = (5,5)
sc.settings.figdir = "HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/figures"

#adata=sc.read("HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_seed_7_res_0.4_clean_ononh_retina_scvi_trg.h5ad")
#sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,palette="tab20", save=f"_Astrocyte_seed_7_res_0.4_clean_ononh_retina_subclass.png")


#####adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg_2000_fl_0.4_res_0.4_ononh_retina_mg_sb_scvi_trg.h5ad")
#####sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,palette="tab20", save=f"_Astrocyte_hvg_2000_fl_0.4_res_0.4_ononh_retina_mg_sb_subclass.png")

#####df=adata.obs.groupby(["subclass","leiden"]).size().unstack(fill_value=0)
#####conf_mat = df / df.sum(axis=1).values[:,np.newaxis]
#####conf_mat.to_csv(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_hvg_2000_fl_0.4_res_0.4_ononh_retina_mg_sb_conf.csv',header=True, index=True)

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/major_hvg_2000_fl_0.4_res_0.4_ononh_retina_mg_sb_scvi_trg.h5ad")
sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,palette="tab20", save=f"_major_hvg_2000_fl_0.4_res_0.4_ononh_retina_mg_sb_subclass.png")

df=adata.obs.groupby(["subclass","leiden"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=1).values[:,np.newaxis]
conf_mat.to_csv(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/major_hvg_2000_fl_0.4_res_0.4_ononh_retina_mg_sb_conf.csv',header=True, index=True)
