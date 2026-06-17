import scanpy as sc
import scvelo as scv
sc.settings.figdir="HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/"
bname="Astrocyte_subclass_new5_clean_hvg_5000_raw"
dir1="HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/"
adata=sc.read("HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new5_clean_hvg_5000_raw_velocity.h5ad")
scv.pl.velocity_embedding_stream(adata, basis='umap', save=f'{dir1}/{bname}_velocity_latent_on.png', color="subclass", legend_loc="right margin")
sc.pl.violin(adata,"latent_time", groupby="subclass", save=f'{bname}_velocity_latent_time.png', rotation=90)
