import scanpy as sc

main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/"
path=f"{main}/5_refine_major/scvi/lattice/ONONH_Astrocyte_rawcount.h5ad"
adata=sc.read_h5ad(path)

adata.obs["subclass"]="Astro_ON"

adata.obs.loc[adata.obs["celltype"].isin(["Astro_retina","Astro_retina_NLGN1+","Astro_ONH"]),"subclass"]="Astro_ret"

adata.write(f"{main}/5_refine_major/scvi/Astrocyte/clean/ONONH_Astrocyte_rawcount_2subclass.h5ad")
