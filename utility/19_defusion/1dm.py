import scanpy as sc
import numpy as np


outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
sc.settings.figdir = f"{outdir}/figures"

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/oligo_opc_seurat_new_scvi.h5ad")

sc.pp.neighbors(adata, use_rep = "X_scVI", n_neighbors=15, n_pcs=30)
sc.tl.diffmap(adata)
sc.pl.diffmap(adata, color='subclass',save="_oligo_opc_seurat_new_scvi_diffmap.png")

sc.pl.embedding(adata, basis='diffmap', color='subclass', components=['1,2'], save="_oligo_opc_seurat_new_scvi_diffmap_two_component.png")

# Use diffusion pseudotime (DPT)

#sc.pl.diffmap(adata, color='dpt_pseudotime', save="_oligo_opc_seurat_new_scvi_pseudotime.png")

####adata.obs['dpt_groups'] = adata.obs['subclass'] 
####sc.pl.dpt_groups_pseudotime(adata, save="_oligo_opc_seurat_new_scvi_pseudotime_groups.png")

####sc.pl.scatter(adata, x='dpt_pseudotime', y='gene_name', color='subclass')

###adata.uns["iroot"] = np.flatnonzero(adata.obs["subclass"] == "OPC_GLIS3+")


adata.uns["iroot"] = np.flatnonzero(adata.obs["subclass"] == "OPC_GLIS3+")[0]

sc.tl.dpt(adata, n_dcs=15)  # or specify root cell manually


sc.pl.umap(adata, color='dpt_pseudotime', save="_oligo_opc_seurat_new_scvi_pseudotime.png")

adata.obs['dpt_rank'] = adata.obs['dpt_pseudotime'].rank()

sc.pl.violin(adata,"dpt_pseudotime",groupby="subclass",save="_oligo_opc_seurat_new_scvi_pseudotime_violin.png")

#################

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"
sc.settings.figdir = f"{outdir}/figures"

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_seurat_new.h5ad")

sc.pp.neighbors(adata, use_rep = "X_scVI", n_neighbors=15, n_pcs=30)
sc.tl.diffmap(adata)
sc.pl.diffmap(adata, color='subclass',save="_oligo_opc_seurat_new_scvi_diffmap.png")

sc.pl.embedding(adata, basis='diffmap', color='subclass', components=['1,2'], save="_oligo_opc_seurat_new_scvi_diffmap_two_component.png")

# Use diffusion pseudotime (DPT)

#sc.pl.diffmap(adata, color='dpt_pseudotime', save="_oligo_opc_seurat_new_scvi_pseudotime.png")

####adata.obs['dpt_groups'] = adata.obs['subclass'] 
####sc.pl.dpt_groups_pseudotime(adata, save="_oligo_opc_seurat_new_scvi_pseudotime_groups.png")

####sc.pl.scatter(adata, x='dpt_pseudotime', y='gene_name', color='subclass')

###adata.uns["iroot"] = np.flatnonzero(adata.obs["subclass"] == "OPC_GLIS3+")


adata.uns["iroot"] = np.flatnonzero(adata.obs["subclass"] == "OPC_GLIS3+")[0]

sc.tl.dpt(adata, n_dcs=15)  # or specify root cell manually


sc.pl.umap(adata, color='dpt_pseudotime', save="_oligo_opc_seurat_new_scvi_pseudotime.png")

adata.obs['dpt_rank'] = adata.obs['dpt_pseudotime'].rank()

sc.pl.violin(adata,"dpt_pseudotime",groupby="subclass",save="_oligo_opc_seurat_new_scvi_pseudotime_violin.png")

