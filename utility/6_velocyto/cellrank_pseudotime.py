import numpy as np
import cellrank as cr
import scanpy as sc
import scvelo as scv
import sys
scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2

dir="/dfs3b/ruic20_lab/junw42"

dir1=f"{dir}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/cellrank"

import os
if not os.path.exists(dir1):
    os.makedirs(dir1)

#adata=sc.read(f'{dir1}/{bname}_seurat_new.h5ad')

adata=sc.read(f"{dir}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat.h5ad")

sc.pp.neighbors(adata, use_rep = "X_scVI", n_neighbors=15, n_pcs=30)
sc.tl.diffmap(adata)

root_ixs =adata.obsm['X_diffmap'][:, 3].argmax()

scv.pl.scatter(
    adata,
    basis="diffmap",
    c=["subclass", root_ixs],
    legend_loc="right",
    components=["2, 3"],
    save=f"{dir1}/Oligodendrocyte_precursor_cell_diffmap_c2c3.png"
)


scv.pl.scatter(
    adata,
    basis="diffmap",
    c=["subclass", root_ixs],
    legend_loc="right",
    components=["1,2"],
    save=f"{dir1}/Oligodendrocyte_precursor_cell_diffmap_c1c2.png"
)

adata.uns["iroot"] = np.flatnonzero(adata.obs["subclass"] == "OPC")[0]

sc.tl.dpt(adata, n_dcs=15)

sc.settings.figdir=dir1

sc.pl.umap(adata, color='dpt_pseudotime', save=f"_opc_subclass_seurat_pseudotime.png")

adata.obs['dpt_rank'] = adata.obs['dpt_pseudotime'].rank()

sc.pl.violin(adata,"dpt_pseudotime",groupby="subclass",save="_opc_subclass_seurat_pseudotime_violin.png")

pk = cr.kernels.PseudotimeKernel(adata, time_key="dpt_pseudotime")
pk.compute_transition_matrix()

print(pk)

pk.plot_projection(basis="umap", recompute=True,save=f"{dir1}/opc_subclass_seurat_pseudotime_proj_umap.png")
