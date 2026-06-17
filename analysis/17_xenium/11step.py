from step import crossModel
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import torch

adata_list=[]

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/"

for i in range(1,2):
    h5ad=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_merge_cutoff20_major_{i}.h5ad"
    adata=sc.read(h5ad)
    csv=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/xenium_spatial_major_{i}.csv"
    spatial=pd.read_csv(csv,index_col=0,header=0)
    adata.obsm["spatial"]=np.array(spatial.loc[adata.obs.index,], dtype=float)
    adata.obs["array_row"]=spatial.loc[adata.obs.index,"x"]
    adata.obs["array_col"]=spatial.loc[adata.obs.index,"y"]
    txt=f"HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_combcoembed.combined_inte_{i}.txt"
    anno=pd.read_csv(txt, index_col=1,header=0, sep="\t")
    adata.obs["harmony_anno"]=anno.loc[adata.obs.index,"harmony_anno"]
    adata.obs["library_key"]=f"slide_{i}"
    adata.obs["sampleid"]=f"slide_{i}"
    adata.obs.index=adata.obs.index.astype(str)+"_"+str(i)
    adata_list.append(adata)

st_adata=ad.concat(adata_list)

sc_adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/ONONH_all_rawcount_subset.h5ad")
st_pth="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/step"

common_gene=st_adata.var.index.intersection(sc_adata.var.index)

sc_adata=sc_adata[:, common_gene].copy()
sc_adata.obs["modality"]="SC"
sc_adata.var_names=sc_adata.var.index

st_adata=st_adata[:, common_gene].copy()
st_adata.obs["modality"]="ST"
st_adata.var_names=st_adata.var.index

#common_gene=common_gene.astype("str")

stepc = crossModel(
    sc_adata=sc_adata,
    st_adata=st_adata,
    n_top_genes=None,
    batch_key='sampleid',
    class_key='majorclass',
    st_batch_key='library_key',
    module_dim=30,
    hidden_dim=64,
    n_modules=32,
    decoder_type='zinb',
    dec_norm=None,
    n_dec_hid_layers=1,
)

#stepc.integrate(
#    epochs=400,
#    batch_size=1024,
#    split_rate=0.2,
#    need_anchors=True  #False, # whether use cell-type annotations to postprocess the co-embedding
#)


stepc.integrate(epochs=400,
               batch_size=2048,
               split_rate=0.2,
               tune_epochs=200,
               need_anchors=False,
               beta1=1e-2)


# Save the model config and weights
#########stepc.save(f'{st_pth}/crossmodel_config')

stepc.impute(layer_key='counts', qc=True)

adata = stepc.adata

sc.pp.neighbors(adata, use_rep='X_rep')
sc.tl.umap(adata)
sc.pl.umap(adata, color=['modality', 'majorclass'], ncols=2, save="_ONONH_sc_st_coembed_step.png")
adata.write(f'{st_pth}/crossmodel.h5ad')
torch.save(stepc._functional.model.state_dict(), st_pth)


