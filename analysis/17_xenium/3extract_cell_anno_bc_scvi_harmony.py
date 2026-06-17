import scanpy as sc
import matplotlib.pyplot as plt
#import sys
import pandas as pd

#cell=sys.argv[1]
#label=sys.argv[2]

#indir="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_crossmap_epi"
#adata=sc.read_h5ad(f'{indir}/{cell}_combined-emb.h5ad')

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi/figures/"

label="sn_st_scvi_harmony"


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi/" #sys.argv[1]

adata=sc.read(f"{indir}/{label}_coembed_flt_dbt.h5ad")

cell=["Oligodendrocyte", "Astrocyte", "Fibroblast"]
#sample=adata.obs["sampleid"].unique()
sample=['ONH_1', 'ONH_2', 'ONH_3', 'PP_1', 'PP_2']

for c in cell:
    for s in sample:
        data=adata.obs.loc[(adata.obs["sampleid"]==s)&(adata.obs["lr_celltype"]==c), "barcode"]
        data.to_csv(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_{s}_{c}")


