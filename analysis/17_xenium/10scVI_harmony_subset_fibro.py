import scanpy as sc
ref="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount.h5ad"
adata=sc.read(ref)

adata1=adata[adata.obs["majorclass"]=="Fibroblast"].copy()

adata1.obs["celltype1"]=[ x.split("_")[0]+"_"+ x.split("_")[1] for x in adata1.obs["celltype"]  ]

adata1.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update/ONONH_all_rawcount_fibro.h5ad")
