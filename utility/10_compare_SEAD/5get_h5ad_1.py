import scanpy as sc
import anndata as ad



#donor={'SRR29007591': 'ALSP',
# 'SRR29007592': 'ALSP',
# 'SRR29007593': 'ALSP',
# 'SRR29007594': 'ALSP',
# 'SRR29007595': 'ALSP',
# 'SRR30562939': 'AD',
# 'SRR30562940': 'Ctr',
# 'SRR30562941': 'Ctr',
# 'SRR30562942': 'Ctr',
# 'SRR30562943': 'Ctr',
# 'SRR30562944': 'Ctr',
# 'SRR30562945': 'AD'}


adata0=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/GSE267301/GSE267301_hvg_10000_epoch_none_seurat_rs_1_scvi_trg.h5ad")
adata_full0=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/GSE267301.h5ad")

adata0=adata0[adata0.obs["sampleid"].isin(['SRR30562940','SRR30562941','SRR30562942','SRR30562943','SRR30562944'])]
adata_full0=adata_full0[adata_full0.obs["sampleid"].isin(['SRR30562940','SRR30562941','SRR30562942','SRR30562943','SRR30562944'])]


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/GSE167494/GSE167494_hvg_10000_epoch_none_seurat_v3_rs_1_sb_scvi_trg.h5ad")
adata_full=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/GSE167494/GSE167494.h5ad")

cell="Astrocyte"
astro0=adata0[adata0.obs["leiden"].isin(["0","4","7","18"])].obs.index

adata_astro0=adata_full0[astro0].copy()
adata_astro0.obs["majorclass"]=cell
adata_astro0.obs["subclass"]="Astro_GSE267301"



astro=adata[adata.obs["leiden"].isin(["0","1","7","8","11","23","5","35"])].obs.index

adata_astro=adata_full[astro].copy()
adata_astro.obs["majorclass"]=cell
adata_astro.obs["subclass"]="Astro_GSE167494"


h5ad=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{cell}/clean/{cell}_subclass_new2.h5ad"

adata1=sc.read(h5ad)
adata1.X=adata1.layers["counts"]

adata_full1=ad.concat([adata_astro,adata1])


adata_full1.write(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/{cell}_PFC_ONONH.h5ad")

del adata_astro
del adata1
del adata_full1

cell="Oligodendrocyte"
astro=adata[adata.obs["leiden"].isin(["2","4","6","16","26"])].obs.index

adata_astro=adata_full[astro].copy()
adata_astro.obs["majorclass"]=cell
adata_astro.obs["subclass"]="Oligo_GSE167494"



h5ad=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{cell}/clean/{cell}_subclass_new.h5ad"

adata1=sc.read(h5ad)
adata1.X=adata1.layers["counts"]

adata_full1=ad.concat([adata_astro,adata1])


adata_full1.write(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/{cell}_PFC_ONONH.h5ad")


del adata_astro
del adata1
del adata_full1


cell="Oligodendrocyte_precursor_cell"


astro=adata[adata.obs["leiden"].isin(["10"])].obs.index


adata_astro=adata_full[astro].copy()
adata_astro.obs["majorclass"]=cell
adata_astro.obs["subclass"]="OPC_GSE167494"


h5ad=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{cell}/clean/{cell}_subclass.h5ad"

adata1=sc.read(h5ad)
adata1.X=adata1.layers["counts"]

adata_full1=ad.concat([adata_astro,adata1])


adata_full1.write(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/scvi/RNA/raw/{cell}_PFC_ONONH.h5ad")


del adata_astro
del adata1
del adata_full1

