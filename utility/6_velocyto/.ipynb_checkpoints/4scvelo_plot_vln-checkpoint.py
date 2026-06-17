import scanpy as sc
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc/oligo_opc_velocity.h5ad")
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc/"
sc.settings.figdir = f"{outdir}"

sc.pl.violin(adata[adata.obs["subclass"]!="OLIGO_NAV2-hi"],keys=["velocity_pseudotime"],groupby="subclass",order=["OPC_GPR17+","OPC","OPC_LAMA2+","OPC_RNF220+","OLIGO_LRRC7-hi","OLIGO_SVEP1-hi","OLIGO","OLIGO_RBFOX1-hi"], save="_oligo_opc_velocity_pseudotime.png")
