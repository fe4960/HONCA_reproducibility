import scanpy as sc
import anndata as ad
adata=sc.read("HCA_ON/data/5_refine_major/scvi/major/major_ON_ONH_new_update_final.h5ad")
adata=adata[~adata.obs["majorclass1"].isin(["Pigmented_cell", "Schwann_cell", "Adipocyte", "RPE"])]
adata=adata[adata.obs["tissue"]=="ONH"]
sc.settings.figdir="HCA_ON/data/5_refine_major/scvi/major/figures/"
#sc.pl.dotplot(adata,"RNASEH2B", groupby="majorclass1", save="_RNASEH2B_exp_ONONH.png", swap_axes=True)

adata.obs["majorclass"]=adata.obs["majorclass1"]

adata=adata[~adata.obs["majorclass"].isin(["Rod", "Cone", "BC", "HC", "AC", "RGC", "MG"])]

adata1=sc.read("/dfs3b/ruic20_lab/junw42/human_ret_anc/data/snRNA/retina_snRNA_majorclass_5k.h5ad")
#adata1=sc.read("/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HRCA/CAP/fix_donor/scrnah5adobssubset2replacebydict/HRCA_snRNA_allcells_normalized.h5ad")
adata1=adata1[~adata1.obs["majorclass"].isin(["RPE", "Astrocyte", "Microglia"])]
#sc.pl.dotplot(adata,"RNASEH2B", groupby="majorclass", save="_RNASEH2B_exp_retina.png", swap_axes=True)


adata2=ad.concat([adata,adata1],  join="inner")

sc.pl.dotplot(adata2,"RNASEH2B", groupby="majorclass", save="_RNASEH2B_exp_retina_5k_ONH.png", swap_axes=True)


#sc.pl.dotplot(adata,"RNASEH2B", groupby="majorclass", save="_RNASEH2B_exp_retina_5k.png", swap_axes=True)

