import scanpy as sc
from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)


#adata=asc.read("/dfs3b/ruic20_lab/jinjingj/human_atlas/sn/sn_TMCB_v2_36601_normalize_clean.h5ad")
adata=sc.read("HCA_ON/data/5_refine_major/scvi/major/clean/major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_rmUnk_clean.h5ad")
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/figures"
mk="ARHGEF3"
sc.pl.violin(adata,mk, groupby="majorclass", save="_{mk}_majorclass_nodot.pdf", jitter=False, rotation=90)
sc.pl.dotplot(adata, mk, groupby="majorclass", save=f"_{mk}_majorclass.pdf", swap_axes=True)

adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean.h5ad')

level=[  "Astro_retina1_PAX5+ME1+", "Astro_retina2_NLGN1+", "Astro_ONH1_GABBR2+SV2B+", "Astro_ONH2_PAX5+GABBR2+", "Astro_ONHON1_SLC4A11+MARCH1+", "Astro_ONHON2_CST3+APOE+",  "Astro_ON1_WNK2+TSHZ2+","Astro_ON2_ACTN1+SERPINA3+", "Astro_ON3_NR4A3+FOS+", "Astro_ON4_AFF3+DMGDH+", "Astro_ON5_DPP10+"]


adata.obs["subclass"] = (
    adata.obs["subclass"]
    .astype("category")
    .cat.set_categories(level, ordered=True)
)
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/figures"

sc.pl.dotplot(adata, mk, groupby="subclass", save=f"_{mk}_astro_sub_class.pdf", swap_axes=True)
sc.pl.violin(adata,mk, groupby="subclass", save="_{mk}_astro_sub_class_nodot.pdf", jitter=False, rotation=90)


#asc.pl.violin(adata,mk, groupby="majortype", save="_PLEKHA7_majorclass_nodot.pdf", jitter=False, rotation=90)
#asc.pl.violin(adata,mk, groupby="majortype", save="_PLEKHA7_majorclass.png")
