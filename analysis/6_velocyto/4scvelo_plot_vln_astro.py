import scanpy as sc
import matplotlib.pyplot as plt
import scvelo as scv
plt.rcParams.update({'font.size': 14})
dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/"
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new4_clean_hvg_2000_raw_velocity.h5ad")
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/"
sc.settings.figdir = f"{outdir}"

dt={"Astro_ONH_GABBR2+":"Astro_ONH1_GABBR2+", "Astro_retina1_PAX5+GABBR2+": "Astro_ONH2_PAX5+GABBR2+", "Astro_retina2_PAX5+ME1+": "Astro_retina1_PAX5+ME1+", "Astro_retina3_NLGN1+": "Astro_retina2_NLGN1+"}
adata.obs["subclass"]=adata.obs["subclass"].replace(dt)
bname="Astrocyte_subclass_new4_clean_hvg_2000_raw"
#scv.pl.velocity_graph(adata, threshold=.1, save=f'{dir1}/{bname}_velocity_graph.png', color="subclass")

scv.pl.velocity_embedding_stream(adata, basis='umap', save=f'{dir1}/{bname}_velocity_latent.png', color="subclass")

adata.obs[["latent_time","subclass"]].groupby("subclass").mean().sort_values("latent_time").to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new4_clean_hvg_2000_raw_latent_time")

sc.pl.violin(adata,keys=["velocity_pseudotime"],groupby="subclass", save="_astro_velocity_pseudotime.png", rotation=90, order=["Astro_ON4_TSHZ2+DPP10+","Astro_ON1_TSHZ2+","Astro_ON3_AFF3+","Astro_ON2_NR4A3+","Astro_ONHON2_SLC4A11+APOE+","Astro_ONHON3_ANKUB1+","Astro_ONHON1_SLC4A11+","Astro_ONH1_GABBR2+","Astro_ONH2_PAX5+GABBR2+","Astro_retina1_PAX5+ME1+","Astro_retina2_NLGN1+"])
#sc.pl.violin(adata,keys=["latent_time"],groupby="subclass", save="_astro_latent_time.png", rotation=90)

sc.pl.violin(adata,keys=["latent_time"],groupby="subclass", save="_astro_latent_time.pdf", rotation=90,  order=["Astro_ON4_TSHZ2+DPP10+","Astro_ON1_TSHZ2+","Astro_ON3_AFF3+","Astro_ON2_NR4A3+","Astro_ONHON2_SLC4A11+APOE+","Astro_ONHON3_ANKUB1+","Astro_ONHON1_SLC4A11+","Astro_ONH1_GABBR2+","Astro_ONH2_PAX5+GABBR2+","Astro_retina1_PAX5+ME1+","Astro_retina2_NLGN1+"])

#sc.pl.violin(adata,keys=["latent_time"],groupby="subclass", save="_astro_latent_time.png", rotation=90, order=["Astro_ON4_TSHZ2+DPP10+","Astro_ON1_TSHZ2+","Astro_ON3_AFF3+","Astro_ON2_NR4A3+","Astro_ONHON2_SLC4A11+APOE+","Astro_ONHON3_ANKUB1+","Astro_ONHON1_SLC4A11+","Astro_ONH_GABBR2+","Astro_retina1_PAX5+GABBR2+","Astro_retina2_PAX5+ME1+","Astro_retina3_NLGN1+"])


#sc.pl.violin(adata[adata.obs["subclass"]!="OLIGO_NAV2-hi"],keys=["velocity_pseudotime"],groupby="subclass",order=["OPC_GPR17+","OPC","OPC_LAMA2+","OPC_RNF220+","OLIGO_LRRC7-hi","OLIGO_SVEP1-hi","OLIGO","OLIGO_RBFOX1-hi"], save="_oligo_opc_velocity_pseudotime.png")
