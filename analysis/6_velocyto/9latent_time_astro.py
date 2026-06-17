import scanpy as sc

from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

rcParams["figure.figsize"] = (5,5)
rcParams.update({'font.size': 14})

#sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc/"
#########adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc/oligo_opc_seurat_new_scvi_hvg_2000_raw_velocity.h5ad")
##########desired_order = ["OLIGO2_LRRC7+", "OLIGO1", "OLIGO1_SVEP1+", "OLIGO1_RBFOX1+"]
########sc.pl.violin(adata[adata.obs["majorclass"]=="Oligodendrocyte"],"latent_time",groupby="subclass", order=desired_order, save="_opc_latent_time.pdf",rotation=90, jitter=False)


##########sc.pl.violin(adata[adata.obs["majorclass"]=="Oligodendrocyte"],"latent_time",groupby="subclass", order=desired_order, save="_opc_latent_time.png")
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/"

#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/veloc/Oligodendrocyte_precursor_cell_subclass_seurat_hvg_2000_raw_velocity_reservse.h5ad")
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/veloc/Astrocyte_subclass_new5_clean_hvg_5000_raw_velocity.h5ad")
#########adata=adata[(adata.obs["sampleid"]!="MMD_23_17976_ONH_RNA")&(adata.obs["sampleid"]!="MMD_23_20181_ONH_RNA")]
#dt={"COP_LAMA2+":"NFOL_UTRN+"}
#adata.obs["subclass"]=adata.obs["subclass"].replace(dt)
#order=["OPC_GLIS3+","OPC","OPC_Cycling", "COP_NFOL", "NFOL_UTRN+", "MFOL"]
#order=["Astro_ON5_DPP10+", "Astro_ON4"]
dt={"Astro_ON4_AFF3+DMGDH+": "Astro_ON4", "Astro_ON5_DPP10+": "Astro_ON5", "Astro_ONHON1_SLC4A11+MARCH1+": "Astro_ONHON1", "Astro_ON1_WNK2+TSHZ2+": "Astro_ON1", "Astro_ON2_ACTN1+SERPINA3+": "Astro_ON2", "Astro_ON3_NR4A3+FOS+": "Astro_ON3", "Astro_ONHON2_CST3+APOE+": "Astro_ONHON2", "Astro_ONH2_PAX5+GABBR2+": "Astro_ONH2", "Astro_ONH1_GABBR2+SV2B+": "Astro_ONH1", "Astro_retina1_PAX5+ME1+": "Astro_retina1", "Astro_retina2_NLGN1+": "Astro_retina2"}
adata.obs["subclass1"]=adata.obs["subclass"].astype("str").replace(dt)
order=["Astro_ON5", "Astro_ON4", "Astro_ON3", "Astro_ON2", "Astro_ON1", "Astro_ONHON1", "Astro_ONHON2", "Astro_ONH1", "Astro_ONH2", "Astro_retina1", "Astro_retina2" ]
#sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/veloc/"
sc.pl.violin(adata, "latent_time",groupby="subclass1", order=order, save="_astro_latent_time_reorder.pdf",rotation=90, jitter=False)
