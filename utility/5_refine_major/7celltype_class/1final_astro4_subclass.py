import scanpy as sc

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple.h5ad")

adata.obs["celltype"].value_counts()

adata.obs["celltype"]=adata.obs["subclass"]

dt={"Astro_ONHON1_SLC4A11+MARCH1+": "Astro_ONHON",
"Astro_ON1_WNK2+TSHZ2+": "Astro_ON",
"Astro_ON2_ACTN1+SERPINA3+": "Astro_ON",
"Astro_ON3_NR4A3+FOS+": "Astro_ON",
"Astro_ONH2_PAX5+GABBR2+": "Astro_ONH2",
"Astro_ONH1_GABBR2+SV2B+": "Astro_ONH1",
"Astro_ONHON2_CST3+APOE+":   "Astro_ONHON",
"Astro_retina1_PAX5+ME1+":         "Astro_retina",
"Astro_retina2_NLGN1+":           "Astro_retina",
"Astro_ON4_AFF3+DMGDH+":            "Astro_ON",
"Astro_ON5_DPP10+": "Astro_ON"         }

adata.obs["subclass"]=adata.obs["celltype"].replace(dt)

adata.obs["subclass"].value_counts()

adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_celltype.h5ad")

dt={"Astro_ONH2": "Astro_ONH", "Astro_ONH1": "Astro_ONH" }

adata.obs["subclass1"]=adata.obs["subclass"].replace(dt)

dt={"Astro_ONHON": "Astro_ON", "Astro_ONH": "Astro_ONHret", "Astro_retina": "Astro_ONHret" }

adata.obs["subclass2"]=adata.obs["subclass1"].replace(dt)

adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_celltype.h5ad")
