import scanpy as sc
import anndata as ad
import pandas as pd

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/SEAD"

df=pd.read_csv("/dfs3b/ruic20_lab/junw42/reference/cellranger/refdata-gex-GRCh38-2020-A_ensemblID_genename",header=None,index_col=0,sep="\t")
dt=df.to_dict()[1]

Astro1=sc.read(f"{dir}/Astrocyte_DLPFC.h5ad")
Astro1=Astro1[Astro1.obs["assay"]=="10x 3' v3"]
Astro1.obs["majorclass"]="Astro_DLPFC"
Astro1.obs["celltype"]="DLPFC_"+Astro1.obs["Supertype"].astype(str)
Astro1.obs["sampleid"]="DLPFC_"+Astro1.obs["Specimen ID"].astype(str)

ensembl_names=pd.DataFrame(Astro1.var_names.values)
Astro1.var.index=ensembl_names[0].replace(dt)
Astro1.var_names=ensembl_names[0].replace(dt)
Astro1.var_names_make_unique()

Astro2=sc.read(f"{dir}/Astrocyte_MTG.h5ad")
Astro2=Astro2[Astro2.obs["assay"]=="10x 3' v3"]
Astro2.obs["majorclass"]="Astro_MTG"
Astro2.obs["celltype"]="MTG_"+Astro2.obs["Supertype"].astype(str)
Astro2.obs["sampleid"]="MTG_"+Astro2.obs["Specimen ID"].astype(str)


ensembl_names=pd.DataFrame(Astro2.var_names.values)
Astro2.var.index=ensembl_names[0].replace(dt)
Astro2.var_names=ensembl_names[0].replace(dt)
Astro2.var_names_make_unique()



Astro=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new2.h5ad')
Astro.obs["majorclass"]="Astro_ONONH"
Astro.obs["celltype"]="ONONH_"+Astro.obs["subclass"].astype(str)


Astro_all=ad.concat([Astro1,Astro2,Astro])

Astro_all.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/merged/Astro.h5ad", compression="gzip")

del Astro_all

Microglia1=sc.read(f"{dir}/Microglia_PVM_DLPFC.h5ad")
Microglia1=Microglia1[Microglia1.obs["assay"]=="10x 3' v3"]
Microglia1.obs["majorclass"]="Micro_DLPFC"
Microglia1.obs["celltype"]="DLPFC_"+Microglia1.obs["Supertype"].astype(str)
Microglia1.obs["sampleid"]="DLPFC_"+Microglia1.obs["Specimen ID"].astype(str)


ensembl_names=pd.DataFrame(Microglia1.var_names.values)
Microglia1.var.index=ensembl_names[0].replace(dt)
Microglia1.var_names=ensembl_names[0].replace(dt)
Microglia1.var_names_make_unique()


Microglia2=sc.read(f"{dir}/Microglia_PVM_MTG.h5ad")
Microglia2=Microglia2[Microglia2.obs["assay"]=="10x 3' v3"]
Microglia2.obs["majorclass"]="Micro_MTG"
Microglia2.obs["celltype"]="MTG_"+Microglia2.obs["Supertype"].astype(str)
Microglia2.obs["sampleid"]="MTG_"+Microglia2.obs["Specimen ID"].astype(str)

ensembl_names=pd.DataFrame(Microglia2.var_names.values)
Microglia2.var.index=ensembl_names[0].replace(dt)
Microglia2.var_names=ensembl_names[0].replace(dt)
Microglia2.var_names_make_unique()



Microglia=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/Microglia_subclass_sb.h5ad')
Microglia.obs["majorclass"]="Micro_ONONH"
Microglia.obs["celltype"]="ONONH_"+Microglia.obs["subclass"].astype(str)

Microglia_all=ad.concat([Microglia1,Microglia2,Microglia])

Microglia_all.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/merged/Micro.h5ad", compression="gzip")

del Microglia_all

Oligo1=sc.read(f"{dir}/Oligodendrocyte_DLPFC.h5ad")
Oligo1=Oligo1[Oligo1.obs["assay"]=="10x 3' v3"]
Oligo1.obs["majorclass"]="Oligo_DLPFC"
Oligo1.obs["celltype"]="DLPFC_"+Oligo1.obs["Supertype"].astype(str)
Oligo1.obs["sampleid"]="DLPFC_"+Oligo1.obs["Specimen ID"].astype(str)

ensembl_names=pd.DataFrame(Oligo1.var_names.values)
Oligo1.var.index=ensembl_names[0].replace(dt)
Oligo1.var_names=ensembl_names[0].replace(dt)
Oligo1.var_names_make_unique()


Oligo2=sc.read(f"{dir}/Oligodendrocyte_MTG.h5ad")
Oligo2=Oligo2[Oligo2.obs["assay"]=="10x 3' v3"]
Oligo2.obs["majorclass"]="Oligo_MTG"
Oligo2.obs["celltype"]="MTG_"+Oligo2.obs["Supertype"].astype(str)
Oligo2.obs["sampleid"]="MTG_"+Oligo2.obs["Specimen ID"].astype(str)

ensembl_names=pd.DataFrame(Oligo2.var_names.values)
Oligo2.var.index=ensembl_names[0].replace(dt)
Oligo2.var_names=ensembl_names[0].replace(dt)
Oligo2.var_names_make_unique()


Oligo=sc.read(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_new.h5ad")

Oligo.obs["majorclass"]="Oligo_ONONH"
Oligo.obs["celltype"]="ONONH_"+Oligo.obs["subclass"].astype(str)

Oligo_all=ad.concat([Oligo1,Oligo2,Oligo])

Oligo_all.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/merged/Oligo.h5ad", compression="gzip")

del Oligo_all
