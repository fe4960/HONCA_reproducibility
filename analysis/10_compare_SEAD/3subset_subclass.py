import scanpy as sc
import anndata as ad
import pandas as pd
import sys

cell=sys.argv[1]
h5ad=sys.argv[2]
dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/SEAD"


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/SEAD/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad")

Astro1=adata[(adata.obs["method"]=="10Xv3.1") & (adata.obs["Subclass"]==cell)]

#adata=adata[ (adata.obs["Subclass"]=="Astrocyte")]


df=pd.read_csv("/dfs3b/ruic20_lab/junw42/reference/cellranger/refdata-gex-GRCh38-2020-A_ensemblID_genename",header=None,index_col=0,sep="\t")
dt=df.to_dict()[1]


ensembl_names=pd.DataFrame(Astro1.var_names.values)
Astro1.var.index=ensembl_names[0].replace(dt)
Astro1.var_names=ensembl_names[0].replace(dt)
Astro1.var_names_make_unique()

Astro1.obs["majorclass"]=f"{cell}_PFC"
Astro1.obs["celltype"]="PFC_"+Astro1.obs["Supertype"].astype(str)
Astro1.obs["sampleid"]="PFC_"+Astro1.obs["sample_name"].astype(str)
Astro1.layers['counts']=Astro1.layers['UMIs'].copy()
del Astro1.layers['UMIs']
del adata


adata1=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/SEAD/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad")


#adata.write(f'{dir}/Astro_PFC.h5ad')


Astro2=adata1[(adata1.obs["method"]=="10Xv3.1") & (adata1.obs["Subclass"]==cell)]

#adata=adata[ (adata.obs["Subclass"]=="Astrocyte")]


df=pd.read_csv("/dfs3b/ruic20_lab/junw42/reference/cellranger/refdata-gex-GRCh38-2020-A_ensemblID_genename",header=None,index_col=0,sep="\t")
dt=df.to_dict()[1]


ensembl_names=pd.DataFrame(Astro2.var_names.values)
Astro2.var.index=ensembl_names[0].replace(dt)
Astro2.var_names=ensembl_names[0].replace(dt)
Astro2.var_names_make_unique()

Astro2.obs["majorclass"]="{cell}_PFC"
Astro2.obs["celltype"]="PFC_"+Astro2.obs["Supertype"].astype(str)
Astro2.obs["sampleid"]="PFC_"+Astro2.obs["sample_name"].astype(str)
Astro2.layers['counts']=Astro2.layers['UMIs'].copy()
del Astro2.layers['UMIs']
del adata1

Astro=sc.read(h5ad)

#Astro=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new2.h5ad')
Astro.obs["majorclass"]=f"{cell}_ONONH"
Astro.obs["celltype"]="ONONH_"+Astro.obs["subclass"].astype(str)


Astro_all=ad.concat([Astro1,Astro2,Astro])

Astro_all.write(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/merged/{cell}.h5ad", compression="gzip")


