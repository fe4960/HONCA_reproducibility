import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
import sys
cell=sys.argv[1]
rcParams["figure.figsize"] = (5,5)
dir="/dfs3b/ruic20_lab/junw42"
#obs=""
sc.settings.figdir=f"{dir}/HCA_ON/data/9_Brain_h5ad/{cell}_PFC_ONONH/figures"

#sc.settings.figdir=f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/figures"
#adata=sc.read(f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/Oligo_GSE267301_merge_hvg10000_epochnone_seurat_v3_rs_0.4_rm3_4_5_none_scvi_trg.h5ad")
adata=sc.read(f"{dir}/HCA_ON/data/9_Brain_h5ad/{cell}_PFC_ONONH/{cell}_PFC_ONONH_hvg_10000_epoch_none_seurat_v3_rs_1_sb_scvi_trg.h5ad")
#bname="Oligo_GSE267301_merge_rm4_rm3_4_5"
bname=f"{cell}_PFC_ONONH_hvg_10000_epoch_none_seurat_v3_rs_0.4_sb"
#adata=sc.read(f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/Oligo_GSE267301_merge_hvg_10000_epoch_none_seurat_v3_rs_1_scvi_trg.h5ad")
adata.obs["source"]="ONONH"
#adata.obs.loc[adata.obs["sampleid"].isin([,"source"]="GSE267301"
adata.obs.loc[adata.obs["sampleid"].str.contains("GSM"),"source"]="GSE167494"

dt={"ONONH": "lightyellow", "GSE267301": "darkblue","GSE167494": "red"}
sc.pl.umap(adata,color="source",save=f"_{bname}_source.png", palette=dt,size=1)
sc.pl.umap(adata,color="subclass",save=f"_{bname}_subclass.png", size=1)

#meta=f"{dir}/HCA_ON/data/scvi/RNA/raw/SraRunTable-2.csv"
#me=pd.read_csv(meta,header=0,index_col=0)

#dt=me.to_dict()

#dt["brain_region"]
adata.obs["region"]="PFC"
#adata.obs["region"]=adata.obs["sampleid"].replace(dt["brain_region"])
#adata.obs["region"]=adata.obs["region"].cat.add_categories(["ON","ONH"])
adata.obs.loc[adata.obs["sampleid"].str.contains("ONH"),"region"]="ONH"
adata.obs.loc[adata.obs["sampleid"].str.contains("_ON_"),"region"]="ON"
#me1=pd.read_csv(f"{dir}/HCA_ON/data/scvi/RNA/raw/Josh_meta",header=None,index_col=0,sep="\t")
#dt1=me1.to_dict()
#adata.obs["region"]=adata.obs["region"].replace(dt1[7])

#for c in adata.obs["region"].value_counts().index:
sc.pl.umap(adata,color="region",save=f"_{bname}_region_PFC.png",groups=["PFC"], size=1, na_color="lightyellow",frameon=False)
#sc.pl.umap(adata,color="region",save=f"_{bname}_region_PRC.png",groups=["Parietal cortex"], size=1, alpha=0.7, na_color="whitesmoke", frameon=False)
#sc.pl.umap(adata,color="region",save=f"_{bname}_region_MTL.png",groups=["Medial temporal lobe"], size=1, alpha=0.7, na_color="whitesmoke", frameon=False)
#sc.pl.umap(adata,color="region",save=f"_{bname}_region_OC.png",groups=["Occipital cortex"], size=1, alpha=0.7, na_color="whitesmoke", frameon=False)
sc.pl.umap(adata,color="region",save=f"_{bname}_region.png", size=1, alpha=0.7, frameon=False)
sc.pl.umap(adata,color="region",save=f"_{bname}_region_ON.png",groups=["ON"], size=1, na_color="lightyellow", frameon=False)
sc.pl.umap(adata,color="region",save=f"_{bname}_region_ONH.png",groups=["ONH"], size=1, na_color="lightyellow", frameon=False)

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

#adata.obs["disease"]="Healthy_eye"
#adata.obs.loc[adata.obs["sampleid"].str.contains("SRR"),"disease"]=adata.obs["sampleid"].replace(donor)
#sc.pl.umap(adata,color="disease",save=f"_{bname}_disease.png",groups=["ALSP","AD","Ctr"], size=1, alpha=0.7, frameon=False)

#obs=pd.read_csv(f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/Oligo_GSE267301_merge_hvg_10000_epoch_none_seurat_v3_rs_0.4_res_0.4.obs.gz",header=0,index_col=0)

obs=pd.read_csv(f"{dir}/HCA_ON/data/9_Brain_h5ad/{cell}_PFC_ONONH/{cell}_PFC_ONONH_hvg_10000_epoch_none_seurat_v3_rs_0.4_sb_res_0.4.obs.gz",header=0,index_col=0)


#obs=pd.read_csv(f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/{obs}",header=0,index_col=0)

adata.obs["leiden1"]=obs.loc[adata.obs.index,"leiden"]

adata.obs[["subclass","leiden1","nCount_RNA"]].groupby(["subclass","leiden1"]).count().to_csv(f"{dir}/HCA_ON/data/9_Brain_h5ad/{cell}_PFC_ONONH/{bname}_leiden_subclass")

#adata[adata.obs["region"].isin(["ON","ONH","Prefrontal cortex"])].obs[["disease","leiden","nCount_RNA"]].groupby(["disease","leiden"]).count().to_csv(f"{dir}/HCA_ON/data/scvi/RNA/raw/Oligo_GSE267301_merge/{bname}_leiden_disease_ONONH_PFC")
