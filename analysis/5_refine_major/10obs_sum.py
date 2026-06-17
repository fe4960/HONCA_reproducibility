import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

rcParams["figure.figsize"] = (5,5)
rcParams.update({'font.size': 14})


meta=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid", sep=",", header=0)
meta=meta.set_index("sampleid")



def generate_txt(outdir, bname):
    adata=sc.read(f'{outdir}/{bname}.h5ad') 
    adata=adata[(adata.obs["sampleid"]!="MMD_23_17976_ONH_RNA")&(adata.obs["sampleid"]!="MMD_23_20181_ONH_RNA")]
    obs=adata.obs.copy()
    dt={"Asian-Indian": "Asian"}
    mask=obs["source"] == "Sanes"
    obs.loc[mask, "race"] = obs.loc[mask, "sampleid"].map(meta["race"])
    obs["race"]=obs["race"].astype("str").replace(dt)
#    obs["race"]=obs["race"].replace(dt)
#    obs["race"] = obs["race"].cat.rename_categories({"Asian-Indian": "Asian"})
    obs[["nCount_RNA","subclass","gender"]].groupby(["subclass","gender"]).count().to_csv(f"{outdir}/{bname}_gender.csv",sep="\t")
    obs[["nCount_RNA","subclass","race"]].groupby(["subclass","race"]).count().to_csv(f"{outdir}/{bname}_race.csv",sep="\t")
    sc.pl.embedding(adata, basis="X_umap", color=["sampleid"],ncols=1,frameon=False,save=f'_{bname}_sampleid.png', palette="tab20")






outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
bname="Astrocyte_subclass_new5_clean"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/"
bname="Oligodendrocyte_precursor_cell_subclass_seurat_cycling"
generate_txt(outdir, bname)

outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"
bname="Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location"
generate_txt(outdir, bname)

outdir=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/'
bname='Endothelial_cell_subclass_sb_seurat_rmRPE'
generate_txt(outdir, bname)


outdir=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Mural_cell/clean/'
bname="Mural_cell_subclass_rmRPE"
generate_txt(outdir, bname)

outdir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/"
bname="Oligodendrocyte_subclass_seurat"
generate_txt(outdir, bname)

outdir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Macrophage/clean/"
bname="Macrophage_subclass_sb_clean_rmRPE"
generate_txt(outdir, bname)

outdir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/"
bname="Microglia_subclass_sb_clean"
generate_txt(outdir, bname)









