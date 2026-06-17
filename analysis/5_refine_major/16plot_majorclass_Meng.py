import scanpy as sc
import pandas as pd
dir="/dfs3b/ruic20_lab/junw42"
f1="hvg2000_nonsb_epochnone_seurat_v3_rs"
bname=f"_major_{f1}_1_clean_ON_ONH_new"
adata=sc.read(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/major_ON_ONH_new_final.h5ad")
adata.obs["majorclass1"]=adata.obs["majorclass1"].replace("Unknown?","Immune_neuron_unk")

sc.settings.figdir = f"{dir}/HCA_ON/data/5_refine_major/scvi/major/figures"

lis=pd.read_csv(f'{dir}/HCA_ON/scripts/5_refine_major/Meng_mk',header=None)
dt=lis[0].values #.set_index(0).to_dict()[1]
sc.pl.dotplot(adata, dt, groupby='majorclass1', save=f"_ON_ONH_ODD.png")


#adata=sc.read()
