import scanpy as sc
import matplotlib.pyplot as plt
#import sys
import pandas as pd

#cell=sys.argv[1]
#label=sys.argv[2]

#indir="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_crossmap_epi"
#adata=sc.read_h5ad(f'{indir}/{cell}_combined-emb.h5ad')

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi/figures/"

label="sn_st_scvi_harmony"


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi/" #sys.argv[1]

adata=sc.read_h5ad(f"{indir}/{label}_coembed.h5ad")


dbt_table=f'{indir}/{label}_dbt.txt.gz'

dbt_full=pd.read_csv(dbt_table,sep="\t",index_col=0, header=0)

df=pd.read_csv(f'{indir}/{label}_score.txt.gz',sep="\t", header=0,index_col=0)

df1=df.max(axis=1)

cell_dbt=df1[df1<0.7].index

cells=adata.obs.index.difference(cell_dbt)
adata=adata[cells,:].copy()
#adata=adata_new

###atac_obs=pd.read_csv(f'{indir}/{cell}_atac_obs.txt.gz',sep="\t",index_col=0)

####adata.obs.loc[atac_obs.index,label]=atac_obs["lr_celltype"]
labels="lr_celltype"
adata.obs[labels]=adata.obs[labels].astype("category")
adata.obs[labels] = adata.obs[labels].cat.add_categories(["Immune_cell"])
#adata.obs[labels]=adata.obs[labels].astype(str)
adata.obs.loc[adata.obs[labels].isin(["T_cell","NK_cell", "Dendritic_cell", "B_cell", "Mast_cell"]), labels] = "Immune_cell"

adata=adata[~adata.obs[labels].isin(["Pigmented_cell","Adipocyte"])]

adata.obs[labels] = adata.obs[labels].cat.remove_unused_categories()
# Get unique conditions
conditions = adata.obs['system'].unique()
# Create figure with subplots for each condition
fig, axes = plt.subplots(1, len(conditions), figsize=(10, 5), sharex=True, sharey=True)

# Loop over conditions and plot t-SNE with CD11b color for each
for i, condition in enumerate(conditions):
    # Subset data for current condition
	adata_sub = adata[adata.obs['system'] == condition]
	sc.pl.umap(adata_sub, color=labels, size=1,  title=condition, show=False, ax=axes[i], legend_loc='on data', frameon=False)
#	sc.pl.embedding(adata_sub, basis="X_umap", color=[labels], size=10,  title=condition, frameon=False, show=False, ax=axes[i]  ,legend_loc='on data') #, palette="tab20")
    
# Save and show the figure
fig.tight_layout()
fig.savefig(f'{indir}/{label}_umap_plots.png', dpi=300)
#plt.show()

#df=pd.DataFrame({"cell" : ["AC","Astrocyte","BC","Cone","Endothelial-cell","Fibroblast","HC","Immune-cell","Macrophage","Melanocyte","MG","Microglia","Mural-cell","Oligodendrocyte","OPC","RGC","RPE","Rod", "Schwann-cell"], "mk" : ["GAD1","SLC4A11","GRM6","ARR3","PTPRB","BICC1","ONECUT1","CD69","F13A1", "MLANA","RGR","CD74","NOTCH3", "CERCAM", "SMOC1", "NEFL","RPE65","PDE6A","MPZ"]})

df=pd.DataFrame({"cell" : ["AC","Astrocyte","BC","Cone","Endothelial_cell","Fibroblast","HC", "Immune_cell", "MG", "Macrophage","Melanocyte","Microglia","Mural_cell","Oligodendrocyte","OPC","RGC","RPE","Rod", "Schwann_cell"], "mk" : ["GAD1","SLC4A11","GRM6","ARR3","PECAM1","BICC1","ONECUT1","PTPRC", "RGR", "F13A1", "MLANA","CD74","NOTCH3", "CERCAM", "SMOC1", "NEFL","RPE65","PDE6A","MPZ"]})


df = df.set_index("cell")
dt=df.to_dict()["mk"]
for i, condition in enumerate(conditions):
	adata_sub = adata[adata.obs['system'] == condition]
	sc.pl.dotplot(adata_sub, dt, groupby="lr_celltype", save=f"_{label}_{condition}.png")


adata.write(f"{indir}/{label}_coembed_flt_dbt.h5ad")
