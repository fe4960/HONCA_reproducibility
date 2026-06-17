import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
#from matplotlib import rcParams
sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/figures/"
########
#outdir="/dfs3b/ruic20_lab/junw42/DRG/data/scvi/raw/major/"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/"
#bname="major_DRG_horse_ultima_merged_hvg10000_epochnone_seurat_rs1_sd7_sb_scvi1_2_rm_none_500gene_anno_flt"
bname="major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_rmUnk"
#adata=sc.read("/dfs3b/ruic20_lab/junw42/DRG/data/scvi/raw/major/major_DRG_horse_ultima_merged_hvg10000_epochnone_seurat_rs1_sd7_sb_scvi1_2_rm_none_500gene_scvi_trg_anno_flt.h5ad")
#adata=sc.read("/dfs3b/ruic20_lab/junw42/DRG/data/scvi/raw/major/major_DRG_horse_ultima_merged_hvg10000_epochnone_seurat_rs1_sd7_sb_scvi1_2_rm_none_500gene_scvi_trg_anno_flt_final_new.h5ad")
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_scvi_trg.h5ad")
adata=adata[adata.obs["celltype"]!="Unknown?"]
adata.write(f"{outdir}/clean/{bname}.h5ad")
#dt={"Adipocytes": "Adipocyte", "Lymphatic_endothelium": "Lymphatic_Endothelium"}

#adata.obs["cell_class"]=adata.obs["cell_class1"].replace(dt)

#adata.obs.to_csv(f"{outdir}/{bname}.obs.gz",sep="\t")

#adata.write(f"/dfs3b/ruic20_lab/junw42/DRG/data/scvi/raw/major/{bname}.h5ad")

#sc.pl.umap(adata,color="cell_class", frameon=False, save="_{bname}.png")


cell_counts = adata.obs['majorclass1'].value_counts()

# 2. Get the sorted categories (most to least frequent)
sorted_classes = cell_counts.index.tolist()
categories=cell_counts.index
# 3. Reorder the categorical
adata.obs['majorclass1'] = pd.Categorical(
    adata.obs['majorclass1'],
    categories=sorted_classes,
    ordered=True
)

plt.rcParams["figure.figsize"] = (6,6)
from matplotlib.colors import to_hex
#colors = cm.get_cmap('tab20')(range(len(categories)))
#colors = cm.get_cmap('tab20')(range(len(categories)))
cmap = cm.get_cmap('tab20')
colors = [to_hex(cmap(i)) for i in range(len(categories))]
sc.pl.umap(adata,color="majorclass1",save=f"_{bname}_umap.png", frameon=False, palette=colors,  size=1)


values = cell_counts
# Create the plot
fig=plt.figure(figsize=(5, 5))
colors = cm.get_cmap('tab20')(range(len(categories)))
plt.bar(categories, values, color=colors)

# Add titles and labels
plt.title('Cell class distribution',fontsize=15)
plt.xlabel('Cell class',fontsize=15)
plt.ylabel('Number of cells',fontsize=15)
plt.xticks(rotation=90)

# Show the plot
plt.tight_layout()
#plt.show()

plt.savefig(f"{outdir}/figures/{bname}_celltype_bar_plot.pdf", format='pdf', bbox_inches='tight')
plt.show()
##########ntrg=20
##########sc.tl.rank_genes_groups(adata, "majorclass1", method="wilcoxon", use_raw=True)
#########sc.pl.rank_genes_groups(adata, n_genes=ntrg,sharey=False, save=f'_{bname}_TRG.png')
#########sc.pl.rank_genes_groups_dotplot(adata, n_genes=ntrg,  save=f'_{bname}_TRG.png', dendrogram=False)
#########sc.pl.rank_genes_groups_dotplot(adata, n_genes=2,  save=f'_{bname}_2TRG.png', dendrogram=False)

########top_genes = pd.DataFrame({
#######    group: adata.uns['rank_genes_groups']['names'][group][:20]
########    for group in adata.uns['rank_genes_groups']['names'].dtype.names
##########    })
########print(top_genes)

#########top_genes.to_csv(f"{outdir}/{bname}_top20.txt",sep="\t")

