import os
import tempfile
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import sys
# Reproducibility
import scvi
import anndata
from scvi.external import SysVI

scvi.settings.seed = 7
print("Last run with scvi-tools version:", scvi.__version__)

ref=sys.argv[1]
query_list=sys.argv[2]
syst=sys.argv[3]
outdir=sys.argv[4]
label=sys.argv[5]
cl=sys.argv[6]
sc.settings.figdir = f"{outdir}/figures"

adata_ref = sc.read(ref)

adata_ref.obs["system"]="snRNA"
adata_ref.obs["ident"]="ref"

#query_list=pd.read(query)

sample=[]

with open(query_list, "r") as sl:
    for line in sl:
        line1=line.strip()
        sample.append(line1)

print(f'done reading list')
adata_list=[]

for sam in sample:
    print(f'read {sam}')
    path=f'{sam}'
    if os.path.exists(path):
        adata=sc.read(path)
    else :
        path=f'{sam}'
        adata=sc.read(path)
#    adata.obs["sample"]=sam
#    print(f'mody index')
#    if sam not in array:
#        adata.obs.index=sam+"_"+adata.obs.index
    adata.obs["system"]=syst
    adata.obs[cl]="unknown"
    adata.obs["ident"]="query"
    adata_list.append(adata)
    print(f'done reading {sam}')

adata_query=anndata.concat(adata_list)

gene=adata_query.var.index.intersection(adata_ref.var.index)

adata_query=adata_query[:,gene]
adata_ref=adata_ref[:,gene]

adata=anndata.concat([adata_query,adata_ref])


adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata,
        n_top_genes=100000,
        batch_key = "sampleid",
        subset=False)
#        span=span )


SysVI.setup_anndata(
    adata=adata,
    batch_key="system",
    categorical_covariate_keys=["sampleid"],
)

# Example showing how to turn on embedding of all categorical covariates
model = SysVI(
    adata=adata,
    embed_categorical_covariates=True,
)

# Example showing how to adjust loss weights
#model.train(
#    plan_kwargs={
#        "kl_weight": 1,
#        "z_distance_cycle_weight": 0
#        # Add additional parameters, such as number of epochs
#    }
#)

# Change the seed
scvi.settings.seed = 7

# Now initialise and train the model
#model = SysVI(adata=adata)
#model.train()

# Initialise the model
#model = SysVI(adata=adata)
# Train
max_epochs = 200
model.train(
    max_epochs=max_epochs, check_val_every_n_epoch=1, plan_kwargs={"z_distance_cycle_weight": 0,  "kl_weight": 1}
)
adata.write(f"{outdir}/{label}_sysVI.h5ad")

# Plot loses
# The plotting code below was specifically adapted to the above-specified model and its training
# If changing the model or training the plotting functions may need to be adapted accordingly

# Make detailed plot after N epochs
epochs_detail_plot = 100

# Losses to plot
losses = [
    "reconstruction_loss_train",
    "kl_local_train",
    "cycle_loss_train",
]
fig, axs = plt.subplots(2, len(losses), figsize=(len(losses) * 3, 4))
for ax_i, l_train in enumerate(losses):
    l_val = l_train.replace("_train", "_validation")
    l_name = l_train.replace("_train", "")
    # Change idx of epochs to start with 1
    l_val_values = model.trainer.logger.history[l_val].copy()
    l_val_values.index = l_val_values.index + 1
    l_train_values = model.trainer.logger.history[l_train].copy()
    l_train_values.index = l_train_values.index + 1
    for l_values, c, alpha, dp in [
        (l_train_values, "tab:blue", 1, epochs_detail_plot),
        (l_val_values, "tab:orange", 0.5, epochs_detail_plot),
    ]:
        axs[0, ax_i].plot(l_values.index, l_values.values.ravel(), c=c, alpha=alpha)
        axs[0, ax_i].set_title(l_name)
        axs[1, ax_i].plot(l_values.index[dp:], l_values.values.ravel()[dp:], c=c, alpha=alpha)

fig.tight_layout()
fig.savefig(f"{outdir}/figures/training_qc_{label}.png")

# Get embedding - save it into X of new AnnData
embed = model.get_latent_representation(adata=adata)
embed = sc.AnnData(embed, obs=adata.obs)
# Make system categorical for plotting below
#embed.obs["system"] = embed.obs["system"].map({0: "mouse", 1: "human"})

# Compute UMAP
sc.pp.neighbors(embed, use_rep="X")
sc.tl.umap(embed)


embed.write(f"{outdir}/{label}_sysVI_embed.h5ad")

# Plot UMAP embedding

# Obs columns to color by
cols = ["system", cl ]

# One plot per obs column used for coloring
fig, axs = plt.subplots(len(cols), 1, figsize=(3, 3 * len(cols)))
for col, ax in zip(cols, axs, strict=False):
    sc.pl.embedding(
        embed,
        "X_umap",
        color=col,
        s=10,
        ax=ax,
        show=False,
        sort_order=False,
        frameon=False,
        save=f"_{col}_{label}.png"
    )


    # Plot UMAP embedding per system
systems = sorted(embed.obs.system.unique())
ncols = len(systems)
# Plot systms side by side
fig, axs = plt.subplots(1, ncols, figsize=(3 * ncols, 3))
for i, system in enumerate(systems):
    ax = axs[i]
    # Plot all cells as background and on top cells from one system colored by cell type
    sc.pl.umap(embed, ax=ax, show=False, s=5, frameon=False)
    sc.pl.umap(
        embed[embed.obs.system == system, :],
        color=cl,
        ax=ax,
        show=False,
        s=10,
        frameon=False,
        title=system,
#        save="_{systems}_{label}.png"
    )
    # Keep legend only on the last plot (assuming this legend contains all categories)
    if i != ncols - 1:
        ax.get_legend().remove()

fig.tight_layout()
fig.savefig(f"{outdir}/figures/umap_{label}_facet.png")
