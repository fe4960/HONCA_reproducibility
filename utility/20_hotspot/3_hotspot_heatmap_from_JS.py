import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot
import joblib
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
dir=sys.argv[1]
label=sys.argv[2]
t=sys.argv[3]
hvg=sys.argv[4]

######hs = joblib.load(f"{dir}/{label}_hs.create_modules_years.pkl")
#######adata = sc.read(f"{dir}/{label}_hs.adata_years.h5ad")


hs = joblib.load(f"{dir}/{label}_hs.create_modules_years_{hvg}.pkl")
adata = sc.read(f"{dir}/{label}_hs.adata_years_{hvg}.h5ad")

#dir = "/dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_velocity/3_Velocity/RPE/Hotsopt/"
#adata = sc.read(
#    "/dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_velocity/3_Velocity/RPE/Hotsopt/RPE_All_ntop2000_hs.adata_age.h5ad"
#)
#hs = joblib.load(
#    "/dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA_NoFetal_velocity/3_Velocity/RPE/Hotsopt/RPE_All_ntop5000_hs.create_modules_age.pkl"
#)
 
adata.obs.age=pd.to_numeric(adata.obs[t], errors='coerce')
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
#sc.pp.scale(adata)
#cols=['Unnamed: 0','age']
#disease=pd.read_csv(f"/dfs3b/ruic20_lab/jianms1/snRNA/Velocity/3_Velocity_RPE_Fetal/2_Velocity/NDims_100/cellrank/snRNA_RPE_Fetal_xpca_ndims_100_100_dynamic_latern_time.csv", header=0, index_col=0,sep=",",usecols=cols)#,keep_default_na=False)
 
#common_dis=adata.obs.index.intersection(disease.index)
#adata=adata[common_dis].copy()
#disease=disease.loc[common_dis] ## [[]] is not need. common is always an Index object
#header=disease.columns.tolist()
#for h in header:
#    adata.obs[h]=disease[h].tolist()
#    adata.obs[h]=adata.obs[h].astype(str)
 
#min_gene=250
#for min_gene in [100,160,200,250, 300, 350]:
for min_gene in [250, 300, 350]:
    print(min_gene)
    modules = hs.create_modules(min_gene_threshold=min_gene, core_only=True, fdr_threshold=0.05)
    module_scores = hs.calculate_module_scores()
    module_scores.to_csv(f"{dir}/{label}_hs_gene_module_scores_min_gene_{min_gene}_{hvg}_{t}.csv")
    #plot loca correlation
    hs.plot_local_correlations()
    fig = plt.gcf()
    fig.set_size_inches(10, 10)
    plt.savefig(
        f"{dir}/{label}_min_gene_{min_gene}_{hvg}_plot_local_correlations_{t}.tiff",
        dpi=300,
        transparent=True,
        bbox_inches="tight",
    )
    #output gene module assignment
    results = hs.results.join(hs.modules)
    results.to_csv(f"{dir}/{label}_hs_gene_module_assignment_min_gene_{min_gene}_{hvg}_{t}.csv")
    #plot top module genes in heatmap sorted by latent time
    for module in range(1, len(modules.unique())):  # Loop from 1 to 5
        # Join results and modules
        #results = hs.results.join(hs.modules)
        # Filter rows where Module equals the current module
        filtered_results = results.loc[results.Module == module]
        # Sort and get the top 10 rows
        top_results = filtered_results.sort_values('Z', ascending=False).head(20)
        # Print or process the results for the current module
        print(f"Top results for Module {module}:")
        print(top_results)
        print("-" * 40)  # Separator for clarity
        valid_genes = [gene for gene in top_results.index if gene in adata.var_names]
        #top_genes = ["C21orf58", "MIS18BP1", "LINC01572", "SMC4", "RTKN2", "ASPM", "KIF14", "APOLD1", "MKI67", "TOP2A"] + ["LINC01414", "PCDH7", "NALCN-AS1", "KIAA1217", "ROBO2", "CNTN5", "MEG8", "FIGN", "PLEKHG1", "ESRRG"]
#        kwargs = {"xticklabels": True, "colorbar": True,  "cbar_pos": "bottom" }
#        kwargs = {"colorbar": True }
        scv.pl.heatmap(adata, var_names=valid_genes, sortby='age',color_map = "RdBu_r", show=False )
# Grab the current figure and axes
        fig = plt.gcf()
        ax = plt.gca()

# Set x-tick positions and labels manually if needed
        x_ticks = np.linspace(0, adata.shape[0]-1, 10)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([str(int(x)) for x in np.linspace(adata.obs["age"].min(), adata.obs["age"].max(), 10)])

# Optionally set label and show
        ax.set_xlabel("Age")
        
        # Get the heatmap image (assumes it's the first image object in the axis)
#        im = ax.get_images()[0]

# Add the colorbar
#        cbar = plt.colorbar(im, ax=ax)
#        cbar.set_label("Expression level", fontsize=12)

#        cbar=cg.cax

#        cbar.tick_params(labelsize=10)
#        cbar.set_ylabel("Normalized expression", fontsize=12)
        plt.savefig(
            f"{dir}/{label}_gene_module_{module}_min_gene_{min_gene}_{hvg}_heatmap_{t}.tiff",
            bbox_inches="tight",
            transparent=True,
            dpi = 300
        )
#        cr.pl.heatmap(
#            adata,
#            model=model,  # use the model from before
#            lineages="Beta",
#            cluster_key="clusters",
#            show_fate_probabilities=True,
#            data_key="magic_imputed_data",
#            genes=beta_drivers.head(40).index,
#            time_key="age",
#            figsize=(12, 10),
#            show_all_genes=True,
#            weight_threshold=(1e-3, 1e-3),
#        )
