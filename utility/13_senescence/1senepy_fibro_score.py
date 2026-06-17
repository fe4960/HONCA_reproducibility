import scanpy as sc
import senepy as sp
from matplotlib import rcParams
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

rcParams["figure.figsize"] = (5,5)

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/sene"
hubs = sp.load_hubs(species = 'Human')
hubs.metadata
#Look at hubs within Lung

adata = sc.read_h5ad('/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_Fibroblast_rawcount.h5ad')

translator = sp.translator(hub = hubs.hubs, data = adata)


fibro_metadata=hubs.metadata[hubs.metadata.cell == 'fibroblast'] #capitalization matters

hubs.merge_hubs(fibro_metadata, new_name = 'fibro_merged', calculate_thresh = False)


adata.obs['fibro_merged_score'] = sp.score_hub(adata, hubs.hubs["fibro_merged"], translator = translator)


hubs.merge_hubs(hubs.metadata, new_name = 'universal_thre', calculate_thresh = True, p_thres = 0.01)
adata.obs['universal_thre_score'] = sp.score_hub(adata, hubs.hubs["universal_thre"], translator = translator)


hubs.merge_hubs(hubs.metadata, new_name = 'universal', calculate_thresh = False)
adata.obs['universal_score'] = sp.score_hub(adata, hubs.hubs["universal"], translator = translator)



#function to add senescent label
def is_putative_sen(x):
    if x >= thresh:
        return 1
    else:
        return 0


#map function to a new row in adata.obs

e = adata.obs.universal_score.mean() #distribution mean
std = adata.obs.universal_score.std() #distribution std
thresh = e + 3*std
thresh


adata.obs['sene_universal'] = adata.obs.universal_score.map(is_putative_sen)

e = adata.obs.universal_thre_score.mean() #distribution mean
std = adata.obs.universal_thre_score.std() #distribution std
thresh = e + 3*std
thresh


adata.obs['sene_universal_thre'] = adata.obs.universal_thre_score.map(is_putative_sen)

e = adata.obs.fibro_merged_score.mean() #distribution mean
std = adata.obs.fibro_merged_score.std() #distribution std
thresh = e + 3*std
thresh



adata.obs['sene_fibro'] = adata.obs.fibro_merged_score.map(is_putative_sen)

main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/sene/"

adata.write(f'{main}/ONONH_Fibroblast_rawcount_sene.h5ad')

adata.obs.to_csv(f"{main}/ONONH_Fibroblast_rawcount_sene.obs.gz")



sene=["universal_score","universal_thre_score","fibro_merged_score"]

for s in sene:
    adata.obs[s].hist(bins=50)

    # Save the figure
    plt.savefig(f"{main}/Fibro_{s}_histogram.png", dpi=300, bbox_inches='tight')

    # Optionally show or close the plot
    plt.show()
    plt.close()  # if you want to close it right after saving
    sc.pl.umap(adata, color=s, vmax = 8, save=f"Fibro_{s}.png")



score=['sene_universal', 'sene_universal_thre', 'sene_fibro']
for s in score:
    sc.pl.umap(adata, color=s, save=f"Fibro_{s}.png")

    # Plot the histogram

    prop = adata.obs.groupby("donor")[s].mean()
    count=adata.obs.groupby("donor")[s].count()
    age=adata.obs[["donor","age_year"]].drop_duplicates()
    age=age.set_index("donor")

    sene=pd.DataFrame({ "cell_count": count, "sen_prop": prop})

    sene = sene.join(age)

    sene1 = sene.loc[sene["cell_count"]>=100,]


# Optionally bin the age if it's numeric and you want grouped boxes
# prop["age_bin"] = pd.cut(prop["age"], bins=[20,30,40,50,60,70,80], right=False)

# Basic boxplot
    sns.boxplot(x="age_year", y="sen_prop", data=sene1)

    plt.xlabel("Age (years)")
    plt.ylabel("Proportion of Putative Senescent Cells")
    plt.title("Putative Senescence vs. Age")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(f"{main}/fibro_{s}_all_box_plot.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close("all")

    g=sns.lmplot(
        data=sene1,
        x="age_year",
        y="sen_prop",
        scatter_kws={"s": 20, "alpha": 0.6},
        line_kws={"color": "red"},
        ci=None           # remove confidence interval for clean look
    )
    g.set_axis_labels("Age (years)", "Proportion of Putative Senescent Cells")
    g.fig.suptitle("Putative Senescence vs. Age", y=1.03)  # Shift title up slightly

# Rotate x-axis tick labels if needed
    for ax in g.axes.flat:
        for label in ax.get_xticklabels():
            label.set_rotation(90)

    plt.tight_layout()
    plt.savefig(f"{main}/fibro_{s}_all_lm_plot.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close("all")

   
##########
    prop = adata.obs.groupby(["donor","celltype"])[s].mean()
    count=adata.obs.groupby(["donor","celltype"])[s].count()
    age=adata.obs[["donor","age_year"]].drop_duplicates()
    age=age.set_index("donor")

    sene=pd.DataFrame({ "cell_count": count, "sen_prop": prop})

    sene = sene.join(age)

    sene1 = sene.loc[sene["cell_count"]>=100,]
    sene1 = sene1.reset_index()

    sns.catplot(
        data=sene1,
        x="age_year",           # X-axis variable
        y="sen_prop",  # Y-axis variable
        col="celltype",    # Facet by 'celltype'
        kind="box",
        col_wrap=3,        # Number of columns per row
        height=4,
        aspect=1
    )

    plt.tight_layout()
    plt.savefig(f"{main}/fibro_{s}_celltype_box_plot.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close("all")

# Facet plot: one plot per cell type
    sns.lmplot(
        data=sene1,
        x="age_year",
        y="sen_prop",
        col="celltype",
        col_wrap=3,       # Wrap plots to multiple rows
        height=4,
        aspect=1,
        scatter_kws={"s": 20, "alpha": 0.6},
        line_kws={"color": "red"},
        ci=None           # remove confidence interval for clean look
    )
    plt.tight_layout()
    plt.savefig(f"{main}/fibro_{s}_celltype_lm_plot.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close("all")    



# Optionally bin the age if it's numeric and you want grouped boxes
# prop["age_bin"] = pd.cut(prop["age"], bins=[20,30,40,50,60,70,80], right=False)

# Basic boxplot

# Facet boxplot: one boxplot per cell type




#adata.obs.groupby('age_year').count()[['cell']]



