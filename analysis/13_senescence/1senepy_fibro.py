import scanpy as sc
#import senepy as sp
from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/figures"

adata=sc.read('/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/ONONH_Fibroblast_rawcount_sene.h5ad')
#adata.obs.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/ONONH_Fibroblast_rawcount_sene.obs.gz")

#adata.obs.universal_score.hist(bins = 50)

#sc.pl.umap(adata, color='universal_score', vmax = 8, save="Fibroblast_universal_score.png")

#adata.obs.groupby('age_year').count()[['cell']]

e = adata.obs.universal_score.mean() #distribution mean
std = adata.obs.universal_score.std() #distribution std
thresh = e + 3*std
thresh

#function to add senescent label
def is_putative_sen(x):
    if x >= thresh:
        return 1
    else:
        return 0


#map function to a new row in adata.obs
adata.obs['putative_sen'] = adata.obs.universal_score.map(is_putative_sen)

#sc.pl.umap(adata, color='putative_sen', save="Fibroblast_putative_sen.png")

#plt.figure(figsize = (3,3))
#ax = sns.stripplot(adata.obs, x = 'Condition', y = 'Hepatocyte_0')
#import seaborn as sns
#from scipy import stats
#import matplotlib.pyplot as plt
#sns.stripplot(adata.obs, x = 'age_year', y = 'universal_score')
import pandas as pd
prop = adata.obs.groupby("donor")["putative_sen"].mean()
count=adata.obs.groupby("donor")["putative_sen"].count()
age=adata.obs[["donor","age_year"]].drop_duplicates()
age=age.set_index("donor")

sene=pd.DataFrame({ "cell_count": count, "sen_prop": prop})

sene = sene.join(age)

sene1 = sene.loc[sene["cell_count"]>=100,]

import seaborn as sns
import matplotlib.pyplot as plt

# Optionally bin the age if it's numeric and you want grouped boxes
# prop["age_bin"] = pd.cut(prop["age"], bins=[20,30,40,50,60,70,80], right=False)

# Basic boxplot
sns.boxplot(x="age_year", y="sen_prop", data=sene1)

plt.xlabel("Age (years)")
plt.ylabel("Proportion of Putative Senescent Cells")
plt.title("Putative Senescence vs. Age")
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()


##########
prop = adata.obs.groupby(["donor","celltype"])["putative_sen"].mean()
count=adata.obs.groupby(["donor","celltype"])["putative_sen"].count()
age=adata.obs[["donor","age_year"]].drop_duplicates()
age=age.set_index("donor")

sene=pd.DataFrame({ "cell_count": count, "sen_prop": prop})

sene = sene.join(age)

sene1 = sene.loc[sene["cell_count"]>=100,]
sene1 = sene1.reset_index()
import seaborn as sns
import matplotlib.pyplot as plt

# Optionally bin the age if it's numeric and you want grouped boxes
# prop["age_bin"] = pd.cut(prop["age"], bins=[20,30,40,50,60,70,80], right=False)

# Basic boxplot

# Facet boxplot: one boxplot per cell type
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
plt.show()

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

#prop.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/ONONH_Fibroblast_rawcount_sene_prop_age",sep="\t")
#compare % of sen cells at 3m vs 21m
#m3_sub = adata[adata.obs.age == '3m']
#m21_sub = adata[adata.obs.age == '21m']

#m3_sub.obs.putative_sen.sum() / len(m3_sub) * 100

#m21_sub.obs.putative_sen.sum() / len(m21_sub) * 100


