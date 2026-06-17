import scanpy as sc
import senepy as sp
from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/figures"
hubs = sp.load_hubs(species = 'Human')
hubs.metadata
#Look at hubs within Lung
hubs.metadata[hubs.metadata.tissue == 'Lung'] #capitalization matters

#hubs.hubs[('Lung', 'fibroblast', 1)][:10] #only showing first 10 with [:10]

#hubs.get_genes(('Lung', 'fibroblast', 1))[:10]

#filt_meta = hubs.metadata[hubs.metadata.tissue == 'Lung']
#filt_meta

#hubs.merge_hubs(filt_meta, new_name = 'Lung_merged')

#len(hubs.hubs['Lung_merged'])

#hubs.hubs['Lung_merged'][:10]

#hubs.merge_hubs(filt_meta, new_name = 'Lung_merged_min2', overlap_threshold = 2)

#len(hubs.hubs['Lung_merged_min2'])

hubs.merge_hubs(hubs.metadata, new_name = 'universal', calculate_thresh = True, p_thres = 0.01)

len(hubs.hubs['universal'])

#hubs.hubs['universal'][:10]

#sorted(hubs.hubs['universal'], key=lambda x: x[1], reverse=True)[:10]

#hubs.literature_markers

#list(hubs.literature_markers)[:10]

#hubs.senGPT[0:10]

#hubs.search_hubs_by_genes(['Cdkn2a', 'Cdkn1a', 'Il6', 'Cxcl13']).head() #remove .head() for all results

#hubs.search_hubs_by_genes(hubs.senGPT).head()

adata = sc.read_h5ad('/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_Fibroblast_rawcount.h5ad')

translator = sp.translator(hub = hubs.hubs, data = adata)

#fib_subset = adata[adata.obs['cell_ontology_class'] == 'fibroblast of lung'].copy()

#fib_subset.obs['sen_score'] = sp.score_hub(fib_subset, hubs.hubs[('Lung', 'fibroblast', 0)])

#adata.obs['universal_score'] = sp.score_all_cells(adata, hubs.hubs["universal"],
#                   identifiers = ['majorclass'])

adata.obs['universal_score'] = sp.score_hub(adata, hubs.hubs["universal"], translator = translator)


adata.write('/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/ONONH_Fibroblast_rawcount_sene.h5ad')
adata.obs.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/ONONH_Fibroblast_rawcount_sene.obs.gz")

adata.obs.universal_score.hist(bins = 50)

sc.pl.umap(adata, color='universal_score', vmax = 8, save="Fibroblast_universal_score.png")

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

sc.pl.umap(adata, color='putative_sen', save="Fibroblast_putative_sen.png")

#plt.figure(figsize = (3,3))
#ax = sns.stripplot(adata.obs, x = 'Condition', y = 'Hepatocyte_0')
#import seaborn as sns
#from scipy import stats
#import matplotlib.pyplot as plt
#sns.stripplot(adata.obs, x = 'age_year', y = 'universal_score')
import pandas as pd
prop=obs[["putative_sen","donor","celltype"]].groupby(["donor","celltype"]).sum()/obs[["putative_sen","donor","celltype"]].groupby(["donor","celltype"]).count()
prop.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/ONONH_Fibroblast_rawcount_sene_prop",sep="\t")
prop=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/ONONH_Fibroblast_rawcount_sene_prop",header=0,index_col=0,sep="\t")
age=obs[["donor","age_year"]].drop_duplicates()
age=age.set_index("donor")
prop["age"]=age.loc[prop.index,"age_year"]
prop.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/ONONH_Fibroblast_rawcount_sene_prop_age",sep="\t")
#compare % of sen cells at 3m vs 21m
#m3_sub = adata[adata.obs.age == '3m']
#m21_sub = adata[adata.obs.age == '21m']

#m3_sub.obs.putative_sen.sum() / len(m3_sub) * 100

#m21_sub.obs.putative_sen.sum() / len(m21_sub) * 100


