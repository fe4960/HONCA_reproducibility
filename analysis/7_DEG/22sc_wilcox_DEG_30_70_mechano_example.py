import scanpy as sc
import pandas as pd
import sys

ct=sys.argv[1]
name=sys.argv[2]
st=sys.argv[3]
gene=sys.argv[4]


sc.settings.figdir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{ct}/clean/figures/"

adata=sc.read(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{ct}/clean/{name}.h5ad")

adata=adata[adata.obs["subclass"]==st]

####adata.obs["donor_age"]=adata.obs["donor"].astype("str")+"_"+adata.obs["age"].astype("str")

#####adata = adata[adata.obs.sort_values('age').index, :]


adata.obs["age_range"]="between31n70"
adata.obs.loc[adata.obs["age_year"]<31,"age_range"]="younger31"
adata.obs.loc[adata.obs["age_year"]>70,"age_range"]="older70"
adata.obs["age_range"]=adata.obs["age_range"].astype("category")


#sc.pl.violin(adata,gene,groupby="donor_age", save=f"{gene}_age_donor.png" )

sc.pl.violin(adata,gene,groupby="age_range", save=f"{gene}_age_range.png" )

sc.pl.violin(adata,gene,groupby="age_range", jitter=False, save=f"{gene}_age_range_no_jitter.png" )
