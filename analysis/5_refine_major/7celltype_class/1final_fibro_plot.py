import scanpy as sc
import pandas as pd

sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/figures/"

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new.h5ad")

sc.pl.dotplot(adata,["PIEZO2","PIEZO1","TRPV4","GJA1","RAPGEF3","KCNK2", "KCNK15", "KCNK4", "KCNK3","KCNK5","KCNK9","KCNK13","KCNK12","KCNK2","KCNK10","KCNK18","KCNK1","KCNK6", "TRPV1", "TRPV2", "TRPV4", "TRPC1", "TRPC3", "TRPC5", "TRPC6", "TRPA1", "PKD2", "TRPM3", "TRPM4", "TRPM7"],groupby="subclass", save="Fibroblast_subclass_new5_clean_mechanochannel.png" )

adata.obs["age_range"]="31-70"
adata.obs.loc[adata.obs["age_year"]<31,"age_range"]="<31"
adata.obs.loc[adata.obs["age_year"]>70,"age_range"]=">70"
sc.pl.violin(adata[adata.obs["subclass"].str.contains("ONHON")],"PIEZO2",groupby="age_range",jitter=False)

import scipy.stats as stats


outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/"

for ct in adata.obs["subclass"].unique():

    bname=f"Fibroblast_subclass_clean_new_{ct}_70_30"

    adata1=adata[(adata.obs["subclass"]==ct)&(adata.obs["age_range"].isin(["<31",">70"]))]

    sc.tl.rank_genes_groups(adata1, "age_range", method="wilcoxon", use_raw=False)

    deg1={ "score": list(pd.DataFrame(adata.uns['rank_genes_groups']['scores'])[">70"]),
    'pvals': list(pd.DataFrame(adata.uns['rank_genes_groups']['pvals'])[">70"]),
    'pvals_adj': list(pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj'])[">70"]),
    'logfoldchanges': list(pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges'])[">70"])}

    data2=pd.DataFrame(data=deg1,index=pd.DataFrame(adata.uns['rank_genes_groups']['names'])[">70"])
    data2.to_csv(f"{outdir}/{bname}_deg",sep="\t")
