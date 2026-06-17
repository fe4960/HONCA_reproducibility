import scanpy as sc
import pandas as pd
import sys

ct=sys.argv[1]
name=sys.argv[2]



sc.settings.figdir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{ct}/clean/figures/"

adata=sc.read(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{ct}/clean/{name}.h5ad")

mk=["P2RX1","P2RX2","P2RX3","P2RX4","P2RX5","P2RX6","P2RX7", "PIEZO2","PIEZO1","TRPV4","GJA1","RAPGEF3","KCNK2", "KCNK15", "KCNK4", "KCNK3","KCNK5","KCNK9","KCNK13","KCNK12","KCNK2","KCNK10","KCNK18","KCNK1","KCNK6", "TRPV1", "TRPV2", "TRPV4", "TRPC1", "TRPC3", "TRPC5", "TRPC6", "TRPA1", "PKD2", "TRPM3", "TRPM4", "TRPM7", "PKD1", "TMEM63A", "TMEM63B", "TMEM63C"]

sc.pl.dotplot(adata, mk, groupby="subclass", save=f"{ct}_subclass_new5_clean_mechanochannel.png" )

adata.obs["age_range"]="between31n70"
adata.obs.loc[adata.obs["age_year"]<31,"age_range"]="younger31"
adata.obs.loc[adata.obs["age_year"]>70,"age_range"]="older70"
adata.obs["age_range"]=adata.obs["age_range"].astype("category")
sc.pl.violin(adata[adata.obs["subclass"].str.contains("ONHON")],"PIEZO2",groupby="age_range",jitter=False)

import scipy.stats as stats


outdir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{ct}/clean/"

for ct1 in adata.obs["subclass"].unique():

    bname=f"{name}_{ct1}_70_30".replace(" ", "_")
    print(f"{outdir}/{bname}_deg")
    if ("younger31" in adata.obs["age_range"].values) and ("older70" in adata.obs["age_range"].values) :
        print("y")
        a = len(adata[(adata.obs["subclass"] == ct1 ) & (adata.obs["age_range"] == "younger31")])
        b = len(adata[(adata.obs["subclass"] == ct1 ) & (adata.obs["age_range"] == "older70")])

        if ( a >= 50 ) and ( b >= 50 ) : 
            print("y1")
            adata1=adata[(adata.obs["subclass"]== ct1 )&(adata.obs["age_range"].isin(["younger31","older70"]))]

            sc.tl.rank_genes_groups(adata1, "age_range", method="wilcoxon", use_raw=False)

            deg1={ "score": list(pd.DataFrame(adata1.uns['rank_genes_groups']['scores'])["older70"]),
            'pvals': list(pd.DataFrame(adata1.uns['rank_genes_groups']['pvals'])["older70"]),
            'pvals_adj': list(pd.DataFrame(adata1.uns['rank_genes_groups']['pvals_adj'])["older70"]),
            'logfoldchanges': list(pd.DataFrame(adata1.uns['rank_genes_groups']['logfoldchanges'])["older70"])}

            data2=pd.DataFrame(data=deg1,index=pd.DataFrame(adata1.uns['rank_genes_groups']['names'])["older70"])
            data2.to_csv(f"{outdir}/{bname}_deg",sep="\t")

