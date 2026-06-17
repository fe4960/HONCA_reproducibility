import scanpy as sc
import sys
import pandas as pd
indir=sys.argv[1]
mk=sys.argv[2]
outdir=sys.argv[3]
bname=sys.argv[4]

#adata=sc.read(f'{outdir}/{bname}_scvi_trg.h5ad')
adata=sc.read(f'{outdir}/{bname}.h5ad')

sc.settings.figdir = f"{outdir}/figures"

#####lis=pd.read_csv(f'{indir}/{mk}',header=None)
#####dt=lis[0].values #.set_index(0).to_dict()[1]

#######sc.tl.dendrogram(adata_query, groupby='leiden', use_rep="X_scVI", use_raw=False, n_pcs=30)
########sc.pl.dendrogram(adata_query, groupby='leiden', save=f"{bname}_dendrogram.png")

########sc.pl.dotplot(adata_query, dt, groupby='leiden', dendrogram=True, save=f"{bname}_marker_sanes_ct.png", use_raw=False)

#########sc.pl.violin(adata,keys=["pANN"], groupby="leiden", rotation=75, save=f"{bname}_pANN.png")
#adata.obs[["nCount_RNA","leiden","sampleid"]].groupby(["leiden","sampleid"]).count().to_csv(f"{outdir}/{bname}_sample_number.txt",sep="\t")

df=pd.read_csv(f"{outdir}/{bname}_sample_number.txt",header=0)
total=df.groupby(["subclass"])["nCount_RNA"].transform("sum")
df["nCount_RNA"]=df["nCount_RNA"]/total
for i in pd.unique(df["sampleid"]):
    plt.barh(df[df["sampleid"]==i]["subclass"], df[df["sampleid"]==i]["nCount_RNA"],  left=df[df["source"]!="Chen"]["nCount_RNA"], color='lightblue' )
plt.legend(["Sanes", "Chen"])
plt.xlabel("Cell proportion")
plt.savefig(f"{outdir}/figures/{bname}_source.png")
