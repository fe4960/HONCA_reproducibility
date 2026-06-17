import scanpy as sc
import sys
import pandas as pd
indir=sys.argv[1]
mk=sys.argv[2]
outdir=sys.argv[3]
bname=sys.argv[4]

adata=sc.read(f'{outdir}/{bname}_scvi_trg.h5ad')
sc.settings.figdir = f"{outdir}/figures"

#####lis=pd.read_csv(f'{indir}/{mk}',header=None)
#####dt=lis[0].values #.set_index(0).to_dict()[1]

#######sc.tl.dendrogram(adata_query, groupby='leiden', use_rep="X_scVI", use_raw=False, n_pcs=30)
########sc.pl.dendrogram(adata_query, groupby='leiden', save=f"{bname}_dendrogram.png")

########sc.pl.dotplot(adata_query, dt, groupby='leiden', dendrogram=True, save=f"{bname}_marker_sanes_ct.png", use_raw=False)

######sc.pl.violin(adata,keys=["pANN"], groupby="leiden", rotation=75, save=f"{bname}_pANN.png")

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)
sc.pl.rank_genes_groups(adata, n_genes=20,sharey=False, save=f'_{bname}_TRG20.png')

sc.pl.rank_genes_groups_dotplot(adata, n_genes=10,  save=f'{bname}_TRG10.png', dendrogram=False)
