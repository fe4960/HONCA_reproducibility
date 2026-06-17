import scanpy as sc
from matplotlib import rcParams
import sys
rcParams["figure.figsize"] = (5,5)
outdir=sys.argv[1]
bname=sys.argv[2]

sc.settings.figdir = f"{outdir}/figures"

#adata=sc.read(f'{outdir}/{bname}_scvi_trg.h5ad')
adata=sc.read(f'{outdir}/{bname}.h5ad')
mk=["BMP4", "BMP5", "BMP7", "PIEZO1", "PIEZO2"]
sc.pl.umap(adata, color=mk, ncols=3, save=f'_{bname}_mechano.png')
#sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,palette="tab20", save=f'_{bname}_subclass.png')

