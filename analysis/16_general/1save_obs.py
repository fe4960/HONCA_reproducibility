import scanpy as sc
import sys
h5ad=sys.argv[1]
adata=sc.read(f"{h5ad}.h5ad")

# check
print(adata.var.columns[adata.var.columns == "_index"])
print(adata.raw.var.columns[adata.raw.var.columns == "_index"])


if "_index" in adata.var.columns:
    adata.var = adata.var.drop(columns=["_index"])

if adata.raw is not None and "_index" in adata.raw.var.columns:
    adata.raw._var = adata.raw.var.drop(columns=["_index"])


adata.obs.to_csv(f"{h5ad}_simple.obs.gz")
adata.var.to_csv(f"{h5ad}_simple.var.gz")


if "counts" in adata.layers.keys():
    adata.X=adata.layers["counts"]
del adata.layers
del adata.raw

adata.write(f"{h5ad}_simple.h5ad")

