import sys
import scanpy as sc
import anndata
import pandas as pd
import os
import numpy as np
from pathlib import Path
#ct=sys.argv[1]
#fn=sys.argv[2]
#od=sys.argv[3]
#lab=sys.argv[4]
#sampleid=sys.argv[5]
main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/"

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final"

dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/"

#dir1=f"{main}/7_DEG/{od}/{ct}/"

#Path(dir1).mkdir(parents=True,exist_ok=True)

#path=f"{main}/5_refine_major/scvi/{ct}/clean/{fn}.h5ad"
path=f"{dir}/ONONH_all_raw_normcount_only.h5ad"
adata=sc.read_h5ad(path)

adata=adata[adata.obs["celltype"]!="Unknown?"]

def rename_gene(adata):
	for i in range(0, adata.var.features.shape[0]):
		if adata.var.features[i] in gene_dict:
			adata.var.features[i] = gene_dict[adata.var.features[i]]
	return adata


def grouped_obs_cpm(adata, group_key, layer=None, gene_symbols=None):
	if layer is not None:
		getX= lambda x: x.layers[layer]
	else:
		getX = lambda x: x.raw.X
	if gene_symbols is not None:
		new_idx = adata.var[idx]
	else:
		new_idx = adata.var_names
	grouped = adata.obs.groupby(group_key)
	out = pd.DataFrame(
		np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
		columns=list(grouped.groups.keys()),
		index=adata.var_names
	)
	for group, idx in grouped.indices.items():
		X = getX(adata[idx])
		total_count=X.sum(dtype=np.float64)
		factor=1000000.0
		out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))/total_count*factor
#		out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))
	return(out)

#celltype=list(adata.obs["celltype"].unique()) #["Astro_ret","Astro_ON"]



lab="celltype"
genexp=grouped_obs_cpm(adata=adata, group_key=lab, layer="counts")
fileout=f'{dir1}/gene_exp_celltype_cpm.txt.gz'

#fileout=f'{dir1}/gene_exp_celltype.txt.gz'
genexp.to_csv(fileout,sep="\t")

lab="majorclass1"
genexp=grouped_obs_cpm(adata=adata, group_key=lab, layer="counts")
fileout=f'{dir1}/gene_exp_majorclass_cpm.txt.gz'

#fileout=f'{dir1}/gene_exp_majorclass.txt.gz'
genexp.to_csv(fileout,sep="\t")


