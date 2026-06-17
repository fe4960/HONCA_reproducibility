import sys
import scanpy as sc
import anndata
import pandas as pd
import os
import numpy as np
from pathlib import Path


dir1=sys.argv[1]
path=sys.argv[2]
label=sys.argv[5]
ct_list=sys.argv[6]
print(ct_list)
#dir1="/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/genexp_sample_raw_atlas"
Path(dir1).mkdir(parents=True,exist_ok=True)

#path="/storage/chen/data_share_folder/jinli/gptprj/Rui_Bo_datashare/HRCAv2_snRNA.h5ad"
#path="/storage/chen/data_share_folder/jinli/gptprj/Rui_Bo_datashare/Hahn/HRCAv2.h5ad"
#path="/storage/chentemp/wangj/human_ret_anc/data/HRCA/HRCA.h5ad"

adata=sc.read_h5ad(path)

#adata.obs["donor"]=adata.obs["donor"].replace("Chen_17D013","Chen_D017_13")

#adata.obs["donor"]=adata.obs["donor"].replace("Chen_D017_13","Chen_17D013")
#adata.obs["donor"]=adata.obs["donor"].replace("Chen_19_D016","Chen_19D016")


def rename_gene(adata):
	for i in range(0, adata.var.features.shape[0]):
		if adata.var.features[i] in gene_dict:
			adata.var.features[i] = gene_dict[adata.var.features[i]]
	return adata


def grouped_obs_cpm(adata, group_key, layer=None, gene_symbols=None):
	if layer is not None:
		getX= lambda x: x.layers[layer]
	else:
                getX = lambda x: x.X
#		getX = lambda x: x.raw.X
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
#		factor=1000000.0
#		out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))/total_count*factor
		out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))
	return(out)



celltype=adata.obs[label].value_counts().index

if ct_list != "none":
    tmp=pd.read_csv(ct_list,header=None)
    celltype=tmp[0].values.astype(str)
#celltype=["AC","BC","Cone","HC","MG","RGC","Rod","Astrocyte","Microglia"] 
cellnum={}

list1=sys.argv[3] #"/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/snRNA_all_sample_celltype_major_num_all_atlas5"

out=open(list1,"w")



for donor in adata.obs["sampleid"].value_counts().keys() :
		adata_sub = adata[adata.obs["sampleid"]==donor]
#		adata_sub.layers["raw"] = adata_sub.X
		genexp=grouped_obs_cpm(adata=adata_sub, group_key=label) #, layer="raw")
		fileout=f'{dir1}/{donor}.txt.gz'
		genexp.to_csv(fileout,sep="\t")
		for cell in celltype :
			if cell in adata_sub.obs[label].value_counts():
				if donor not in cellnum:
					cellnum[donor]={}
				cellnum[donor][cell] = adata_sub.obs[label].value_counts()[cell]
				num_tmp=adata_sub.obs[label].value_counts()[cell]
				out.write(f'{donor}\t{cell}\t{num_tmp}\n')
			else:
				if donor not in cellnum:
					cellnum[donor]={}
				cellnum[donor][cell] = 0
				num_tmp=0
				out.write(f'{donor}\t{cell}\t{num_tmp}\n')


out.close()
df=pd.DataFrame.from_dict(cellnum)

list2=sys.argv[4] #"/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/snRNA_all_sample_celltype_major_num_atlas"

df.to_csv(f'{list2}')

