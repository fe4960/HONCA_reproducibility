import sys
import scanpy as sc
import anndata
import pandas as pd
import os
import numpy as np
from pathlib import Path
ct=sys.argv[1]
fn=sys.argv[2]
od=sys.argv[3]
lab=sys.argv[4]
main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/"
#dir1="/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/genexp_sample_raw_atlas"
#dir1=f"{main}/7_DEG/genexp_sample_raw_atlas"
#dir1=f"{main}/7_DEG/genexp_donor_raw_atlas"
#dir1=f"{main}/7_DEG/genexp_donor_subclass_final/Astrocyte/"
dir1=f"{main}/7_DEG/{od}/{ct}/"

Path(dir1).mkdir(parents=True,exist_ok=True)

#path=f"{main}/5_refine_major/scvi/major/ON_ONH_clean_new.h5ad"
#path=f"{main}/5_refine_major/scvi/major/major_ON_ONH_new_final.h5ad"
#path=f"{main}/5_refine_major/scvi/lattice/ONONH_Astrocyte_rawcount.h5ad"
#path=f"{main}/5_refine_major/scvi/Astrocyte/clean/ONONH_Astrocyte_rawcount_2subclass.h5ad"
#path=f"{main}/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean.h5ad"
path=f"{main}/5_refine_major/scvi/{ct}/clean/{fn}.h5ad"

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
#		factor=1000000.0
#		out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))/total_count*factor
		out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))
	return(out)

#celltype=["Oligodendrocyte","Astrocyte","Fibroblast","Rod","Microglia","Endothelial_cell","BC","Oligodendrocyte_precursor_cell","Mural_cell","MG","AC","RPE","Macrophage","Melanocyte","HC","Cone","RGC","Schwann_cell"] 
celltype=list(adata.obs[lab].unique()) #["Astro_ret","Astro_ON"]
cellnum={}

#list1=f"{main}/7_DEG/genexp_donor_subclass/snRNA_donor_num_cell_Astrocyte" #snRNA_donor_majorclass_num_atlas_cell"
list1=f"{main}/7_DEG/{od}/snRNA_donor_num_cell_{ct}" #snRNA_donor_majorclass_num_atlas_cell"

out=open(list1,"w")



for donor in adata.obs["donor"].value_counts().keys() :
		adata_sub = adata[adata.obs["donor"]==donor]
		adata_sub.layers["raw"] = adata_sub.X
		genexp=grouped_obs_cpm(adata=adata_sub, group_key=lab, layer="counts")
#		genexp=grouped_obs_cpm(adata=adata_sub, group_key="subclass", layer="counts")
#		genexp=grouped_obs_cpm(adata=adata_sub, group_key="subclass", layer="raw")
		fileout=f'{dir1}/{donor}.txt.gz'
		genexp.to_csv(fileout,sep="\t")
		for cell in celltype :
			if cell in adata_sub.obs[lab].value_counts():
    #			if cell in adata_sub.obs["subclass"].value_counts():
				if donor not in cellnum:
					cellnum[donor]={}
				cellnum[donor][cell] = adata_sub.obs[lab].value_counts()[cell]
				num_tmp=adata_sub.obs[lab].value_counts()[cell]
				out.write(f'{donor}\t{cell}\t{num_tmp}\n')
			else:
				if donor not in cellnum:
					cellnum[donor]={}
				cellnum[donor][cell] = 0
				num_tmp=0
				out.write(f'{donor}\t{cell}\t{num_tmp}\n')


out.close()
df=pd.DataFrame.from_dict(cellnum)

#list=f"{main}/7_DEG/genexp_donor_subclass/snRNA_donor_num_Astrocyte" #snRNA_donor_majorclass_num_atlas"
list=f"{main}/7_DEG/{od}/snRNA_donor_num_{ct}" #snRNA_donor_majorclass_num_atlas"

df.to_csv(f'{list}')

