import scanpy as sc
import scvelo as scv
import os
import pandas as pd
import anndata as ad
import sys

h5ad=sys.argv[1]
indir=sys.argv[2]
#adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/oligo_opc.h5ad")

adata=sc.read(f'{indir}/{h5ad}.h5ad')

def find_file(i):
    if i.find("GSM") > -1: 
        dir1=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE236566/cellranger/{i}/{i}/outs/'
    else:
        dir1=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pri/{i}/'
    files = os.listdir(dir1)
    for file1 in files:
        if file1.find("loom")>0:
            return(f'{dir1}/{file1}')

main="/dfs3b/ruic20_lab/junw42/"

samplelist=f"{main}/HCA_ON/data/0_sample_list/RNA_sample_list_meta.gz"
#samplelist=f"{main}/HCA_ON/data/0_sample_list/RNA_sample_list_meta_header"

sl=pd.read_csv(samplelist,sep="\t",header=0,index_col=9)

loom_list=[]
n=0
#loom_full=ad.AnnData()
for i in adata.obs["sampleid_legacy"].value_counts().index:
    n=n+1
    file1=find_file(i)
    adata_sp=sc.read(file1)
    j=sl.loc[i,sl.columns[0]]
    adata_sp.obs.index = [j + "_" + x.split(":", 1)[1][:-1] + "-1" for x in adata_sp.obs.index]
    idx=adata_sp.obs.index.intersection(adata.obs.index)
    adata_sp=adata_sp[idx].copy()

    adata_sp.var_names_make_unique()
    loom_list.append(adata_sp)
#    if n == 2:
#        break

loom_full=ad.concat(loom_list)
loom_full.write(f"{indir}/{h5ad}_sp_loom.h5ad")
#loom_full.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/oligo_opc_sp_loom.h5ad")
