import pandas as pd

data=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta",sep=",",index_col=7, header=0)

path=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/ON_chen_crpath.txt",sep="\t",index_col=9,header=0)

tmp=["crpath", "crpath_legacy",  "fastqpath"]

idx=data.index.intersection(path.index)

data1=data[data.index.isin(idx)]

for t in tmp:
    data1.loc[idx,t]=path.loc[idx,t]

data1.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/onoh_pri_meta_datapath",sep="\t")


path1=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/ONONH_meta_pub.txt",sep="\t",index_col=0,header=0)

data

idx=data.index.intersection(path1.index)

data2=data[data.index.isin(idx)]

tmp=["crpath",   "fastqpath"]

for t in tmp:
    data2.loc[idx,t]=path1.loc[idx,t]

data2.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/onoh_pub_meta_datapath",sep="\t")

