import pandas as pd

meta=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/Table_HCA_donorinfo.txt",index_col=0,header=0,sep="\t")
time=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/interval_to_time.txt_simple_format",index_col=0,header=None,sep="\t")
df=time.to_dict()[1]
meta["pmit"]=meta["pmit"].replace(df)
meta["pmit"]=meta["pmit"].replace("Unknown","NotApplicable")
meta["globe"]=meta["globe"].replace("Unknown","NotApplicable")
meta.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/Table_HCA_donorinfo.txt_format",sep="\t")

