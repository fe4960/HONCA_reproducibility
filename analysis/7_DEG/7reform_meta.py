import pandas as pd
meta=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta",header=0,index_col=0)
gsm=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_race",header=None,index_col=0,sep="\t")
#dt=gsm.to_dict()[0]

meta.loc[gsm.index,"race"]=gsm[1]

meta.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race")
