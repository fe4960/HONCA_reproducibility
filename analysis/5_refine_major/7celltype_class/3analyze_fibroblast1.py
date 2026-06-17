import pandas as pd
data=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_res_0.2_clean_sd7_res_0.2.obs.gz")
data["celltype"]=data["celltype"].astype("str")
df=data.groupby(["leiden","celltype"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=0).values[np.newaxis,:]
conf_mat.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_res_0.2_clean_sd7_res_0.2_conf_mat.csv")

data=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_res_0.3_clean_sd7_res_0.3.obs.gz")
df=data.groupby(["leiden","celltype"]).size().unstack(fill_value=0)
conf_mat = df / df.sum(axis=0).values[np.newaxis,:]
conf_mat.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_res_0.3_clean_sd7_res_0.3_conf_mat.csv")
