import pandas as pd
import numpy as np

file1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/ONONH_new_final.obs.gz"

obs=pd.read_csv(file1,sep="\t",header=0,index_col=0)

df=obs.loc[obs["majorclass"].isin(["Macrophage","Microglia"]),["gender","age_year","majorclass","tissue","race","sampleid"]].groupby(["sampleid","majorclass"]).count()
df1=obs.loc[obs["majorclass"].isin(["Macrophage","Microglia"]),["gender","age_year","majorclass","tissue","race","sampleid"]].groupby(["sampleid"]).count()
dfp = df / df1
dfp.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/Micro_Macro_prop_age")
dfp=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/Micro_Macro_prop_age",header=0)
obs1=obs[["gender","age_year","tissue","race","sampleid"]].drop_duplicates()
obs1=obs1.set_index("sampleid")
dfp["age"]=list(obs1.loc[dfp["sampleid"],"age_year"])
del dfp["tissue"]
del dfp["race"]
del dfp["gender"]
dfp["tissue"]=list(obs1.loc[dfp["sampleid"],"tissue"])
dfp["race"]=list(obs1.loc[dfp["sampleid"],"race"])
dfp["gender"]=list(obs1.loc[dfp["sampleid"],"gender"])
dfp.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/Micro_Macro_prop_age_anno")
#conf_mat.to_csv(f'{dir}/HCA_ON/data/5_refine_major/scvi/major/{bname}_source_conf_mat_table.csv',header=True, index=True)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#from sklearn.linear_model import LinearRegression

# Generate random data
#np.random.seed(42)
#x = np.random.rand(50) * 10  # Random values for x
#y = 2.5 * x + np.random.randn(50) * 5  # Linear relation with some noise

# Create a DataFrame
#df = pd.DataFrame({'x': x, 'y': y})
df_mic=dfp.loc[(dfp["majorclass"]=="Macrophage") and (dfp["tissue"]=="ON") , ]
# Scatter plot with regression line using Seaborn
plt.figure(figsize=(8, 5))
sns.regplot(x='age', y='age_year', data=df_mic, scatter_kws={'color': 'blue'}, line_kws={'color': 'red'})
plt.xlabel("Age")
plt.ylabel("# of macro / (# of micro and macro)")
plt.title("Macrophage proportion along age")
plt.savefig("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/Micro_Macro_prop_age.png")
