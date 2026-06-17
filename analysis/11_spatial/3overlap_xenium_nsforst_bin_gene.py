import pandas as pd
import ast
df=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/nsforest/nsforest_celltypecelltype_results.csv",header=0,index_col=0)
xiem=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/Xenium_hSkin_v1_metadata_gene",header=None)
file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/Xenium_hSkin_v1_metadata_ns_bin_gene"
out=open(file,"w")
for i in df.index:
    df1=ast.literal_eval(df.loc[i,"binary_genes"])
    for g in df1:
        if g in list(xiem[0]):
            out.write(f"{i}\t{g}\n")
out.close()

