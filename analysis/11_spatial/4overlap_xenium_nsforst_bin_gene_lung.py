import pandas as pd
import ast
import json

# Define the path to your JSON file
#file_path = "/dfs3b/ruic20_lab/rawdata/Xenium/20241018__Human_Eye_MDF_FFPE_lung_probe/MDF/output-XETG00221__0029333__24-1065_OS_PP__20241018__000131/gene_panel.json"

# Open and read the JSON file
#with open(file_path, "r") as file:
#    data = json.load(file)


df=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/nsforest/nsforest_celltypecelltype_results.csv",header=0,index_col=0)
#xiem=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/Xenium_hSkin_v1_metadata_gene",header=None)
xiem=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/Xenium_hLung_v1_metadata.csv",header=0)
file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/Xenium_hLung_v1_metadata_ns_bin_gene"
out=open(file,"w")
for i in df.index:
    df1=ast.literal_eval(df.loc[i,"binary_genes"])
    for g in df1:
        if g in list(xiem["Gene"]):
            out.write(f"{i}\t{g}\n")
out.close()


file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/Xenium_hLung_v1_metadata_ns_marker_gene"
out=open(file,"w")
for i in df.index:
    df1=ast.literal_eval(df.loc[i,"NSForest_markers"])
    for g in df1:
        if g in list(xiem["Gene"]):
            out.write(f"{i}\t{g}\n")
out.close()

