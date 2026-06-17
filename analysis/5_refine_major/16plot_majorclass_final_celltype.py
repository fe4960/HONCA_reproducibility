import pandas as pd
import scanpy as sc
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

rcParams["figure.figsize"] = (5,5)


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_all_normcount.h5ad")

dir="/dfs3b/ruic20_lab/junw42"
#f1="hvg2000_nonsb_epochnone_seurat_v3_rs"
#bname=f"_major_{f1}_1_clean_ON_ONH_new"
bname=f"major_ON_ONH_clean_new"
#adata=sc.read(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/major_ON_ONH_new_final.h5ad")
df=pd.DataFrame(adata.obs["celltype"].value_counts())
#df["cell type"]=df.index
#df["cell number"]=df["celltype"]
#df["cell type"]=df["cell type"].replace("Oligodendrocyte_precursor_cell","OPC")
#df["major class"]=df["major class"].replace("Unknown?","Immune_neuron_unk")



#df = df.sort_values(by='cell number', ascending=False)
#fig=plt.figure(figsize=(10, 8))
#sns.barplot(x='cell number', y="cell type", data=df, order=df["cell type"],palette="tab20")
#plt.show()
#fig.savefig(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/figures/{bname}_celltype.png", dpi=300)

df.to_csv(f"{dir}/HCA_ON/data/5_refine_major/scvi/major/figures/{bname}_celltype.csv", sep = "\t" )

#adata=sc.read_h5ad("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_all_normcount.h5ad")
