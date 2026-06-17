import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

go="enrichr_MSigDB_Hallmark_2020"

dir1="HCA_ON/data/20_hotspot/Oligodendrocyte_precursor_cell_permut_hs_min_gene_200_hvg_5000_age_m"

for i in range(1,6):
    file=f"{dir1}_{i}_bg/{go}.txt"
    df=pd.read_csv(file,header=0,sep="\t")
# Example DataFrame with a color column
#df = pd.DataFrame({
#    'Category': ['A', 'B', 'C', 'D'],
#    'Value': [10, 24, 36, 18],
#    'Group': ['X', 'Y', 'X', 'Z']  # group determines color
#})
    df["Log10P"]=-np.log10(df["P.value"])
    df=df.head(5).copy()
    df["Log10P"]=df["Log10P"].astype("float")
# Use hue or set palette manually
    sns.barplot(data=df, x='Log10P', y='Term', hue='Log10P', dodge=False,
            palette='viridis')

    plt.title(f"Module {i}")
    plt.xlabel("Log10P")
    plt.ylabel("GO")
    plt.tight_layout()
    plt.show()
    plt.savefig(f"{file}_{i}_top5.pdf")
