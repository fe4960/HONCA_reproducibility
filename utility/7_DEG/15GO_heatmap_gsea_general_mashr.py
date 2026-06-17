import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
import math
#cell=["Rod","Cone","BC","RGC","MG","AC","HC"]
annot_kws={'fontsize': 11 }
#mol=["KEGG","msigdbr","GO_BP"]
ct=sys.argv[1]
fn2=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/ct_list_{ct}"
cell2=pd.read_csv(fn2,header=None)
cell=cell2[0].values.astype("str")
mol=["MSigDB_Hallmark_2020"]
dem=sys.argv[2]
t=["up","down"]
n=int(sys.argv[3])
for m in mol:
    cell1={}
    cell1["up"]={}
    cell1["down"]={}
    go={}
    for c in cell:
        for t1 in t:
            fn=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/{ct}_subclass_DEG/{ct}/{c}_DEG_GO_age_{t1}_mashr/enrichr_{m}.txt'
            if os.path.exists(fn):
                with open(fn,"r") as fn1:
                    next(fn1)
                    count=0
                    label=t1
                    if label not in go:
                        go[label]={}
                    for line in fn1:
                        info=line.split("\t")
                        count+=1
                        if count <= n:
                            if info[0] not in cell1[label]:
                                cell1[label][info[0]]={}
                        if c not in go[label]:
                            go[label][c]={}
                        if info[0] not in go[label][c]:
                            go[label][c][info[0]]={}
                        go[label][c][info[0]]=info[3]

    for l in go:
        fn=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/{ct}_subclass_DEG/{ct}/{ct}_GO_top{n}_{m}_{dem}_{l}_mashr'  # msigdbr top 10
        with open(fn, "w") as fn1:
            line=[]
            for c in cell1[l]:
                line.append(str(c))
		
            line1="\t".join(line)
            fn1.write(f'{line1}\n')

            for key in go[l]:
                line=str(key)
                for c in cell1[l]:
                    if c in go[l][key]:
                        line+="\t"+str(-math.log10(float(go[l][key][c])))
                    else:
                        line+="\t"+"0" #"1"
                fn1.write(f'{line}\n')
        cor=pd.read_csv(fn,sep="\t",header=0,index_col=0)
        df=cor.T
        plt.figure(figsize=(9,9))
        plt.title(f"{ct} MSigDB_Hallmark_2020 {dem} {l}")
        # Perform hierarchical clustering on columns
        col_linkage = linkage(df.T, method='average')

# Plot heatmap with clustered columns
        g=sns.clustermap(df, col_linkage=col_linkage, cmap='inferno') #, standard_scale=1)
	# using the upper triangle matrix as mask
#        g=sns.heatmap(cor2, annot=False,  annot_kws=annot_kws )
#        g.set_xticklabels(g.get_xticklabels(), rotation=30)
        plt.savefig(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/{ct}_subclass_DEG/{ct}/{ct}_GO_top{n}_{m}_{dem}_{l}_mashr.pdf',dpi=500)
        plt.close("all") 



						
