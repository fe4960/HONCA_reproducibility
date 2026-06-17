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
od=sys.argv[4]
label=sys.argv[5]
od1=sys.argv[2]
gn=sys.argv[6]
fn2=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/{od1}/ct_list_{label}"
#/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_trg500_enrichR/OLIGO_non_specific/enrichr_GO_Biological_Process_2023.txt
#fn2=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/ct_list_{ct}"
cell2=pd.read_csv(fn2,header=None)
cell=cell2[0].values.astype("str")
#mol=["MSigDB_Hallmark_2020"]
mol=["GO_Biological_Process_2023", "MSigDB_Hallmark_2020"]
#t=["up","down"]
#t=["down"]
t=["up"]
n=int(sys.argv[3])
for m in mol:
    cell1={}
    cell1["up"]={}
#    cell1["down"]={}
    go={}
    for c in cell:
        for t1 in t:
#            fn=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{ct}/clean/{od}_trg500_enrichR/{c}/enrichr_{m}.txt'
            fn=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{ct}/clean/{od}_trg{gn}_enrichR/{c}/enrichr_{m}.txt'

#            fn=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/{ct}_subclass_DEG/{ct}/{c}_DEG_GO_age_{t1}/enrichr_{m}.txt'
            if os.path.exists(fn):
                with open(fn,"r") as fn1:
                    next(fn1)
                    count=0
                    label=t1
                    if label not in go:
                        go[label]={}
                    for line in fn1:
                        info=line.split("\t")
                        item=info[0].split()
                        if m == "GO_Biological_Process_2023":
                            item1=" ".join(item[0:-1])
                        else:
                            item1=info[0]
                        if float(info[3]) < 0.1:

#                        if float(info[3]) < 0.05:
                            count+=1
                            if count <= n:
                                if item1 not in cell1[label]:
                                    cell1[label][item1]={}
                        if c not in go[label]:
                            go[label][c]={}
                        if item1 not in go[label][c]:
                            go[label][c][item1]={}
#                        go[label][c][item1]=info[3]
                        go[label][c][item1]=info[7]

    for l in go:
        fn=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{ct}/clean/{od}_trg{gn}_enrichR/{ct}_GO_top{n}_{m}_{l}'  # msigdbr top 10
        fn2=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{ct}/clean/{od}_trg{gn}_enrichR/{ct}_GO_top{n}_{m}_{l}_column'  # msigdbr top 10

        #        fn=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/{ct}_subclass_DEG/{ct}/{ct}_GO_top{n}_{m}_{dem}_{l}'  # msigdbr top 10
#        fn2=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/{ct}_subclass_DEG/{ct}/{ct}_GO_top{n}_{m}_{dem}_{l}_column'  # msigdbr top 10
        fn3=open(fn2, "w")        
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
#                        line+="\t"+str(-math.log10(float(go[l][key][c])))
                   #     if -math.log10(float(go[l][key][c])) > 5:
                   #         val=str(5)
                   #     else:
                   #         val=str(-math.log10(float(go[l][key][c])))
                        val=go[l][key][c]
                        fn3.write(f"{key}\t{c}\t{val}\n")
                        line+="\t"+val
                    else:
                        line+="\t"+"0" #"1"
                        fn3.write(f"{key}\t{c}\t0\n")                        
                fn1.write(f'{line}\n')
        fn3.close()
        cor=pd.read_csv(fn,sep="\t",header=0,index_col=0)
        df=cor.T
#        plt.figure(figsize=(1,4))
#        if l == "down": 
        plt.figure(figsize=(4,2))
#####        plt.title(f"{ct} GO_Biological_Process_2023 {l}")
        plt.title(f"{ct} MSigDB_Hallmark_2020 {l}")

        # Perform hierarchical clustering on columns
        if len(cor.columns) > 1:
            col_linkage = linkage(df.T, method='average')
            #g=sns.clustermap(df, col_linkage=col_linkage, cmap='YlGnBu')
            g=sns.clustermap(df, col_linkage=col_linkage, cmap='viridis') #BuPu')
            plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha="right", fontsize=20)
            plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=20) #, fontsize=7) 
            g.fig.subplots_adjust(bottom=0.2)
#            g.fig.subplots_adjust(right=0.6)
        else:
            g=sns.heatmap(df,  cmap='YlGnBu')
            g.set_xticklabels(g.get_xticklabels(), rotation=45)
# Plot heatmap with clustered columns
        if l == "up":
            if len(cor.columns) > 1:
                g=sns.clustermap(df, col_linkage=col_linkage, cmap='YlOrRd')
                plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha="right", fontsize=20)
                plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=20) #, fontsize=7) 
                g.fig.subplots_adjust(bottom=0.2)
 #               g.fig.subplots_adjust(right=0.6)
            else:
                g=sns.heatmap(df, cmap='YlOrRd')
                g.set_xticklabels(g.get_xticklabels(), rotation=45)
        #inferno') #, standard_scale=1)
	# using the upper triangle matrix as mask
#        g=sns.heatmap(cor2, annot=False,  annot_kws=annot_kws )
#        g.set_xticklabels(g.get_xticklabels(), rotation=30)
#        plt.savefig(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/{od}/{ct}_subclass_DEG/{ct}/{ct}_GO_top{n}_{m}_{dem}_{l}_limit10.pdf',dpi=500)
        plt.savefig(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{ct}/clean/{od}_trg{gn}_enrichR/{ct}_GO_top{n}_{m}_{l}.pdf',dpi=500)
        plt.close("all") 



						
