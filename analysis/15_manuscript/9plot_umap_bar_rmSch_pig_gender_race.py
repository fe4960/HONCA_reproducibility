#import scanpy as sc
import pandas as pd
import numpy as np
#from matplotlib import cm
#import matplotlib.pyplot as plt
#from matplotlib import rcParams
#sc.settings.figdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean/figures/"
########
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/clean"
bname="major_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_rmUnk"

bname=f"{bname}_clean"

obs=pd.read_csv(f"{outdir}/{bname}.obs.gz",sep="\t", header=0, index_col=0)

obs[["nCount_RNA","majorclass1","gender"]].groupby(["gender","majorclass1"]).count().to_csv(f"{outdir}/{bname}_gender.txt",sep="\t")

obs[["nCount_RNA","majorclass1","race"]].groupby(["race","majorclass1"]).count().to_csv(f"{outdir}/{bname}_race.txt",sep="\t")




