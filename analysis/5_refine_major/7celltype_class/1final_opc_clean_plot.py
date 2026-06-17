import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches

rcParams["figure.figsize"] = (5,5)
celltype="Oligodendrocyte_precursor_cell"
indir=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{celltype}/"
outdir=f"{indir}/clean"
sc.settings.figdir = f"{outdir}/figures"

import scanpy as sc
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_new_clean.h5ad")
mk={"OL lineage":["SOX10","OLIG1","OLIG2","PLP1","MYRF","MOG","MAG"],"OPC":["PTPRZ1","PDGFRA","VCAN"],"COP":["BMP4","GPR17","BCAS1"],"NFOL":["RRAS2","PROM1"],"MFOL":["THBS3"],"MOL1":["FOS"],"MOL2":["KLK6","HOPX"],"MOL5/6":["PTGDS","IL33"],"DA1":["SERPINA3","C4B"],"DA2":["CDKN1A","TNFRSF12A"],"IFN":["IRF9","IFIT1"]}

sc.pl.dotplot(adata,mk,groupby="subclass",save="Oligodendrocyte_subclass_new_clean_opc_oligo_mk.png")

adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass.h5ad")
sc.pl.dotplot(adata,mk,groupby="subclass",save="Oligodendrocyte_precursor_cell_subclass_new_clean_opc_oligo_mk.png")

