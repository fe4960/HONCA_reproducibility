import scanpy as sc
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass.h5ad")
mk={"Dural":["CXCL12","FXYD5","MGP","OGN","TGFBI"]}
sc.settings.figdir = f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/figures"
sc.pl.dotplot(adata,mk,groupby="subclass", save="Fibroblast_subclass_dura.png")
mk={"dura mater":["MGP", "GJA1", "FXYD5","COL18A1"],"arachnoid mater":["CLDN11", "TBX18", "TAGLN"],"pia mater":["LAMA2","S100A6","NGFR"]}
sc.pl.dotplot(adata,mk,groupby="subclass", save="Fibroblast_subclass_dura_detail.png")

mk={"dura mater":["FOXP1", "MGP", "FXYD5"],"arachnoid mater":["NNAT", "ALDH1A2", "PTGDS"],"pia mater":["LAMA1","CXCL12","COL15A1"],"Arachnoid barrier":["CDH1","CLDN11","TJP1"]}
sc.pl.dotplot(adata,mk,groupby="subclass", save="Fibroblast_subclass_dura_detail_wang_elife.png")

sc.settings.figdir = f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/figures"
mk={"Fb-d1": "CDH18", "Fb-d3":"SLC47A1","Fb-d2":"EPHA3","fb-a":"IGSF8","fb-p":"TMEM132C"}
sc.pl.dotplot(adata,mk,groupby="subclass", save="Fibroblast_subclass_dura_detail_wang_elife_top1.png")


import pandas as pd
adata2=adata
obs=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_res_0.6_clean_sb_res_0.6.obs.gz",header=0,index_col=0)
adata2.obs["leiden1"]=obs.loc[adata2.obs.index,"leiden"].astype("str")
mk={"fb x":["EBF2", "SHISA6"],"sclera":["FBLN1"],"dura mater":["FOXP1", "MGP", "FXYD5","SLC47A1"],"arachnoid mater":["NNAT", "ALDH1A2", "PTGDS"],"pia mater":["LAMA1","CXCL12","COL15A1"],"Arachnoid barrier":["CDH1","CLDN11","TJP1"],"Arachnoid barrier":["CDH1","SEMA3D","COL25A1","CXADR","SORBS1","WNT16"]}
sc.pl.dotplot(adata2,mk,groupby="leiden1", save="Fibroblast_subclass_dura_detail_wang_elife_sb_res_0.6.png")
