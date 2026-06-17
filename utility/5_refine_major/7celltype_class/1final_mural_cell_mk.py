import scanpy as sc
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Mural_cell/clean/Mural_cell_subclass.h5ad")
mk={"aSMC":["ACTA2","MYH11","CNN1","TAGLN","ELN","S1PR1","KCNK3"],"vSMC":["DES","MRC1","CD74","ACTA2","VIM","MCAM","EDF1","PDGFRB"],"pericytes":["GRM8","PDGFRB"]}
sc.settings.figdir = f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Mural_cell/clean/figures"
sc.pl.dotplot(adata,mk,groupby="subclass",save="Mural_cell_subclass_mk_35158371.png")

mk={"aSMC":["ACTA2","MYH11"],"vSMC":["MRC1","CD74"],"pericytes":["GRM8","PDGFRB"]}
sc.pl.dotplot(adata,mk,groupby="subclass",save="Mural_cell_subclass_mk_35158371_real_short.png")

mk={"contractile":["MYH11","TAGLN","ACTA2","CNN1","SMTN","LMOD1","KCNMB1","KCNA5","SYNM","MYOCD"],"Synthetic":["COL1A1","MYH10","SPP1","KCNA3","KLF4","LGALS3"]}
sc.pl.dotplot(adata,mk,groupby="subclass",save="Mural_cell_subclass_mk_Chakraborty.png")
