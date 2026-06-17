import scanpy as sc
c="Oligodendrocyte_precursor_cell"
adata=sc.read(f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{c}/clean/{c}_subclass.h5ad')

mk={"OPC":["PDGFRA","CSPG4","OLIG2"], "COP": ["GPR17","VCAN","NEU4","BMP4","FYN","SOX6","PTPRZ1","PDCD4"],"NFOL":["TCF7L2","CEMIP2"], "MFOL":["MBP","OPALIN","CTPS1"],"MOL":["MBP", "PLP1","MOG","GRM3","KLK6","PTGDS","ANXA5","HOPX"]}
sc.settings.figdir=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{c}/clean/figures'
sc.pl.dotplot(adata,mk,groupby="subclass",save=f"_{c}_subclass_ABC.png")

#mk={"generic":["OLIG1", "OLIG2", "SOX10", "MBP", "MOG"], "OPC":[ "PTPRZ1", "CA10","PDGFRA", "CSPG4"], "COP":["SOX6", "GPR17", "BMPER","PRICKLE1"]}
mk={"generic":["OLIG1", "OLIG2", "SOX10", "MBP","MOG"],"OPC":[ "PTPRZ1", "SOX6", "CA10", "PDGFRA"],"COP":[ "CSPG4", "GPR17", "BMPER", "PRICKLE1"],"Oligo1": ["SLC5A11", "RBFOX1", "RASGRF2", "RASGRF1", "ABCA6", "ACTN2", "PLEKHG1", "ACSBG1",  "LGALS1", "ELOVL2"], "Oligo2":[ "PLXDC2", "LINC01608",  "QDPR", "DYSF", "OPALIN", "LAMA2", "SVEP1"] , "Oligo3":[ "FRY", "SGCZ", "KCNK10", "RRAS2"], "Oligo4":[ "HSPH1", "ELOVL5", "HSPD1", "FOS", "DNAJA4"]}
sc.pl.dotplot(adata,mk,groupby="subclass",save=f"_{c}_subclass_cellReport_36001972_MS.png")

mk={"generic":["OLIG1", "OLIG2", "SOX10", "MBP", "MOG"], "OPC":[ "PTPRZ1", "CA10","PDGFRA", "CSPG4"], "COP":["SOX6", "GPR17", "BMPER","PRICKLE1"], "Oligo1":[ "SLC5A11", "RASGRF2", "ACTN2", "RASGRF1", "ABCA6", "PLEKHG1", "RBFOX1", "ACSBG1", "LGALS1", "ELOVL2"], "Oligo2":[ "QDPR", "OPALIN", "PLXDC2",  "LINC01608", "DYSF", "SVEP1", "LAMA2"]}

sc.pl.dotplot(adata,mk,groupby="subclass",save=f"_{c}_subclass_cellReport_36001972_AD.png")

mk={"OL lineage":["SOX10","OLIG1","OLIG2","PLP1","MYRF","MOG","MAG"], "OPC":["PTPRZ1","PDGFRA","VCAN"],"COP":["BMP4","GPR17","BCAS1"],"NFOL":["RRAS2","PROM1"],"MFOL":["THBS3","NHERF2"],"MOL1":"FOS", "MOL2":["KLK6","HOPX"], "MOL5/6":["PTGDS","IL33"],"DA1":["SERPINA3","C4B"],"DA2":["CDKN1A","TNFRSF12A"],"IFN":["IRF9","IFIT1"]}

sc.pl.dotplot(adata,mk,groupby="subclass",save=f"_{c}_subclass_cellReport_36001972_mouse.png")


c="Oligodendrocyte"

adata=sc.read(f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_new.h5ad")
sc.settings.figdir=f'/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/{c}/clean/figures'
sc.pl.dotplot(adata,mk,groupby="subclass",save=f"_{c}_subclass_new_ABC.png")

#PALM2AKAP2 is PALM2 and AKAP2 fusion but marker genes should be PALM2

mk={"generic":["OLIG1", "OLIG2", "SOX10", "MBP","MOG"],"OPC":[ "PTPRZ1", "SOX6", "CA10", "PDGFRA"],"COP":[ "CSPG4", "GPR17", "BMPER", "PRICKLE1"],"Oligo1": ["SLC5A11", "RBFOX1", "RASGRF2", "RASGRF1", "ABCA6", "ACTN2", "PLEKHG1", "ACSBG1",  "LGALS1", "ELOVL2"], "Oligo2":[ "PLXDC2", "LINC01608", "QDPR", "DYSF", "OPALIN", "LAMA2", "SVEP1"] , "Oligo3":[ "FRY", "SGCZ", "KCNK10"
    , "RRAS2"], "Oligo4":[ "HSPH1", "ELOVL5", "HSPD1", "FOS", "DNAJA4"]}

sc.pl.dotplot(adata,mk,groupby="subclass",save=f"_{c}_subclass_new_cellReport_36001972.png")


mk={"generic":["OLIG1", "OLIG2", "SOX10", "MBP", "MOG"], "OPC":[ "PTPRZ1", "CA10","PDGFRA", "CSPG4"], "COP":["SOX6", "GPR17", "BMPER","PRICKLE1"], "Oligo1":[ "SLC5A11", "RASGRF2", "ACTN2", "RASGRF1", "ABCA6", "PLEKHG1", "RBFOX1", "ACSBG1", "LGALS1", "ELOVL2"], "Oligo2":[ "QDPR", "OPALIN", "PLXDC2",  "LINC01608", "DYSF", "SVEP1", "LAMA2"]}

sc.pl.dotplot(adata,mk,groupby="subclass",save=f"_{c}_subclass_new_cellReport_36001972_AD.png")

mk={"OL lineage":["SOX10","OLIG1","OLIG2","PLP1","MYRF","MOG","MAG"], "OPC":["PTPRZ1","PDGFRA","VCAN"],"COP":["BMP4","GPR17","BCAS1"],"NFOL":["RRAS2","PROM1"],"MFOL":["THBS3","NHERF2"],"MOL1":"FOS", "MOL2":["KLK6","HOPX"], "MOL5/6":["PTGDS","IL33"],"DA1":["SERPINA3","C4B"],"DA2":["CDKN1A","TNFRSF12A"],"IFN":["IRF9","IFIT1"]}


mk={"OL lineage":["SOX10","OLIG1","OLIG2","PLP1","MYRF","MOG","MAG"], "OPC":["PTPRZ1","PDGFRA","VCAN"],"COP":["BMP4","GPR17","BCAS1"],"NFOL":["RRAS2","PROM1"],"MFOL":["THBS3","NHERF2"],"MOL1":"FOS", "MOL2":["KLK6","HOPX"], "MOL5/6":["PTGDS","IL33"],"DA1":["SERPINA3","C4B"],"DA2":["CDKN1A","TNFRSF12A"],"IFN":["IRF9","IFIT1"]}

sc.pl.dotplot(adata,mk,groupby="subclass",save=f"_{c}_subclass_cellReport_36001972_mouse.png")


