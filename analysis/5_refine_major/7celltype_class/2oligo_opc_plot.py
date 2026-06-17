import scanpy as sc
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams["font.sans-serif"] = ["Nimbus Sans"]
rcParams['font.size'] = 15

celltype="Oligodendrocyte"
indir="/dfs3b/ruic20_lab/junw42"
adata=sc.read(f"{indir}/HCA_ON/data/5_refine_major/scvi/{celltype}/clean/oligo_opc_seurat_new_scvi_cycling.h5ad")
adata=adata[(adata.obs["sampleid"]!="MMD_23_17976_ONH_RNA")&(adata.obs["sampleid"]!="MMD_23_20181_ONH_RNA")]

sc.settings.figdir=f"{indir}/HCA_ON/data/5_refine_major/scvi/{celltype}/clean/"

mk=["PIEZO2", "PIEZO1", "RAPGEF3", "KCNK1", "KCNK2", "KCNK10", "TRPV1", "TRPC1", "TRPC3", "PKD2", "TRPM7"]

dt={"COP_LAMA2+":"NFOL_UTRN+"}
adata.obs["subclass"]=adata.obs["subclass"].replace(dt)

sc.pl.dotplot(adata,mk, groupby="subclass", categories_order=["OPC","OPC_GLIS3+","OPC_Cycling", "COP_NFOL","NFOL_UTRN+", "MFOL","OLIGO2_LRRC7+","OLIGO1_SVEP1+","OLIGO1_RBFOX1+","OLIGO1"], save="_opc_oligo_mechano.png")



adata=sc.read(f"{indir}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat.h5ad")

celltype="Oligodendrocyte_precursor_cell"
adata=sc.read(f"{indir}/HCA_ON/data/5_refine_major/scvi/{celltype}/clean/{celltype}_subclass_seurat.h5ad")

outdir=f"{indir}/HCA_ON/data/5_refine_major/scvi/{celltype}/clean/"
sc.settings.figdir = f"{outdir}/figures"


mk={"OL Lingeage": ["SOX10","OLIG1","OLIG2","PLP1","MYRF","MOG","MAG"],"OPC":["PTPRZ1","PDGFRA","VCAN","BMP4"], "COP":["GPR17","BCAS1"], "NFOL":[ "RRAS2","PROM1"], "MFOL":["THBS3","SLC9A3R2"],"MOL1":["FOS"],"MOL2":["KLK6","HOPX"],"MOL5/6":["PTGDS","IL33"],"DA1":["SERPINA3","C4B"],"DA2":["CDKN1A","TNFRSF12A"],"IFN":["IRF9","IFIT1"]}

sc.pl.dotplot(adata,mk,groupby="subclass",save="oligo_opc_mm_mk.png")

mk1={"OPC":["OLIG1", "OLIG2", "SOX10", "MBP", "MOG", "PTPRZ1", "SOX6", "CA10", "PDGFRA"],"COP": ["CSPG4", "GPR17", "BMPER", "PRICKLE1"],"Oligo1":["SLC5A11", "RBFOX1","RASGRF1" , "RASGRF2", "ABCA6", "PLEKHG1", "ACTN2", "ACSBG1", "LGALS1", "ELOVL2"], "Oligo2":["PLXDC2","LINC01608", "QDPR", "DYSF", "OPALIN", "LAMA2"],"Oligo3":[ "SVEP1", "FRY", "SGCZ", "KCNK10",  "RRAS2"], "Oligo4":["HSPH1","ELOVL5", "HSPD1", "FOS", "DNAJA4"]}

sc.pl.dotplot(adata,mk1,groupby="subclass", save="oligo_opc_hg_ms_mk.png")

mk2={ "OPC": ["OLIG1", "OLIG2", "SOX10", "MBP", "MOG", "PTPRZ1", "CA10", "PDGFRA", "CSPG4"], "COP": ["SOX6", "GPR17", "BMPER", "PRICKLE1"], "Oligo1":[ "SLC5A11", "RASGRF2", "ACTN2", "RASGRF1", "ABCA6", "PLEKHG1", "RBFOX1", "ACSBG1", "LGALS1", "ELOVL2"],"Oligo2":["QDPR", "OPALIN", "PLXDC2", "LINC01608", "DYSF", "SVEP1", "LAMA2"]}

sc.pl.dotplot(adata,mk2, groupby="subclass", save="oligo_opc_hg_ad_mk.png")

mk={"other":["GAP43","LAMA2","EBF1","FOXP1"],"OPC":["PDGFRA","CSPG4","OLIG2"],"COP":["GPR17","VCAN","NEU4","BMP4","FYN","SOX6","PTPRZ1","PDCD4"],"NFOL":["TCF7L2","CEMIP2"],"MFOL":["MBP","OPALIN","CTPS1"],"MOL":["MBP","PLP1","MOG","GRM3","KLK6","PTGDS","ANXA5","HOPX"]}

sc.pl.dotplot(adata,mk, groupby="subclass", save="oligo_opc_andrew_mk.png")

from matplotlib import rcParams
rcParams["figure.figsize"] = (8,5)


#sc.pl.embedding(adata_query, basis="X_umap", color=["leiden"],legend_loc="on data", ncols=1,frameon=False,save=f'_{bname}_leiden_onData.png', palette="tab20")
bname="oligo_opc_hvg10000_epochnone_seurat_rs_0.4_clean"
adata=sc.read(f"{indir}/HCA_ON/data/5_refine_major/scvi/{celltype}/clean/{bname}_sb_none_scvi_trg.h5ad")

sc.pl.embedding(adata, basis="X_umap", color=["subclass"],legend_loc="on data", ncols=1,frameon=False,save=f'_{bname}_subclass_onData.png', palette="tab20")

rcParams["figure.figsize"] = (16,5)

adata=sc.read(f"{indir}/HCA_ON/data/5_refine_major/scvi/{celltype}/clean/{celltype}_subclass_seurat.h5ad")
#mk={"OPC Lineage": ["PTPRZ1","PDGFRA","VCAN","SOX6", "CA10","DCC","NRXN1","TNR","C1orf61"],"OPC_GLIS3+":["GLIS3","GAP43","UBASH3B"],"COP":["TCF7L2"], "COP_LAMA2":["UTRN","UACA","LAMA2","EBF1"],"NFOL":["BCAS1"], "COP_NFOL": ["GPR17","BMPER", "PRICKLE1"],"MFOL":["RNF220","ST18","PLP1","MBP"],"MOL":["MOG","GRM3","PTGDS"],"OLIGO_LRRC7+":["LRRC7"],"Oligo1":["SLC5A11","RASGRF2","ACSBG1","ACTN2","RASGRF1"], "OLIGO_SVEP1+":["SVEP1","PCDH9-AS2"],"OLICO_RBFOX1+":["AFF3","SNX24","RBFOX1"]}
mk={"OPC Lineage": ["PTPRZ1","PDGFRA","VCAN","SOX6", "CA10","DCC","NRXN1","TNR","C1orf61"],"OPC_GLIS3+":["GLIS3","GAP43","UBASH3B"],"COP":["TCF7L2"], "COP_LAMA2":["UTRN","UACA","LAMA2","EBF1"],"NFOL":["BCAS1"], "COP_NFOL": ["GPR17","BMPER", "PRICKLE1"],"MFOL":["RNF220","ST18","PLP1","MBP"],"MOL":["MOG","GRM3","PTGDS"],"OLIGO2_LRRC7+":["LRRC7"],"OLIGO1":["SLC5A11","RASGRF2","ACSBG1","ACTN2","RASGRF1"], "OLIGO1_SVEP1+":["SVEP1","PCDH9-AS2"],"OLIGO1_RBFOX1+":["AFF3","SNX24","RBFOX1"]}

sc.pl.dotplot(adata,mk, groupby="subclass", save="oligo_opc_own_mk_var.pdf",categories_order=["OPC","OPC_GLIS3+","COP_LAMA2+","COP_NFOL","MFOL","OLIGO2_LRRC7+","OLIGO1_SVEP1+","OLIGO1_RBFOX1+","OLIGO1"], standard_scale="var")

#mk={"OPC Lineage": ["PTPRZ1","PDGFRA","VCAN","SOX6", "CA10","DCC","NRXN1","TNR","C1orf61"],"OPC_GLIS3+":["GLIS3","GAP43","UBASH3B"], "OPC_Cycling":["CENPP", "PRIM2", "HELLS"], "COP_NFOL": ["GPR17","BMPER", "PRICKLE1", "TCF7L2"], "COP_LAMA2":["UTRN","UACA","LAMA2","EBF1"],  "MFOL":["RNF220","ST18","PLP1","MBP"],"MOL":["MOG","GRM3","PTGDS"],"OLIGO2_LRRC7+":["LRRC7"],"OLIGO1":["SLC5A11","RASGRF2","ACSBG1","ACTN2","RASGRF1"], "OLIGO1_SVEP1+":["SVEP1","PCDH9-AS2"],"OLIGO1_RBFOX1+":["AFF3","SNX24","RBFOX1"]}

mk={"OPC Lineage": ["PTPRZ1","PDGFRA","SOX6", "DCC"],"OPC_GLIS3+":["GLIS3","GAP43"], "OPC_Cycling":["CENPP", "PRIM2"], "COP_NFOL": ["GPR17","TCF7L2"], "NFOL_UTRN+":["UTRN","EBF1"],  "MFOL":["RNF220","MBP"],"MOL":["MOG","GRM3"],"OLIGO2_LRRC7+":["LRRC7"],"OLIGO1":["SLC5A11","ACSBG1"], "OLIGO1_SVEP1+":["SVEP1","PCDH9-AS2"],"OLIGO1_RBFOX1+":["SNX24","RBFOX1"]}
dt={"COP_LAMA2+":"NFOL_UTRN+"}
adata.obs["subclass"]=adata.obs["subclass"].replace(dt)
sc.pl.dotplot(adata,mk, groupby="subclass", save="oligo_opc_own_mk_var_cycling_rmRPE.pdf",categories_order=["OPC","OPC_GLIS3+","OPC_Cycling", "COP_NFOL","NFOL_UTRN+", "MFOL","OLIGO2_LRRC7+","OLIGO1_SVEP1+","OLIGO1_RBFOX1+","OLIGO1"], standard_scale="var")

mk={"Myelinating Oligos":["MOG","GRM3","PTGDS"], "OLIGO2_LRRC7+":["DCC", "LRRC7"],"OLIGO1":["SLC5A11","RASGRF2","ACSBG1","ACTN2","RASGRF1"], "OLIGO1_SVEP1+":["SVEP1","PCDH9-AS2"],"OLIGO1_RBFOX1+":["AFF3","SNX24","RBFOX1"]}

sc.pl.dotplot(adata,mk, groupby="subclass", save="oligo_only_own_mk_var.pdf",categories_order=["OLIGO1","OLIGO1_SVEP1+","OLIGO1_RBFOX1+", "OLIGO2_LRRC7+"]) #, standard_scale="var")

mk={"OPC Lineage": ["PTPRZ1","PDGFRA","VCAN","SOX6", "CA10","DCC","NRXN1","TNR","C1orf61"],"OPC_GLIS3+":["GLIS3","GAP43","UBASH3B"],"COP":["TCF7L2"], "COP_LAMA2":["UTRN","UACA","LAMA2","EBF1"], "COP_NFOL": ["GPR17","BMPER", "PRICKLE1"],"MFOL":["RNF220","ST18","PLP1","MBP"]}

sc.pl.dotplot(adata,mk, groupby="subclass", save="opc_only_own_mk_var.pdf",categories_order=["OPC","OPC_GLIS3+","COP_LAMA2+","COP_NFOL","MFOL"])
bname="Oligodendrocyte_precursor_cell_subclass_seurat"
sc.pl.violin(adata,keys=["pANN", "pct_counts_mt","n_genes_by_counts", "total_counts"], groupby="subclass", rotation=75, save=f"{bname}_QC.png", multi_panel=True)

#, standard_scale="var")

mk={"OPC Lineage": ["PTPRZ1","PDGFRA","VCAN","SOX6", "CA10","DCC","NRXN1","TNR","C1orf61"], "COP": ["GPR17","BMPER", "PRICKLE1"], "NFOL":["BCAS1","TCF7L2"],  "MFOL":["PLP1","MBP"],"MOL":["MOG","GRM3","PTGDS"],"hOligo1":["SLC5A11","RASGRF2","ACSBG1","ACTN2","RASGRF1"], "hOligo2":["PLXDC2"]}
sc.pl.dotplot(adata,mk, groupby="subclass", save="oligo_opc_own_mk.png",categories_order=["OPC","OPC_GLIS3+","COP_LAMA2+","COP_NFOL","MFOL","OLIGO2_LRRC7+","OLIGO1_SVEP1+","OLIGO1_RBFOX1+","OLIGO1"])

