#import scvelo as scv
import sys
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import LinearSVC
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
import matplotlib.pyplot as plt
from sklearn.metrics import RocCurveDisplay, roc_curve, auc, roc_auc_score
from itertools import cycle
from sklearn.linear_model import LogisticRegression
from os.path import exists
import anndata as ad

from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)

#from matplotlib import scParams

#scParams("figuere.sizefig"=[5,5])

indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/sysvi/" #sys.argv[1]
h5ad=f"{indir}/sn_st_scvi_harmony.h5ad"
label="sn_st_scvi_harmony"
label1="majorclass"

#sch5ad=f"{indir}/mouse_CNS_immune.h5ad" #sys.argv[2]
#######snh5ad="/dfs3b/ruic20_lab/junw42/DRG/data/scvi/raw/major/major_DRG_3V4_6_5_25_exercise_ultima_merged_hvg10000_epochnone_seurat_rs1_sd7_sb_rm_clu3_kp_clu_11_scvi1_2_immune_scvi_trg.h5ad" #sys.argv[3]
#snh5ad="/dfs3b/ruic20_lab/junw42/DRG/data/scvi/raw/major/major_DRG_mouse_early_GARP_Lee_8632_flt_immune_clean_final_full.h5ad"
#label=sys.argv[4]
#label="mouse_immnue_sc_sn_8632"

####label="mouse_immnue_sc_sn"
#label1="labels2"
sc.settings.figdir=f"{indir}/figures/"
#label1="labels2"
#sc_adata=sc.read(f'{sch5ad}')
#sc_adata.var_names=sc_adata.var.gene_symbols.astype(str)
#sc_adata.var_names_make_unique()
#sc_adata.obs["tech"]="scRNA"
#sc_adata.obs[label1]=sc_adata.obs["author_cell_type"]


#sc_adata.obs["sampleid"]=sc_adata.obs["sample_id"]
#sc_adata.obs["donor"]=sc_adata.obs["donor"]
#sc_adata.obs[]
#sn_adata=sc.read(f'{snh5ad}')
#sn_adata.obs[label1]="unk"
#sn_adata.obs["tech"]="snRNA"
#gene=sc_adata.var.index.intersection(sn_adata.var.index)
#sc.pp.pca(adata, random_state=42)

#sc.pl.pca(adata, color="species", title="Species", save=f"_{label}_pca_spe.png", frameon=False)
#sc.pl.pca(adata, color="labels2", title="Cell Type", save=f"_{label}_pca_ct.png", frameon=False)
#adata=ad.concat([sn_adata[:, gene],sc_adata[:, gene]])
#sc.external.pp.harmony_integrate(adata, key='tech', basis="X_scVI", adjusted_basis='X_scVI_harmony')

#sc.pp.neighbors(adata, random_state=42, use_rep='X_pca_harmony')
#sc.pp.neighbors(adata, random_state=42, use_rep='X_scVI_harmony')

#sc.tl.umap(adata, random_state=42)

#sc.pl.umap(adata, color="tech", title="tech", save=f"_{label}_scVI_tech_harmony.png" , frameon=False )

#sc.pl.umap(adata, color=label1, title="Cell Type", save=f"_{label}_scVI_ct_harmony.png", frameon=False )
adata=sc.read(h5ad)
#adata.obs.loc[adata.obs["system"]=="xenium"].index=adata.obs.loc[adata.obs["system"]=="xenium","sampleid"].astype(str)+"_"+adata.obs.loc[adata.obs["system"]=="xenium"].index
#idx = adata.obs["system"] == "xenium"
#adata.obs.loc[idx].index = (
#    adata.obs.loc[idx, "sampleid"].astype(str) + "_" +
#    adata.obs.loc[idx].index
#)
rna=adata[adata.obs["system"]=="snRNA"]
rna.obs["barcode"]=rna.obs.index

rna.obs[label1]=rna.obs[label1].astype("category")
rna.obs[label1] = rna.obs[label1].cat.add_categories(["Immune_cell"])

rna.obs.loc[rna.obs[label1].isin(["T_cell","NK_cell", "Dendritic_cell", "B_cell", "Mast_cell"]), label1] = "Immune_cell"

rna.obs[label1] = rna.obs[label1].cat.remove_unused_categories()


atac=adata[adata.obs["system"]=="xenium"]
atac.obs["barcode"]=atac.obs.index
atac.obs.index=atac.obs["sampleid"].astype(str)+"_"+atac.obs.index
#rna0=sc.read(f'{indir}/{cell}_rna-emb.h5ad')
#atac0=sc.read(f'{indir}/{cell}_atac-emb.h5ad')

#dbt_table="/storage/chenlab/Users/junwang/human_meta/data/snATAC_clean_atac1perc_lr/emb_dbt.txt.gz"
#dbt_table=f'{indir}/{cell}_emb_dbt.txt.gz'

#rna=rna0
#atac=atac0
#rna_cell=rna0.obs.index
#atac_cell=atac0.obs.index

#rna0.obs.index=rna_cell
#atac0.obs.index=atac_cell

#if exists(dbt_table):
#	dbt=pd.read_csv(dbt_table,sep="\t",index_col=0)
#	rna_cell=rna0.obs.index.intersection(dbt.index)
#	atac_cell=atac0.obs.index.intersection(dbt.index)
#	rna=rna0[rna_cell,].copy()
#	atac=atac0[atac_cell,].copy()


random_state=42 #np.random.RandomState(1)

X=rna.obsm["X_scVI_harmony"]

#X=rna.obsm["X_glue"]
y=rna.obs[label1]
#X_atac=atac.obsm["X_glue"]
X_atac=atac.obsm["X_scVI_harmony"]

n_samples, n_features=X.shape
n_classes=len(np.unique(y))
#if n_classes==2 :
#	n_classes=1

target_names=np.unique(y)
##split the X_train and X_test into 50:50
(X_train, X_test, y_train, y_test) = train_test_split(X,y, test_size=0.5, stratify=y, random_state=42)

#######use SVM oneVsRestClassifier
#y_score_svm1=OneVsRestClassifier(LinearSVC(random_state=0)).fit(X_train, y_train).predict(X_test)
#y_score_svm=LabelBinarizer().fit_transform(y_score_svm1)

#y_score_svm_pred1=OneVsRestClassifier(LinearSVC(random_state=0)).fit(X, y).predict(X_atac)
#y_score_svm_pred=LabelBinarizer().fit_transform(y_score_svm_pred1)

#######use logistical regression
classifier= LogisticRegression(random_state=42, max_iter=1000)

#if n_classes==1:



y_score_lr=classifier.fit(X_train, y_train).predict_proba(X_test)
y_score_lr_pred=classifier.fit(X, y).predict_proba(X_atac)

#if n_classes==1 :
#	y_score_lr=y_score_lr[:,1]
#	y_score_lr_pred=y_score_lr_pred[:,1]
######binarize y and transform into one_hot_test
label_binarizer = LabelBinarizer().fit(y_train)
y_onehot_test = label_binarizer.transform(y_test)
y_onehot_test.shape 

####ROC curve using the OvR macro-average
# store the fpr, tpr, and roc_auc for all averaging strategies
fpr_lr, tpr_lr, roc_auc_lr = dict(), dict(), dict()
#fpr_svm, tpr_svm, roc_auc_svm = dict(), dict(), dict()
n=y_onehot_test.shape[1]
for i in range(n):
#	fpr_svm[i],tpr_svm[i], _ = roc_curve(y_onehot_test[:,i], y_score_svm[:,i])
#	roc_auc_svm[i] = auc(fpr_svm[i], tpr_svm[i])
	fpr_lr[i],tpr_lr[i], _ = roc_curve(y_onehot_test[:,i], y_score_lr[:,i])
	roc_auc_lr[i] = auc(fpr_lr[i], tpr_lr[i])



fpr_grid=np.linspace(0.0, 1.0, 1000)

#interpolate all ROC curves at these points
#mean_tpr_svm=np.zeros_like(fpr_grid)
mean_tpr_lr=np.zeros_like(fpr_grid)

for i in range(n_classes):
#	mean_tpr_svm += np.interp(fpr_grid, fpr_svm[i], tpr_svm[i])
	mean_tpr_lr += np.interp(fpr_grid, fpr_lr[i], tpr_lr[i])

#	mean_tpr += np.interp(fpr_grid, fpr[i], tpr[i])

# Average it an compute AUC 

#mean_tpr_svm /=n_classes

mean_tpr_lr /=n_classes


#fpr_svm["macro"] = fpr_grid
fpr_lr["macro"] = fpr_grid

#tpr_svm["macro"] = mean_tpr_svm
tpr_lr["macro"] = mean_tpr_lr

#roc_auc_svm["macro"] = auc(fpr_svm["macro"], tpr_svm["macro"])
roc_auc_lr["macro"] = auc(fpr_lr["macro"], tpr_lr["macro"])



fig, ax = plt.subplots(figsize=(6,6))

#plt.plot(
#	fpr_svm["macro"],
#	tpr_svm["macro"],
#	label=f'SVM macro-ave ROC curve (AUC = {roc_auc_svm["macro"]:.2f})',
#	color="navy",
#	linestyle=":",
#	linewidth=4,
#)


plt.plot(
	fpr_lr["macro"],
	tpr_lr["macro"],
	label=f'LR macro-ave ROC curve (AUC = {roc_auc_lr["macro"]:.2f})',
	color="red",
	linestyle=":",
	linewidth=4,
)


#colors =cycle(["aqua", "darkorange", "cornflowerblue"])
#for class_id, color in zip(range(n_classes), colors):
#	RocCurveDisplay.from_predictions(
#		y_onehot_test[:,class_id],
#		y_score_lr[:,class_id],
#		name=f'LR ROC curve for {target_names[class_id]}',
#		color=color,
#		ax=ax,
#	)

#colors=cycle(["red","blue","green"])
#colors=range(3)

#for class_id, color in zip(range(n_classes), colors):
#	RocCurveDisplay.from_predictions(
#		y_onehot_test[:,class_id],
#		y_score_svm[:,class_id],
#		name=f'SVM ROC curve for {target_names[class_id]}',
#		color=color,
#		ax=ax,
#	)

plt.plot([0,1],[0,1],"k--", label="ROC curve for chance level (AUC=0.5)")
plt.axis("square")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Extension of Receiver Operating Characteristic\nto One-vs-Rest multiclass")
plt.legend()
plt.savefig(f'{indir}/figures/{label}_auc.png')



#atac.obsm["svm_class"] =y_score_svm_pred
#df=pd.DataFrame(y_score_svm_pred, columns=range(y_score_svm_pred.shape[1]))
#atac.obs["svm_celltype"]=target_names[df.idxmax(axis=1)]

atac.obsm["lr_class"] =y_score_lr_pred
df=pd.DataFrame(y_score_lr_pred, columns=range(y_score_lr_pred.shape[1]))
atac.obs["lr_celltype"]=target_names[df.idxmax(axis=1)]

df1=df.max(axis=1)

df.to_csv(f'{indir}/{label}_score.txt.gz',sep="\t")

atac.obs.to_csv(f'{indir}/{label}_obs.txt.gz',sep="\t")
atac.obs.iloc[df1[df1<0.9].index].to_csv(f'{indir}/{label}_dbt.txt.gz',sep="\t")
#cm=atac.obs.groupby(["lr_celltype","svm_celltype"]).size().unstack(fill_value=0)
#cm.to_csv(f'{indir}/{cell}_lr_svm.txt.gz',sep="\t")

#combined0=sc.read(f'{indir}/{label}_combined-emb.h5ad')
#combined.obs["svm_celltype"]=[]
#combined.obs["lr_celltype"]=[]

#new_cell=rna_cell.append(atac_cell)

combined=ad.concat([atac,rna]) #combined0[new_cell].copy()

#combined.obs.loc[atac.obs.index,"svm_celltype"]=atac.obs["svm_celltype"]
combined.obs.loc[atac.obs.index,"lr_celltype"]=atac.obs["lr_celltype"]

#combined.obs.loc[rna.obs.index,"svm_celltype"]=rna.obs[label]
combined.obs.loc[rna.obs.index,"lr_celltype"]=rna.obs[label1]

combined=combined[combined.obs["lr_celltype"]!="NA"]

sc.pl.embedding(combined, basis="X_umap", color=["lr_celltype"],ncols=2,frameon=False,save=f'_{label}_combined_pred.png') #, palette="tab20")

combined.write(f"{indir}/{label}_coembed.h5ad")

#sc.pl.embedding(combined, basis="X_umap", color=["svm_celltype","lr_celltype"],ncols=2,frameon=False,save=f'_{cell}_combined_pred.png') #, palette="tab20")
