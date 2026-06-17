import pandas as pd

data_full=pd.DataFrame()

for i in range(0,9):
    data=pd.read_csv(f"software/deepsas/outputs/Data_{i}_/Senescent_Tables/Cell_Table1_SnC_scores.csv", header=0, index_col=1)
    data_full=pd.concat([data, data_full])


obs=pd.read_csv("HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_simple1.obs.gz", header=0, index_col=0)

obs=obs.loc[data_full.index]
obs["ifsnc"]=data_full["ifSnCs"]
obs["celltype"]=data_full["cell_type"]
#prop=obs.groupby(["celltype", "donor", "age_year"])["ifsnc"].mean()
prop=obs.groupby(["celltype", "sampleid", "age_year"])["ifsnc"].mean()

#prop.to_csv("software/deepsas/outputs/snc_age_ct_prop")
prop.to_csv("software/deepsas/outputs/snc_age_ct_prop_sampleid")

prop=obs.groupby(["sampleid", "age_year"])["ifsnc"].mean()

prop.to_csv("software/deepsas/outputs/snc_age_class_prop_sampleid")

