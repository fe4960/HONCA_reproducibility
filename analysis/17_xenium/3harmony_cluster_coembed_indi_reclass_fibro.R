#meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_simple.obs.gz"
#meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple.obs.gz"
meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_simple.obs.gz"


data=read.table(meta,header=T,sep=",")
data$subclass1=NA
data$subclass1[grepl("arachnoid",data$subclass)]="Fibro_arachnoid"
data$subclass1[grepl("dura",data$subclass)]="Fibro_dura"
data$subclass1[grepl("perivascular",data$subclass)]="Fibro_peri"
data$subclass1[grepl("pia",data$subclass)]="Fibro_pia"
data$subclass1[grepl("RPEchoroid",data$subclass)]="Fibro_RPEchoroid"
data$subclass1[grepl("sclera",data$subclass)]="Fibro_sclera"


write.table(data,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_simple_sort.obs",quote=F,sep=",",row.names=F)
