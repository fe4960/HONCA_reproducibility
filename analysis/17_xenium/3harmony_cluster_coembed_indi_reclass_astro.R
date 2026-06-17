#meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean_simple.obs.gz"
meta="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple.obs.gz"
data=read.table(meta,header=T,sep=",")
data$subclass1=NA
data$subclass1[grepl("retina",data$subclass)]="astro_retina"
data$subclass1[grepl("ONHON",data$subclass)]="astro_ONHON"
data$subclass1[(!grepl("ONHON",data$subclass))& grepl("ONH1",data$subclass)]="astro_ONH1"
data$subclass1[(!grepl("ONHON",data$subclass))& grepl("ONH2",data$subclass)]="astro_ONH2"
data$subclass1[(!grepl("ONHON",data$subclass))& grepl("ON",data$subclass)& ( !grepl("ONH",data$subclass))]="astro_ON"


data$subclass2=NA
data$subclass2[grepl("retina",data$subclass)]="astro_retina"
data$subclass2[grepl("ONHON",data$subclass)]="astro_ON"
data$subclass2[(!grepl("ONHON",data$subclass))& grepl("ONH",data$subclass)]="astro_ONH"
data$subclass2[(!grepl("ONHON",data$subclass))& grepl("ON",data$subclass)& ( !grepl("ONH",data$subclass))]="astro_ON"

write.table(data,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple_sort.obs",quote=F,sep=",",row.names=F)
