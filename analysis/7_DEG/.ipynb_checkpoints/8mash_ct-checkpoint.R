library(mashr)
library(ashr)
set.seed(1)
#trait = c("raceAfrican", "raceAsian", "raceHispanic", "age", "gender", "ageGender")
trait=c("age")
cell=c("all")
#cell=c("Oligodendrocyte")
#cell=c("Endothelial_cell","Astrocyte","Fibroblast")
for(ct in cell){
#genelist=read.table("/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/genexp_sample_raw_celltype_atlas_sample/cpm_1_gene_ct.txt",header=T)
dir=paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/",ct,"_subclass_DEG/")
for(tr in trait){
Bhat=read.table(paste0(dir,"/",tr,"_DEG_res_cpm_",tr,"_beta_ancFix"),header=T)
Shat=read.table(paste0(dir,"/",tr,"_DEG_res_cpm_",tr,"_ste_ancFix"),header=T)

Bhat=data.matrix(Bhat)
Shat=data.matrix(Shat)
#Z=data.matrix(Z)

#Bhat=Bhat[genelist$x,]
#Shat=Shat[genelist$x,]
#Z=Z[genelist$x,]
##########
data   = mash_set_data(Bhat,Shat, zero_Bhat_Shat_reset=0.000001)

#data   = mash_set_data(Z) #Bhat,Shat)
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)
if(ct == "Astrocyte"){
U.pca = cov_pca(data,3,subset=strong)
}else{
U.pca = cov_pca(data,5,subset=strong)
}
print(names(U.pca))

U.ed = cov_ed(data, U.pca, subset=strong)
m.ed = mash(data, U.ed)
print(get_loglik(m.ed),digits = 10)

U.c=cov_canonical(data)  
m=mash(data, c(U.c,U.ed))
print(get_loglik(m),digits = 10)
saveRDS(m,file=paste0(dir,"/",tr,"_DEG_res_cpm_",tr,"_m_ct_ancFix.rds"))
#saveRDS(m,file=paste0("/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/",tr,"_DEG_res_cpm_",tr,"_m_float.rds"))

tbl=get_lfsr(m)
write.table(tbl,file=paste0(dir,"/",tr,"_DEG_res_cpm_",tr,"_lfsr_ct_ancFix"),quote=F,sep="\t",row.names=F)

tbl=get_pm(m)
write.table(tbl,file=paste0(dir,"/",tr,"_DEG_res_cpm_",tr,"_pm_ct_ancFix"),quote=F,sep="\t",row.names=F)

tbl=get_psd(m)
write.table(tbl,file=paste0(dir,"/",tr,"_DEG_res_cpm_",tr,"_psd_ct_ancFix"),quote=F,sep="\t",row.names=F)

head(get_significant_results(m))

print(length(get_significant_results(m)))

tbl=get_pairwise_sharing(m)
write.table(tbl,file=paste0(dir,"/",tr,"_DEG_res_cpm_",tr,"_cell_eff_ct_ancFix"),quote=F,sep="\t",row.names=F)


tbl=get_pairwise_sharing(m, factor=0)
write.table(tbl,file=paste0(dir,"/",tr,"_DEG_res_cpm_",tr,"_cell_dir_ct_ancFix"),quote=F,sep="\t",row.names=F)
}
}
