library(ggplot2)
data=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta",sep=",",header=T)
dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/meta_sum/figures/"

df=data %>% group_by(c("majorclass1","sampleid","tissue","source"))

#df <- df %>% 
#  mutate(prop = v / sum(df$v), 
#         label= scales::percent(prop, accuracy = 0.1)) %>% 
#  arrange(prop) %>% 
#  mutate(ypos = cumsum(prop)- 0.5*prop, 
#         SubSegment = factor(SubSegment, levels=SubSegment[order(-(prop))], ordered=TRUE))



#ggplot(df, aes(x = "", y = prop, fill = SubSegment)) +
#  geom_bar(width = 1, stat = "identity", color="white", alpha=0.8) +
#  coord_polar("y", start = 0) +
#  theme_void() + 
#  geom_text(aes(y = ypos, label = label), size=3, color = "white") +
#  scale_fill_brewer(palette="Set1")

#p=ggplot(data=data,aes(x="",y=gender.Freq,fill=gender))+geom_bar(width=1, stat="identity", color="white")+coord_polar("y",start=0)+theme(text = element_text(size=15),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_blank())+guides(fill=guide_legend(title="Gender"))+scale_fill_brewer(direction=-1,palette="Set3")+xlab("")+ylab("")
#pdf(paste0(dir,"ononh_final_new_gender.pdf"),width=15)
#print(p)
#dev.off()

data1=unique(data[,c("donor","age_year")])
p=ggplot(data=data1,aes(x=age_year))+geom_histogram(binwidth=1)+theme(text = element_text(size=25),panel.background=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_blank())+guides(fill=guide_legend(title="Age distribution"))+xlab("Age")+ylab("Number of donors")
pdf(paste0(dir,"ononh_final_new_age.pdf"),width=15)
print(p)
dev.off()


data_c_on=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass_proportion_sample_ON_Chen.csv",sep=",",header=T)
data_s_on=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass_proportion_sample_ON_Sanes.csv",sep=",",header=T)
data=rbind(data_c_on,data_s_on)
ct=c("T_cell", "NK_cell", "B_cell", "Mast_cell","Dendritic_cell","Immune_neuron_unk")
my_palette <- c("purple", "orange", "darkseagreen1","grey","blue", "thistle","lightblue1", "yellow","darkgreen","red","darkolivegreen2","deepskyblue","lightyellow","pink","chartreuse","lavenderblush","aquamarine","gold","coral1","beige","darkseagreen3")
for(c in ct){
data$majorclass1=gsub(c,"Immu",data$majorclass1)
}
#data=data[order(data$sampleid),]
pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass_proportion_sample_ON.pdf",width=15,height=10)
p=ggplot(data=data,aes(x=sampleid,y=nCount_RNA,fill=majorclass1))+geom_bar(position="fill", stat="identity")+theme(text = element_text(size=25),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_blank(), axis.text.x =element_blank())+scale_fill_manual(values=my_palette)+ylab("Cell class proportion")
print(p)
dev.off()

data_c_on=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass_proportion_sample_ONH_Chen.csv",sep=",",header=T)
data_s_on=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass_proportion_sample_ONH_Sanes.csv",sep=",",header=T)
data=rbind(data_c_on,data_s_on)
ct=c("T_cell", "NK_cell", "B_cell", "Mast_cell","Dendritic_cell","Immune_neuron_unk")
my_palette <- c("purple", "orange", "darkseagreen1","grey","blue", "thistle","lightblue1", "yellow","darkgreen","red","darkolivegreen2","deepskyblue","lightyellow","pink","chartreuse","lavenderblush","aquamarine","gold","coral1","beige","darkseagreen3")
for(c in ct){
data$majorclass1=gsub(c,"Immu",data$majorclass1)
}
#data=data[order(data$sampleid),]
pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/ONONH_new_final_majorclass_proportion_sample_ONH.pdf",width=15,height=10)
p=ggplot(data=data,aes(x=sampleid,y=nCount_RNA,fill=majorclass1))+geom_bar(position="fill", stat="identity")+theme(text = element_text(size=25),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_blank(), axis.text.x =element_blank())+scale_fill_manual(values=my_palette)+ylab("Cell class proportion")
print(p)
dev.off()
