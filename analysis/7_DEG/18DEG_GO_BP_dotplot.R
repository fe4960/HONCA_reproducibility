library(ggplot2)
library(stringr)
cell=c("all","Astrocyte","Endothelial_cell","Fibroblast","Oligodendrocyte")
label=c("up","down")
t="GO_Biological_Process_2023"
m="age"
dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass"
for(c in cell){
#	for(l in label){
		l="up"
		
		file=paste0(dir,"/",c,"_subclass_DEG/",c,"/",c,"_GO_top10_",t,"_",m,"_",l,"_column")
		data=read.table(file,header=F,sep="\t")
		data$sig="lightpink"
                data[data$V3>=-log(0.05,base=10),]$sig="red"
                p=ggplot(data=data,aes(x=V1,y=V2))+geom_point(aes(size=V3),color=data$sig)+theme_minimal()+theme(text=element_text(size=20),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+scale_size_continuous(name=expression(paste('-',log[10],'P')),range=c(2,15))+scale_y_discrete(labels = function(y) str_wrap(y, width = 40))
		if((c=="Fibroblast")&&(l=="down")){
		ht=height=length(data$V1)/9+5
		}else if(c == "all"){
 		ht=height=length(data$V1)/11+1
		}else{
		ht=height=length(data$V1)/9+4
		}
		pdf(paste0(dir,"/",c,"_subclass_DEG/",c,"/",c,"_GO_top10_",t,"_",m,"_",l,"_column.pdf"),width=12,height=ht)
		print(p)
		dev.off()

		l="down"
		file=paste0(dir,"/",c,"_subclass_DEG/",c,"/",c,"_GO_top10_",t,"_",m,"_",l,"_column")
		data=read.table(file,header=F,sep="\t")
		data$sig="lightblue"
                data[data$V3>=-log(0.05,base=10),]$sig="blue"
                p=ggplot(data=data,aes(x=V1,y=V2))+geom_point(aes(size=V3),color=data$sig)+theme_minimal()+theme(text=element_text(size=20),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+scale_size_continuous(name=expression(paste('-',log[10],'P')),range=c(2,15))+scale_y_discrete(labels = function(y) str_wrap(y, width = 40))
if((c=="Fibroblast")&&(l=="down")){
		ht=height=length(data$V1)/9+5
		}else if(c == "all"){
 		ht=height=length(data$V1)/11+1
		}else{
		ht=height=length(data$V1)/9+4
		}
		pdf(paste0(dir,"/",c,"_subclass_DEG/",c,"/",c,"_GO_top10_",t,"_",m,"_",l,"_column.pdf"),width=12,height=ht)
		print(p)
		dev.off()



#	}
}
