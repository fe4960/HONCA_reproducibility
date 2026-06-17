library(ggplot2)
library(dplyr)
library(tidyr)


mtcars=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor_Asian",header=T)

n_cyl <- mtcars %>%
  group_by(gender) %>%
  summarise(N = n()) %>%
  mutate(gender = as.factor(gender))


n_cyl2 <- n_cyl %>%
  arrange(desc(gender)) %>%
  mutate(pos = cumsum(N) - N/2,
         prop = paste(100*round(N/sum(N), 2), "%"))

pdf(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor_gender_color.pdf"))
ggplot(n_cyl2, aes(x = "x", y = N, fill = gender)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_text(aes(label = prop), position = position_stack(vjust = 0.5)) +  # Center-align the vertical positions of the labels on each stack
  coord_polar(theta = "y", direction = -1, clip = "off") +
  scale_y_continuous(breaks = n_cyl2$pos, labels = c("Male     ", "Female")) + 
#  scale_fill_brewer(name = NULL, guide = F, palette="Set3",direction=-1) +
  scale_fill_manual(name = NULL, guide = F, values=c("yellow","slateblue")) +  

  #  scale_fill_manual(name = NULL, guide = F, values=c("yellow","mediumpurple")) +  
  labs(title = "Number of Donor \n (n = 83)") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        axis.text.x = element_text(size = 25),
	text= element_text(size = 25))

dev.off()


n_cyl <- mtcars %>%
  group_by(race) %>%
  summarise(N = n()) %>%
  mutate(race = as.factor(race))


n_cyl2 <- n_cyl %>%
  arrange(desc(race)) %>%
  mutate(pos = cumsum(N) - N/2,
         prop = paste(100*round(N/sum(N), 2), "%"))

pdf(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor_race_color.pdf"))
ggplot(n_cyl2, aes(x = "x", y = N, fill = race)) +
  geom_bar(stat = "identity", position = "stack") + 
  geom_text(aes(label = prop), position = position_stack(vjust = 0.5)) +  # Center-align the vertical positions of the labels on each stack
  coord_polar(theta = "y", direction = -1, clip = "off") +
#  scale_y_continuous(breaks = n_cyl2$pos, labels = c("Male     ", "Female")) + 
#  scale_fill_brewer(name = NULL, guide = F, palette="Set3") +  
	  scale_fill_manual(name = NULL, guide = F, values=c("yellow","slateblue","salmon1","olivedrab3")) +

	  #	  scale_fill_manual(name = NULL, guide = F, values=c("yellow","mediumpurple","salmon","olivedrab3")) +
  labs(title = "Number of Donor \n (n = 83)") + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 25),
        axis.text.x = element_text(size = 25),
	text= element_text(size = 25))

dev.off()

pdf(paste0("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor_age_color.pdf"),width=8,height=4)
ggplot(mtcars, aes(x=age_year))+geom_histogram(fill="slateblue4", bins = 50)+theme(text = element_text(size=30),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"))
dev.off()


cellnum=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/major_ON_ONH_clean_new.csv",sep=",",header=T)

cellnum1=cellnum[!(cellnum$major.class%in%c("BC","Rod","Cone","HC","AC","RGC","Adipocyte","Pigmented_cell","RPE","Schwann_cell")),]

library(ggplot2)
library(ggbreak)

# Sample data
#df <- data.frame(
#  category = c("A", "B", "C", "D"),
#  value = c(10, 15, 100, 12)  # 'C' is an outlier
#)
cellnum1=cellnum1[order(cellnum1$cell.number,decreasing = TRUE),]
pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/major_ON_ONH_clean_new_cellnum_bar_color.pdf",width=9,height=8)

#pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/major_ON_ONH_clean_new_cellnum_bar.pdf",width=9,height=8)
# Standard bar plot with axis break
ggplot(cellnum1, aes(x = factor(major.class,levels=cellnum1$major.class), y = cell.number, fill = major.class)) +
  geom_bar(stat = "identity") +
#  scale_fill_manual(values = c("steelblue", "tomato", "gold", "forestgreen")) +
  scale_y_break(c(130000, 200000), scales = 0.5) + theme(text = element_text(size=30),legend.position="none", panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(), panel.border=element_blank(),axis.line=element_line(colour="black"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+scale_fill_brew(palette="paired")  # Truncates Y-axis between 20 and 90
#  theme_minimal() +
#  labs(title = "Truncated Bar Plot", y = "Value", x = "Category")
dev.off()

