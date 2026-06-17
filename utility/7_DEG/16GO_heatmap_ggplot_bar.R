library(tidyr)
library(dplyr)

shorten_names <- function(names, width = 15) {
  sapply(names, function(x) {
    paste(strwrap(x, width = width), collapse = "\n")
  })
}



up=read.table("HCA_ON/data/7_DEG/genexp_donor_subclass/all_subclass_DEG/all/Oligodendrocyte_precursor_cell_DEG_GO_age_up/enrichr_GO_Biological_Process_2023.txt",header=T,sep="\t")
up=up %>% separate(Term,c("GO"),sep="\\s+\\(")
up=up[up$Adjusted.P.value<0.1,]
# Compute -log10(p-value) for bar length
up$logp <- -log10(up$P.value)
up$GO <- shorten_names(up$GO, width = 35)
up$p=-log10(up$P.value)
# Reorder terms for plotting
up$GO <- factor(up$GO, levels = up$GO[order(up$logp)])
up$type="Up"


down=read.table("HCA_ON/data/7_DEG/genexp_donor_subclass/all_subclass_DEG/all/Oligodendrocyte_precursor_cell_DEG_GO_age_down/enrichr_GO_Biological_Process_2023.txt",header=T,sep="\t")
down=down %>% separate(Term,c("GO"),sep="\\s+\\(")
down=down[down$Adjusted.P.value<0.1,]
# Compute -log10(p-value) for bar length
down$logp <- log10(down$P.value)

# Reorder terms for plotting
down$GO <- factor(down$GO, levels = down$GO[order(down$logp)])
down$type="Down"
write.table(down,"HCA_ON/data/7_DEG/genexp_donor_subclass/all_subclass_DEG/all/Oligodendrocyte_precursor_cell_DEG_GO_age_down/enrichr_GO_Biological_Process_2023.txt_simple_padj1",sep="\t",quote=F,row.names=F)
down1=read.table("HCA_ON/data/7_DEG/genexp_donor_subclass/all_subclass_DEG/all/Oligodendrocyte_precursor_cell_DEG_GO_age_down/enrichr_GO_Biological_Process_2023.txt_simple_padj1_short.txt",sep="\t", header=T)
down1$GO <- shorten_names(down1$GO, width = 35)
#down1$GO <- factor(down1$GO, levels = down1$GO[order(down1$logp, decreasing = T)])
down1=down1[order(down1$logp,decreasing = T),]
down1$p=log10(down1$P.value)
go_results=rbind(up,down1)

go_results$GO <- factor(go_results$GO, levels = go_results$GO)
# Apply the function
#go_results$GO <- shorten_names(go_results$GO, width = 30)

#go_results=go_results
# Plot using ggplot2
library(ggplot2)
pdf("HCA_ON/data/7_DEG/genexp_donor_subclass/all_subclass_DEG/all/Oligodendrocyte_precursor_cell_DEG_GO_age_up_down_ggplot_bar.pdf",width=11,height=10)
ggplot(go_results, aes(x = GO, y = logp,fill=factor(type))) +
  geom_bar(stat = "identity", aes(fill = factor(go_results$type))) +
  coord_flip() +
  theme_minimal(base_size = 20) +
  labs(
    x = "GO Term",
    y = expression(-log[10](p-value)),
#    title = "GO Enrichment Analysis"
  )+scale_fill_brewer(palette="Set2")
dev.off()
