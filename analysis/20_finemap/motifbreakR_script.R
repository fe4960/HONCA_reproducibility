library(motifbreakR)
library(MotifDb)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
#genomes
#search.genome = "BSgenome.Hsapiens.UCSC.hg19"

#snps <- snps.from.file(file = "sc_human_retina/data/motif/sc_eQTL_FDR10_20_50_caQTL_FDR05_10_20_ASCA_FDR10_20.bed",
#                                 search.genome = BSgenome.Hsapiens.UCSC.hg19,
#                                 format = "bed")

#snps.indels <- snps.from.file(file = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/13_GWAS/rs7636574.bed",
snps.indels <- snps.from.file(file = "/dfs3b/ruic20_lab/junw42/HCA_ON/data/13_GWAS/rs10811645.bed",

                                 search.genome = BSgenome.Hsapiens.UCSC.hg19,
                                 format = "bed", indels = FALSE)


#snps.mb.frombed[nchar(snps.mb.frombed$REF) > 1 | nchar(snps.mb.frombed$ALT) > 1]

#snps.indels <- snps.from.file(file = "sc_human_retina/data/motif/sc_eQTL_FDR10_20_50_caQTL_FDR05_10_20_ASCA_FDR10_20.vcf.gz",
#                                 search.genome = BSgenome.Hsapiens.UCSC.hg19,
#                                 format = "vcf")



#data(hocomoco)
data(motifbreakR_motif)
results <- motifbreakR(snpList = snps.indels, filterp = TRUE,
                       pwmList = motifbreakR_motif,
                       threshold = 1e-4,
                       method = "ic",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25))


#rs11024102 <- results[results$geneSymbol == "LHX2"]
#rs11024102 <- results[mcols(results)$geneSymbol == "LHX2"]
#rs11024102 <- results[results$geneSymbol == "LHX2", ]

#idx <- !is.na(mcols(results)$geneSymbol) & mcols(results)$geneSymbol == "LHX2"
#rs11024102 <- results[idx]


#pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/13_GWAS/rs7636574.motifbreakR.pdf", height=25)
#plotMB(results = results, rsid = "chr3:56871432:A:G", effect = "strong")
#dev.off()


#                       BPPARAM = BiocParallel::bpparam())
#df=data.frame(results)
df <- as.data.frame(results, row.names = NULL)
#df <- as.data.frame(results)
df[] <- lapply(df, as.character)

write.table(df,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/13_GWAS/rs10811645.motifbreakR",sep="\t",quote=F)
#write.table(df,file="/dfs3b/ruic20_lab/junw42/HCA_ON/data/13_GWAS/rs7636574.motifbreakR",sep="\t",quote=F)
