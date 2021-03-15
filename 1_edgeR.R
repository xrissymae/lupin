###----------------------###

#title: "1. Transform Counts"

###----------------------###

library(edgeR)
library(ggplot2)
library(dplyr)
library(ggrepel)


# 1_P Transforming Count Files to logCPM/logFC
setwd("~/Documents/Research/Bioinformatics/Illumina/")
files <- c("counts/B_PPLUS.fq.gz.quant.counts", 
           "counts/C_PMIN.fq.gz.quant.counts",
           "counts/D_PMIN.fq.gz.quant.counts")

cpm_var <- as.numeric(".5")
rep_total <- as.numeric("1")

bcv = 0.1
data <- readDGE(files)
group <- c(rep("+P", 1), rep("-P", 2))
dge <- DGEList(counts=data, group=group)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
keep <- rowSums(cpm(dge)>cpm_var) >= rep_total
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")

et <- exactTest(dge, pair=c("+P", "-P"), dispersion = bcv^2)
etp <- topTags(et, n=100000)

write.csv(etp$table, "20210314_-P.csv")

# 1_P MA plot Parameters logFC >1, PVal < .05, FDR <.05

options(stringsAsFactors = FALSE)
results1 = read.csv("1_-P(test).csv")

results1 = mutate(results1, 
                  Genes=ifelse(results1$logFC>1
                               & results1$PValue<.05 
                               & results1$FDR<5e-30, 
                               "Potential Genes of Interest", 
                               ifelse(results1$FDR<0.05, 
                                      "Differentially Expressed", 
                                      "Non-Differentially Expressed")))

MA1 = ggplot(results1, aes(logCPM, logFC)) + geom_point(aes(col=Genes)) + scale_color_manual(values=c("#51A9A2", "#22475D", "#FB8610"))

#volcano plot
volcano1 <- ggplot(results1, aes(logFC, -log2(PValue))) + geom_point(aes(col=Genes))+ scale_color_manual(values=c("#51A9A2", "#22475D", "#FB8610"))

#volcano + top20 n123
target <- read.csv("top_p_n123.csv")
png(paste0("1_Volcano_-P", format(Sys.time(), "%y%m%d.%H%M.png")))
volcano1 + geom_point(data=filter(results1, 
                                  X %in% target$gene ), color="black") + ggtitle("-P vs. +P: Volcano Plot")
dev.off()
# Put points on that graph indicating that stuff from the venn diagram 

png(paste0("1_MA_-P_", format(Sys.time(), "%y%m%d.%H%M.png")))
MA1 + geom_point(data=filter(results1, 
                             X %in% target$gene ), color="black") + ggtitle("+P vs. -P: Differentially Expressed Genes")
dev.off()

# 2_-Fe

setwd("~/Desktop/Bioinformatics/Illumina/Trimmed/")
files2 <- c("counts/B_PPLUS.fq.gz.quant.counts", 
            "counts/E_FEMIN.fq.gz.quant.counts",
            "counts/F_FEMIN.fq.gz.quant.counts")

bcv = 0.1
data <- readDGE(files2)
group <- c(rep("+P", 1), rep("-Fe", 2))
dge <- DGEList(counts=data, group=group)
dge <- estimateCommonDisp(dge, dispersion=bcv^2)
dge <- estimateTagwiseDisp(dge)

et <- exactTest(dge, pair=c("+P", "-Fe"), dispersion = bcv^2)
etp <- topTags(et, n=100000)
write.csv(etp$table, "2_-Fe.csv")

# -Fe MA Plot

results2 = read.csv("2_-Fe.csv")

results2 = mutate(results2, 
                  Genes=ifelse(results2$logFC>1
                               & results2$PValue<0.05 
                               & results2$FDR<5e-30, 
                               "Potential Genes of Interest", 
                               ifelse(results2$FDR<0.05, 
                                      "Differentially Expressed", 
                                      "Non-Differentially Expressed")))

MA2 = ggplot(results2, aes(logCPM, logFC)) + geom_point(aes(col=Genes)) + scale_color_manual(values=c("#51A9A2", "#22475D", "#FB8610"))

# -Fe Volcano Plot

volcano2 <- ggplot(results2, aes(logFC, -log2(PValue))) + geom_point(aes(col=Genes))+ scale_color_manual(values=c("#51A9A2", "#22475D", "#FB8610"))

target <- read.csv("top_fe_n123.csv")

png(paste0("2_Volcano_-Fe_", format(Sys.time(), "%y%m%d.%H%M.png")))
volcano2 + geom_point(data=filter(results2, 
                                  X %in% target$gene ), color="black") + ggtitle("-Fe vs. +P: Volcano Plot")
dev.off()

png(paste0("2_MA_-Fe_", format(Sys.time(), "%y%m%d.%H%M.png")))
MA2 + geom_point(data=filter(results2, 
                             X %in% target$gene ), color="black") + ggtitle("-Fe vs. +P : Differentially Expressed Genes")
dev.off()

### + Suc Files
setwd("~/Desktop/Bioinformatics/Illumina/Trimmed/")
files <- c("counts/B_PPLUS.fq.gz.quant.counts", 
           "counts/G_SPLUS.fq.gz.quant.counts",
           "counts/H_SPLUS.fq.gz.quant.counts")

bcv = 0.1
data <- readDGE(files)
group <- c(rep("+P", 1), rep("+Suc", 2))
dge <- DGEList(counts=data, group=group)
dge <- estimateCommonDisp(dge, dispersion=bcv^2)
dge <- estimateTagwiseDisp(dge)

et <- exactTest(dge, pair=c("+P", "+Suc"), dispersion = bcv^2)
etp <- topTags(et, n=100000)
write.csv(etp$table, "3_+Suc.csv")

# +Suc MA Plot

results3 = read.csv("3_+Suc.csv")

results3 = mutate(results3, Genes=ifelse(results3$logFC>1
                                         & results3$PValue<0.05 
                                         & results3$FDR<5e-30, 
                                         "Potential Genes of Interest", 
                                         ifelse(results3$FDR<0.05, 
                                                "Differentially Expressed", 
                                                "Non-Differentially Expressed")))

MA3 = ggplot(results3, aes(logCPM, logFC)) + geom_point(aes(col=Genes)) + scale_color_manual(values=c("#51A9A2", "#22475D", "#FB8610"))

# +Suc Volcano Plot

volcano3 <- ggplot(results3, aes(logFC, -log2(PValue))) + geom_point(aes(col=Genes))+ scale_color_manual(values=c("#51A9A2", "#22475D", "#FB8610"))

target <- read.csv("top_s_n123.csv")

png(paste0("3_Volcano_+Suc_", format(Sys.time(), "%y%m%d.%H%M.png")))
volcano3 + geom_point(data=filter(results3, 
                                  X %in% target$gene ), color="black") + ggtitle("+Suc vs. +P: Volcano Plot")
dev.off()

png(paste0("3_MA_+Suc_", format(Sys.time(), "%y%m%d.%H%M.png")))
MA3 + geom_point(data=filter(results3, X %in% target$gene ), color="black") + ggtitle("+Suc vs. +P : Differentially Expressed Genes")
dev.off()


###--- Heatmaps ---###

##P
genetable <- read.csv("~/Desktop/Bioinformatics/Illumina/top_p_n123.csv")
genetable <- genetable[, c(1,5:8)]
genetable$Gene <- paste(genetable[,1], genetable[,2], sep=": ")
genetable <- genetable[,2:5]
genematrix <- as.matrix(genetable[2:4])
rownames(genematrix) <- genetable[,1]

png("1_top20_p_heatmap.png")
pheatmap(genematrix, main = "Heatmap: Top 20 -P DEGs")
dev.off()

##Fe
genetable2 <- read.csv("~/Desktop/Bioinformatics/Illumina/Trimmed/top_fe_n123.csv")
genetable2 <- genetable2[, c(1,5:8)]
genetable2$Gene <- paste(genetable2[,1], genetable2[,2], sep=": ")
genetable2 <- genetable2[,2:5]

head(genetable2)

genematrix2 <- as.matrix(genetable2[2:4])
rownames(genematrix2) <- genetable2[,1]

png("2_top20_fe_heatmap.png")
pheatmap(genematrix2, main = "Heatmap: Top 20 -Fe DEGs")
dev.off()

##Suc
genetable <- read.csv("~/Desktop/Bioinformatics/Illumina/Trimmed/top_s_n123.csv")
genetable <- genetable[, c(1,5:8)]
genetable$Gene <- paste(genetable[,1], genetable[,2], sep=": ")
genetable <- genetable[,2:5]

head(genetable)

genematrix <- as.matrix(genetable[2:4])
rownames(genematrix) <- genetable[,1]

png("3_top20_suc_heatmap.png")
pheatmap(genematrix, main = "Heatmap: Top 20 +Suc DEGs")
dev.off()

