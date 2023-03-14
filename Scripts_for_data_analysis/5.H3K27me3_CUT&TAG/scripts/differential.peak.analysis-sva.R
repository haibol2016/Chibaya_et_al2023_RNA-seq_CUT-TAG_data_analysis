# BiocManager::install(c("DESeq2", "edgeR", "sva", "limma", "genefilter"))
BiocManager::install('ashr')

library("sva")
library("ashr")
library("EnhancedVolcano")
library("DESeq2")
library("edgeR")
library("qsmooth")
library("quantro")
library("limma")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library("ggfortify")
library("genefilter")
library("AnnotationDbi")
library("WriteXLS")
library("gplots")
library("doParallel")

#setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\UMass_work\Marcus\cut_N_tag-03112022)")

## DESeq2 analysis
count <- read.delim("10272022.CutTag.count.table.txt", 
                    header = TRUE, as.is =TRUE, check.names = FALSE)

colnames(count) <- gsub("-", "_", colnames(count), perl = TRUE)

meta_data <- data.frame(SampleID = colnames(count),
                       Group  = factor(gsub("_?\\d+$", "", colnames(count))))

rownames(meta_data) <- meta_data$SampleID
all(colnames(count) == rownames(meta_data))

## plot raw count per sample
pdf("Fig 1. Number of reads assigened to consensus peaks.pdf", 
    width = 6, height = 4)
barplot(colSums(count)/1000000, width = 0.5, 
        ylim=c(0, 200),
        ylab = "Reads (Million)", 
        col = rep(rainbow(4), times = c(4,3,4,2)),
        cex.axis = 0.5,
        cex.names = 0.5,
        names.arg = colnames(count), las=2)
dev.off()

## order the metadata according to colnames of count table

all(colnames(count) == rownames(meta_data))

## filter count table to remove the lowly expressed genes
cpms = cpm(count)
keep = rowSums(cpms > 3) >= 2   #141,194 peaks
count <- count[keep, ]



#### sample similarity
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count),
                              colData = meta_data,
                              design = ~0 + Group)

dds <- estimateSizeFactors(dds)

dds <- estimateDispersions(dds, fitType="parametric", 
                           maxit=1000)

### exploratory analysis
vsd <- vst(dds, blind = TRUE)
sampleDists <- dist(t(assay(vsd)))

#### Heatmap showing sample distances
#### samples SP4 and SP8 are very different from SP5-7 of the same treatment
distancePlot <- function(sampleDists, sampleNames)
{
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- sampleNames
  colnames(sampleDistMatrix) <- sampleNames
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  sampleDistMatrix 
}

## one sample PIL_C5 is an outlier, so it should be removed and reanalyze it.
pdf("Fig 4.Heatmap showing original sample distances.pdf", width = 6, height = 5)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = vsd$SampleID)
dev.off()


## PCA using adjusted value
pc <- prcomp(t(assay(vsd)), scale = T, center = T, retx = TRUE)
pdf(file = "Fig 5. PCA.plot.pdf", width = 4.0, height = 2.5)
autoplot(pc, data = meta_data, colour = 'Group', size = 2)
dev.off()

normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
write.table(normCounts, "Table 1. Normalized.count.table.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

# using sva to reduce hidden variations
get.sva <- function (expr.data=NULL, meta.data=NULL)
{
  mod <- model.matrix(~ Group, data=meta.data)
  mod0 <- model.matrix(~ 1, data=meta.data)
  
  num.sva <-svaseq(as.matrix(expr.data), mod, mod0)$n.sv
  sv <- svaseq(as.matrix(expr.data), mod, mod0,n.sv=num.sva)$sv
  
  colnames(sv)<- paste0("sv",1:num.sva)
  
  meta.data.sva <-cbind(meta.data,sv )
  
  meta.data.sva
}

## adjust expression for hidden variations for EDA plots
get.adj <- function(expr.data=NULL, design=NULL,  meta.data=NULL)
{
  ## adjusted expression using voom()
  v <- voom(expr.data, design=design)
  fit <- lmFit(v, design=design)
  adj.exp <- v$E -fit$coefficients[, 11:14] %*% t(design[, 11:14])
  adj.exp
}

## get sva: 4 hidden variables
meta.sva <- get.sva(expr.data= normCounts, meta.data=meta_data)
model <- model.matrix(~ 0 + Group + sv1 + sv2 + sv3 +sv4,
                      data=meta.sva)

## adjust for sva
adj.exp <- get.adj(expr.data=normCounts, design=model,  meta.data=meta.sva)
write.table(adj.exp, "Table 2. sva-adjusted.peak.intensity.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

#### Heatmap showing sample distances after adjusting for hidden variation
sampleDists <- dist(t(adj.exp))
pdf("Fig 2.2 Heatmap showing SVA-adjusted sample distances.pdf",
    width = 6, height = 5)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = vsd$Group)
dev.off()

## PCA using adjusted value
pc2 <- prcomp(t(adj.exp), scale = T,
              center = T, retx = TRUE)
pdf(file = "Fig 3.2. PCA.plot.SVA.adjusted.samples.pdf",
    width = 12.0, height = 10)
autoplot(pc2, data = meta_data,
         colour = 'Group', size = 4)
dev.off()

## SVA is much better
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count),
                              colData = meta.sva,
                              design = ~0 + Group + sv1 + sv2 + sv3 +sv4)

dds <- estimateSizeFactors(dds)

dds <- estimateDispersions(dds, fitType="parametric",
                           maxit=1000)
normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
write.table(normCounts, "Table 1.1 Normalized.count.table.SVA-adjusted.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

model <- model.matrix(~0 + Group + sv1 + sv2 + sv3 + sv4, 
                      data = meta.sva)
contrast.matrix <- matrix(c(1, -1, rep(0,12),   ## LIL_TP-LIL_V
                            0, 0, 1, -1, rep(0,10),   ## LIP_TP-LIP_V
                            0, 0, 0, 0, 1, -1, rep(0,8),   ## PIL_TP-PIL_V
                            0, 0, 0, 0, 0, 0, 1, -1, rep(0,6),   ## PIP_TP-PIP_V
                            0, 0, 0, 0, 0, 0, 1, 0, -1, rep(0,5), ## PIP_TP-TP
                            0, 0, 0, 0, 0, 0, 0, 1, 0, -1, rep(0,4), ## PIP_V-V
                           -1, 1, 0, 0, 0, 0, 1, -1, rep(0,6),  ## (PIP_TP-PIP_V)-(LIL_TP-LIL_V)
                            0, 0,-1, 1, 0, 0, 1, -1, rep(0,6),  ## (PIP_TP-PIP_V)-(LIP_TP-LIP_V)
                            0, 0, 0, 0,-1, 1, 1, -1, rep(0,6),  ## (PIP_TP-PIP_V)-(PIL_TP-PIL_V)
                            1,-1,-1, 1, 0, 0, 0, 0, rep(0,6),  ## (LIL_TP-LIL_V)-(LIP_TP-LIP_V)
                            1,-1, 0, 0,-1, 1, 0, 0, rep(0,6),  ## (LIL_TP-LIL_V)-(PIL_TP-PIL_V)
                            0, 0, 0, 0, 0, 0, 1, -1, -1, 1, rep(0,4)  ## (PIP_TP-PIP_V)-(TP-V)
                            
), nrow = 12, byrow = TRUE)

rownames(contrast.matrix) <- c("LIL_TP-LIL_V",
                               "LIP_TP-LIP_V",
                               "PIL_TP-PIL_V",
                               "PIP_TP-PIP_V",
                               "PIP_TP-TP",
                               "PIP_V-V",
                               "(PIP_TP-PIP_V)-(LIL_TP-LIL_V)",
                               "(PIP_TP-PIP_V)-(LIP_TP-LIP_V)",
                               "(PIP_TP-PIP_V)-(PIL_TP-PIL_V)",
                               "(LIL_TP-LIL_V)-(LIP_TP-LIP_V)",
                               "(LIL_TP-LIL_V)-(PIL_TP-PIL_V)",
                               "(PIP_TP-PIP_V)-(TP-V)")

dds <- nbinomWaldTest(dds, modelMatrix = NULL,
                      betaPrior=FALSE,
                      maxit = 50000, 
                      useOptim = TRUE, 
                      quiet = FALSE, 
                      useT = FALSE, useQR = TRUE)

output_DESeq <- function(i, dds, contrast.matrix, 
                         threshold, shrink)
{
  res <- results(dds, alpha = 0.01,
                 contrast=contrast.matrix[i,], 
                 lfcThreshold=log2(threshold),
                 format = "DataFrame",
                 altHypothesis="greaterAbs") 
  if (shrink)
  {
    res <- lfcShrink(dds,
                     contrast = contrast.matrix[i,],
                     res = res,
                     format = "DataFrame",
                     lfcThreshold = log2(threshold),
                     type = "ashr")
  }
  res <- data.frame(res)
  counts <- normCounts
  
  res <- merge(res, counts, by.x = "row.names", 
               by.y = "row.names", all.x = TRUE)
  colnames(res)[1] <- "Peak_id"
  res 
}

## apply raw FLC
DESeq_out <- lapply(1:nrow(contrast.matrix), 
                    output_DESeq, dds = dds, 
                    contrast.matrix = contrast.matrix, 
                    threshold = 1, shrink = FALSE)

WriteXLS(x = DESeq_out, 
         ExcelFileName = "Table 3.2.DESeq.differential.peaks.xlsx",
         row.names = FALSE, 
         SheetNames = rownames(contrast.matrix))

## apply shrunken LFC
DESeq_out_2 <- lapply(1:nrow(contrast.matrix), 
                      output_DESeq, dds = dds, 
                      contrast.matrix = contrast.matrix, 
                      threshold = 1, shrink = TRUE)

WriteXLS(x = DESeq_out_2, 
         ExcelFileName = "Table 3.3.DESeq.differential.peaks.shrunken.xlsx",
         row.names = FALSE, 
         SheetNames = rownames(contrast.matrix))

save.image(file = "H3K27me3.CUT-Tag.DESeq2.RData")

pdf(file = "Fig 7.Volcano plots showing DAR shrunken LFC.pdf", 
    height = 6, width = 6)
null <- lapply(seq_along(DESeq_out_2), function(.x)
{
  res <- DESeq_out_2[[.x]]
  res <- res[!is.na(res$padj), ]
  res$padj <- ifelse(res$padj == 0, min(res$padj[res$padj != 0]), res$padj)
  p <- EnhancedVolcano(res,
                       lab = NA,
                       x = 'log2FoldChange',
                       y = 'padj',
                       ylab = bquote(~-Log[10]~italic(FDR)),
                       subtitle = NULL,
                       caption = paste0(sum(abs(res$log2FoldChange) >= 1 & 
                                              res$padj < 0.05), " of ", 
                                        nrow(res)),
                       axisLabSize = 10,
                       titleLabSize = 12,
                       legendLabSize = 9,
                       title = rownames(contrast.matrix)[.x],
                       pCutoff = 0.05,
                       FCcutoff = 1.0,
                       pointSize = 1.0,
                       border = "full",
                       labSize = 3.0) + theme(plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust = 0.5)) 
  print(p)
})
dev.off()


pdf(file = "Fig 7.Volcano plots showing DAR raw LFC.pdf", 
    height = 6, width = 6)
null <- lapply(seq_along(DESeq_out), function(.x)
{
  res <- DESeq_out[[.x]]
  res <- res[!is.na(res$padj), ]
  res$padj <- ifelse(res$padj == 0, min(res$padj[res$padj != 0]), res$padj)
  p <- EnhancedVolcano(res,
                       lab = NA,
                       x = 'log2FoldChange',
                       y = 'padj',
                       ylab = bquote(~-Log[10]~italic(FDR)),
                       subtitle = NULL,
                       caption = paste0(sum(abs(res$log2FoldChange) >= 1 & 
                                              res$padj < 0.05), " of ", 
                                        nrow(res)),
                       axisLabSize = 10,
                       titleLabSize = 12,
                       legendLabSize = 9,
                       title = rownames(contrast.matrix)[.x],
                       pCutoff = 0.05,
                       FCcutoff = 1.0,
                       pointSize = 1.0,
                       border = "full",
                       labSize = 3.0) + theme(plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust = 0.5)) 
  print(p)
})
dev.off()

## bed file for peaks
# peaks <- strsplit(DESeq_out_2[[1]][, 1], ":")
# peaks <- as.data.frame(do.call("rbind", peaks))
# peaks[, 2] <- as.integer(peaks[, 2]) - 1
# peaks[, 3] <- as.integer(peaks[, 3])
# peaks <- peaks[order(peaks$V1, peaks$V2), ]
# write.table(peaks, file = "H3K27me3.48400.peaks.3col.bed",
#             sep = "\t", quote = FALSE,
#             row.names = FALSE, col.names = FALSE)
# 
# peaks$V4 <- paste(peaks[,1], peaks[,2], peaks[,3], sep = ":")
# peaks$V5 <- "1000"
# peaks$V6 <- "*"
# write.table(peaks, file = "H3K27me3.48400.peaks.6col.bed",
#             sep = "\t", quote = FALSE,
#             row.names = FALSE, col.names = FALSE)


## add annotation to Differential peaks
peaks_ann <- read.delim("Annotated.H3K27me3.q.0.05.consensus.peaks.txt",
                        header = TRUE, as.is = TRUE, sep = "\t")
#peaks_ann[,2] <- peaks_ann[,2] -2
#peaks_ann[ ,3] <- peaks_ann[ ,3] -2

# write.table(peaks_ann, 
#             file = "Annotated.H3K27me3.q.0.05.consensus.peaks.txt",
#             sep = "\t",
#             row.names = FALSE, col.names = TRUE, quote = TRUE)
rownames(peaks_ann) <- paste(peaks_ann[ ,1], peaks_ann[, 2], peaks_ann[ ,3],
                             sep = ":")

DESeq_out <- lapply(DESeq_out, function(.x){
  merge(.x[!is.na(.x$padj) & .x$padj < 0.1, 1:6], peaks_ann[, c(7,15:17)], by.x = "Peak_id",
        by.y = "row.names", all.x = TRUE)
})
                    

WriteXLS(x = DESeq_out, 
         ExcelFileName = "Table 3.2.Annotated.DESeq.differential.peaks.xlsx",
         row.names = FALSE, 
         SheetNames = rownames(contrast.matrix))

DESeq_out_2 <- lapply(DESeq_out_2, function(.x){
  merge(.x[!is.na(.x$padj) & .x$padj < 0.05 & abs(.x$log2FoldChange) >= 1, 1:6], peaks_ann[, c(7, 13, 15:17)], by.x = "Peak_id",
         by.y = "row.names", all.x = TRUE)
})


WriteXLS(x = DESeq_out_2, 
         ExcelFileName = "Table 3.3.Annotated.DESeq.differential.peaks.shrunken.xlsx",
         row.names = FALSE, 
         SheetNames = rownames(contrast.matrix))


## number of up and down peaks
load("DESeq2-normalized.new.analysis/H3K27me3.CUT-Tag.DESeq2.RData")

null <- lapply(DESeq_out_2, function(.x){
    up <- .x[.x$log2FoldChange >= 1 & !is.na(.x$padj) &
               .x$padj <= 0.05, ]
    down <- .x[.x$log2FoldChange <= -1 & !is.na(.x$padj) &
                 .x$padj <= 0.05, ]
    c(up=nrow(up), down = nrow(down))
})

up_down_DARs <- do.call(rbind, null)
rownames(up_down_DARs) <- rownames(contrast.matrix)




## output peaks for GREAT analysis
up_down_regions <- lapply(DESeq_out_2, function(.x){
  up <- unique(.x[.x$log2FoldChange >=1, c(1, 8,10)])
  down <- unique(.x[.x$log2FoldChange <= -1, c(1, 8,10)])
  list(up = up, down = down)
})

save(up_down_regions, file = "004.up.down.DAR.list.for.GREAT.RData")

## output gene lists for ORA

up_down_genes <- lapply(DESeq_out_2, function(.x){
  up <- unique(.x[.x$log2FoldChange >=1, c(8,10)])
  down <- unique(.x[.x$log2FoldChange <= -1, c(8,10)])
  list(up = up, down = down)
})

save(up_down_genes, file = "004.up.down.DAR.associated.genes.list.for.ORA.RData")