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

setwd("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\differential_accessibility")

## DESeq2 analysis
count <- read.delim("q005.all.samples.merged.peaks.final.refmt.count.table", 
                    header = TRUE, as.is =TRUE, check.names = FALSE)

colnames(count) <- gsub("_\\d+$", "", colnames(count), perl = TRUE)

## myh9 gene is very open
summary(count$Length)
# ex_peakd <- count[count$Length > 10000, ]
pdf(file = "Fig 0. Width distr. Consensus peaks.pdf", 
    width = 5, height = 3)
ggplot(count, aes(x = Length)) +
    geom_density(fill = "#e9ecef", alpha = 0) +
    scale_x_continuous(trans = "log10") + 
    xlab("Consensus Peak Width (bp)")
dev.off()

count <- count[, -c(1,2,4:8)]
rownames(count) <- count[, 1]
count <- count[, -1]

meta_data <- read.delim("../2. SampleKey_ATAC seq samples_Marcus.txt",
                        header = TRUE)

meta_data$SampleID <-  gsub("\\s+","_", 
                            meta_data$Customer.Sample.ID, perl = TRUE)
meta_data <- within(meta_data, {
    TumorType <- as.factor(gsub("(.+?)\\s+.+","\\1", 
                      meta_data$Customer.Sample.ID, perl = TRUE))
    Treatment <- as.factor(gsub(".+?\\s+(.+?)[0-9]","\\1", 
                      meta_data$Customer.Sample.ID, perl = TRUE))
    Group <- as.factor(gsub("\\d+", "", meta_data$SampleID, 
                            perl = TRUE))
})

## test for interaction
meta_data <- within(meta_data, {
    Environment <-  as.factor(gsub("[PL]I([PL])", "\\1", 
                                   meta_data$TumorType, perl = TRUE))
    TumorType <- as.factor(gsub("([PL])I[PL]", "\\1", 
                                meta_data$TumorType, perl = TRUE))
})

rownames(meta_data) <- meta_data$SampleID
meta_data <- meta_data[!rownames(meta_data) %in% "PIL_C5", ]


## plot raw count per sample
pdf("Fig 1. Number of reads assigened to consensus peaks.pdf", 
    width = 6, height = 4)
barplot(colSums(count)/1000000, width = 0.5, 
        ylim=c(0, 50),
        ylab = "Reads (Million)", 
        col = "gray",
        cex.axis = 0.5,
        cex.names = 0.5,
        names.arg = colnames(count), las=2)
dev.off()

## order the metadata according to colnames of count table
count <- count[, rownames(meta_data)]
all(colnames(count) == rownames(meta_data))

## filter count table to remove the lowly expressed genes
cpms = cpm(count)
keep = rowSums(cpms > 0.5) >= 2   #141,194 peaks
cpms <- log2(cpms[keep, ] + 1)
count <- count[keep, ]

pdf("Fig 2.Raw.log2cpm.before.qsmooth.pdf", width = 5, height = 4)
matdensity(as.matrix(cpms),
           groupFactor = meta_data$Group,
           ylab= "Density",
           xlab = "log2(cpm + 1)",
           main = "before filtering")
legend('topright', legend = levels(meta_data$Group), 
       col = 1:nlevels(meta_data$Group),lty = 1)
dev.off()

## accounting for sequencing depth differences
rm("cpms", "keep")

lib_sizes <- colSums(count)
geometric_mean <- exp(mean(log(lib_sizes)))
size_factors <- lib_sizes / geometric_mean
counts_1 <- sweep(count, 2, size_factors, FUN = "/")
rm("count")

pdf("Fig 3.Exploratory.data.analysis.before.qsmooth.pdf",
    width = 8, height = 4)
op <- par(mar = c(4, 4, 2, 2))
matboxplot(log2(counts_1+1),
           groupFactor = meta_data$Group, 
           main = "before filtering",
           cex = 0.2, ylab = "log2(counts + 1)")
par(op)

matdensity(log2(as.matrix(counts_1 +1)),
           groupFactor = meta_data$Group,
           ylab= "Density",
           xlab = "log2(counts + 1)",
           main = "before filtering")
legend('topright', legend = levels(meta_data$Group), 
       col = 1:nlevels(meta_data$Group),lty = 1)
dev.off()

## test the median difference across groups.
registerDoParallel(cores=2)
qtestPerm <- quantro(log2(as.matrix(counts_1 + 1)), groupFactor = meta_data$Group, B = 1000)
quantroPlot(qtestPerm)  # quantroPvalPerm:  < 0.001 

qs_rnaseq <- qsmooth(object = log2(counts_1 + 1), 
                     group_factor = meta_data$Group, 
                     norm_factors = NULL, window = 0.05)

qsmoothPlotWeights(qs_rnaseq)

pdf("Fig 3. Exploratory.data.analysis.after.qsmooth.pdf", width = 8, height = 4)
matdensity(as.matrix(qsmoothData(qs_rnaseq)),
           groupFactor = meta_data$Group,
           ylab= "Density",
           main = "after qsmooth normalization",
           xlab = "log2(counts + 1)")
legend('topright', legend = levels(meta_data$Group),
       col = 1:nlevels(meta_data$Group), lty = 1, cex = 0.5)

op <- par(mar=c(4, 4, 2, 2))
matboxplot(as.matrix(qsmoothData(qs_rnaseq)),
           groupFactor = meta_data$Group, 
           main = "after qsmooth normalization",
           cex = 0.3, ylab = "log2(counts + 1")
par(op)

dev.off()

counts_qs  <- round(2^(as.matrix(qsmoothData(qs_rnaseq))) - 1)

write.table(counts_qs , 
            "Table 1. Qsmooth-Normalized.count.table.txt", 
            sep = "\t",
            quote = FALSE, row.names = TRUE)

write.table(as.matrix(qsmoothData(qs_rnaseq)), 
            "Table 2. Qsmooth-Normalized.log2count.table.txt", 
            sep = "\t",
            quote = FALSE, row.names = TRUE)



rm("counts_1", "qs_rnaseq")



#### sample similarity
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts_qs),
                              colData = meta_data,
                              design = ~0 + Group)

"sizeFactors"(dds) <- rep(1, 31)

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
pdf("Fig 4. Qsmooth heatmap showing original sample distances.pdf", width = 6, height = 5)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = vsd$SampleID)
dev.off()


## PCA using adjusted value
pc <- prcomp(t(assay(vsd)), scale = T, center = T, retx = TRUE)
pdf(file = "Fig 5. PCA.plot.after.squantile.pdf", width = 4.0, height = 2.5)
autoplot(pc, data = meta_data, colour = 'Group', size = 2)
dev.off()

# effect model
# design_2 <- model.matrix(~ 0+ TumorType*Treatment*Environment, 
#                         data = meta_data)

design_1 <- model.matrix(~ 0 + Group, 
                         data = meta_data)
contrast.matrix <- matrix(c(0, 0, 0, 0, -1, 1, 1, -1, #PIP vs. PIL
                            -1, 1, 0, 0, 0, 0, 1, -1, #PIP vs. LIL
                            0, 0, -1, 1, 0, 0, 1, -1, #PIP vs. LIP
                            -1, 1, 0, 0, 1, -1, 0, 0, #PIL vs. LIL
                            0, 0, -1, 1, 1, -1, 0, 0, #PIL vs. LIP
                            1, -1, -1, 1, 0, 0, 0,0), #LIL vs. LIP
                          nrow = 6, byrow = TRUE)
rownames(contrast.matrix) <- c("PIP-PIL", "PIP-LIL",
                               "PIP-LIP", "PIL-LIL",
                               "PIL-LIP", "LIL-LIP")
                              

## interaction
# contrast.matrix <- 
#     matrix(c(0, 0, 0, 0, 1, 0, 0, 0, ## TumorXTrt
#              0, 0, 0, 0, 0, 1, 0, 0, ## TumorXEnv
#              0, 0, 0, 0, 0, 0, 1, 0, ## TrtXEnv
#              0, 0, 0, 0, 0, 0, 0, 1,  ## TumorXTrtXEnv
#              
#              0, 0, -1, 0, 0, 0, 0, 0,  ## LIL_C - LIL_V
#              0, 0, -1, 0, 0, 0,-1, 0, ## LIP_C - LIP_V
#              0, 0, -1, 0,-1, 0, 0, 0,  ## PIL_C - PIL_V
#              0, 0, -1, 0,-1, 0,-1,-1, ## PIP_C - PIP_V
#              
#              0, -1, 0, 0, -1, 0, 0, 0,  ## LIL_V - PIL_V
#              0, -1, 0, 0, -1,-1, 0,-1, ## LIP_V - PIP_V
#              0, 0, 0, 1, 0, 1, 1, 1,  ## PIP_V - PIL_V
#              0, 0, 0, 1, 0, 0, 1, 0, ## LIP_V - LIL_V
#              
#              0,-1, 0, 0, 0, 0, 0, 0, ## LIL_C - PIL_C
#              0,-1, 0, 0, 0,-1, 0, 0, ## LIP_C - PIP_C
#              0, 0, 0, 1, 0, 1, 0, 0,  ## PIP_C - PIL_C
#              0, 0, 0, 1, 0, 0, 0, 0  ## LIP_C - LIL_C
#     ),  
#     nrow = 16, byrow = TRUE)

# rownames(contrast.matrix) <- c("Tumor X Trt", "Tumor X Env",
#                                "Trt X Env", "Tumor X Trt X Env",
#                                
#                                "LIL_C - LIL_V", "LIP_C - LIP_V",
#                                "PIL_C - PIL_V", "PIP_C - PIP_V",
#                                
#                                "LIL_V - PIL_V", "LIP_V - PIP_V",
#                                "PIP_V - PIL_V", "LIP_V - LIL_V",
#                                
#                                "LIL_C - PIL_C", "LIP_C - PIP_C",
#                                "PIP_C - PIL_C", "LIP_C - LIL_C")

rm("pc", "op", "qtestPerm", "res", "res_shrink", "small.count", "vsd")
dds <- nbinomWaldTest(dds, modelMatrix = NULL,
                      betaPrior=FALSE,
                      maxit = 50000, useOptim = TRUE, 
                      quiet = FALSE, 
                      useT = FALSE, useQR = TRUE)
## peak informatio

peaks <- read.delim("q005.all.samples.merged.peaks.final.refmt.count.table", 
           header = TRUE, as.is =TRUE, check.names = FALSE)[, c(1:6,8)]

output_DESeq <- function(i, dds, contrast.matrix, 
                         threshold, shrink)
{
    print(i)
    

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
    res <- as.data.frame(res)
    res <- res[!is.na(res$padj) & res$padj < 0.05 
               & abs(res$log2FoldChange) >= log2(1.5), ]
    res <- merge(peaks, res, by.x = "Geneid", 
                 by.y = "row.names", all.y = TRUE)
    
    res 
}

## apply raw FLC
DESeq_out <- lapply(1:nrow(contrast.matrix), 
                    output_DESeq, dds = dds, 
                    contrast.matrix = contrast.matrix, 
                    threshold = 1, shrink = FALSE)


file <- ".DAR.raw.LFC.xlsx"

for (i in seq_along(DESeq_out))
{
    WriteXLS(x = DESeq_out[[i]], 
             ExcelFileName = paste0(rownames(contrast.matrix)[i],file),
             row.names = FALSE, 
             SheetNames = rownames(contrast.matrix)[i])
}
save.image(file = "ATAC-seq.DAA.raw.LFC.RData")


## apply shrunken LFC
DESeq_out_2 <- lapply(1:nrow(contrast.matrix), 
                    output_DESeq, dds = dds, 
                    contrast.matrix = contrast.matrix, 
                    threshold = 1, shrink = TRUE)

file <- ".DAR.shrunken.LFC.xlsx"

for (i in seq_along(DESeq_out_2))
{
    WriteXLS(x = DESeq_out_2[[i]], 
             ExcelFileName = paste0(rownames(contrast.matrix)[i],file),
             row.names = FALSE, 
             SheetNames = rownames(contrast.matrix)[i])
}

save(DESeq_out_2, file = "ATAC-seq.DAA.shrunken.LFC.RData")

## volcano plot showing DAR
load("ATAC-seq.DAA.raw.LFC.RData")
#load("ATAC-seq.DAA.shrunken.LFC.RData")
pdf(file = "Fig 7.Volcano plots showing DAR shrunken LFC.pdf", height = 8, width = 4)
null <- lapply(seq_along(DESeq_out_2), function(.x)
{
    res <- DESeq_out_2[[.x]]
    res$padj <- ifelse(res$padj == 0, min(res$padj[res$padj != 0]), res$padj)
    p <- EnhancedVolcano(res,
                    lab = NA,
                    x = 'log2FoldChange',
                    y = 'padj',
                    ylab = bquote(~-Log[10]~italic(FDR)),
                    xlim = c(min(range(res$log2FoldChange))-0.5, 
                             max(range(res$log2FoldChange)) + 0.5),
                    subtitle = NULL,
                    caption = paste0(sum(abs(res$log2FoldChange) >= 1 & 
                                             res$padj < 0.05), " of ", 
                                     nrow(res)),
                    axisLabSize = 10,
                    titleLabSize = 12,
                    legendLabSize = 9,
                    ylim = c(0, max(-log10(res$padj), 
                                    na.rm = TRUE) + 1),
                    title = rownames(contrast.matrix)[.x],
                    pCutoff = 0.05,
                    FCcutoff = 1.0,
                    pointSize = 1.0,
                    border = "full",
                    labSize = 3.0) + theme(plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust = 0.5)) 
    print(p)
})
dev.off()


install.packages("heatmaply")
library("heatmaply")
library(RColorBrewer)
library("pheatmap")
library(readxl)

## heatmap showing senescence genes
setwd("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\differential_accessibility")
load("ATAC-seq.DAA.raw.LFC.RData")
vsd <- vst(dds, blind = TRUE)
vsd_val <- assay(vsd)

setwd("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\ChIPseeker.out")

peaks_count <- read_xlsx("Senesce.geneset.annotated.ATAC-seq.consensus.peaks.normalized.counts.xlsx",
                         sheet = 1, guess_max = 21474836)

senescence_genes <- read.delim(file.path("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020", 
                                          "All.mouse.senescence.geneset.refmt.txt"), 
                               header = TRUE, check.names =  FALSE)
peaks_count <- peaks_count[peaks_count$geneId %in% rownames(senescence_genes), ]

peaks_count <- peaks_count[, 1:49]
peaks_count <- peaks_count[, c(1:3, 7, 13,16,18:49)]
peaks_count <- peaks_count[!is.na(peaks_count$PIP_V1), ]
peaks_count <- as.data.frame(peaks_count)

vsd_val<- vsd_val[rownames(vsd_val) %in% peaks_count$Geneid, ]
vsd_val <- as.data.frame(vsd_val)
peaks_count <- merge(peaks_count[, 1:7], vsd_val, by.x = "Geneid", 
      by.y = 'row.names', all = TRUE)

rownames(peaks_count) <- with(peaks_count, 
                              paste0(Symbol,"|",seqnames, ":", start, "-", end))

peaks_count <- peaks_count[, c(6, 8:38)]

peaks_count <- peaks_count[, order(colnames(peaks_count))]


## plot heatmap for each gene sets

col_annot <- data.frame(group = gsub("\\d+", "", colnames(peaks_count)[-1]))
rownames(col_annot) <- colnames(peaks_count)[-1]
col_annot$group <- as.factor(col_annot$group)
annotation_colors <- list(group = c("LIL_C" = "#911eb4", "LIL_V" = "#dcbeff",
                                    "LIP_C" = "#3cb44b", "LIP_V" = "#bfef45",
                                    "PIL_C" = "#4363d8", "PIL_V" = "#42d4f4",
                                    "PIP_C" = "#f58231","PIP_V" = "#ffe119"))

out <- "C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\Heatmap.2" 
if(!dir.exists(out))
{
    dir.create(out, recursive = TRUE)
}


for (i in seq_along(senescence_genes[, -27]))
{
    geneset <-senescence_genes[senescence_genes[, i] != "", i, drop = FALSE]
    dim(geneset)
    temp_count <- peaks_count[peaks_count$geneId %in% rownames(geneset), -1]
    temp_count[] <- as.numeric(as.matrix(temp_count))
    # temp_mat <- t(scale(t(temp_count)))
    # 
    # mat <- temp_count
    # mat[] <- lapply(colnames(mat), function(colname) {
    #     paste0(mat[, colname], ", ", colname)
    # })
    print(colnames(senescence_genes)[i])
    print(dim(temp_count))
    
    out_file <- file.path(out, 
                          paste0(colnames(senescence_genes)[i], 
                                 ".pdf"))
    
  # heatmaply(
  #       as.matrix(temp_mat),
  #       seriate = "none",
  #       distfun = "spearman",
  #       hclust_method = "ward.D2",
  #       column_text_angle = 90,
  #       colors = rev(colorRampPalette(brewer.pal(3, "RdBu"))(256)),
  #       row_dend_left = TRUE,
  #       labRow = rownames(temp_mat),
  #       fontsize_row = 0.5,
  #       col_side_colors = col_annot,
  #       custom_hovertext = mat,
  #       file = out_file
  #   ) %>% orca(width=100, height = 1200)
  #   browseURL(out_file)
    
    pdf(out_file,
        height = 6*nrow(temp_count)*0.0139 + 0.5 , width = 6)
    heat <- pheatmap(as.matrix(temp_count),
                     cluster_rows = TRUE,
                     cluster_cols = FALSE,
                     clustering_method = "ward.D2",
                     clustering_distance_rows = "correlation",
                     scale="row",
                     cellheight = 6,
                     border_color = NA,
                     treeheight_col = 50,
                     treeheight_row = 0,
                     cellwidth = 6,
                     fontsize_row = 5,
                     annotation_col = col_annot,
                     gaps_col = c(5,8,12,15,18,22,27),
                     annotation_colors = annotation_colors,
                     color =  rev(colorRampPalette(brewer.pal(3, "RdBu"))(256)),
                     show_colnames = FALSE,
                     show_rownames = TRUE)
    print(heat)
    dev.off()
}


## plot fold changes between pairwise comparison: PIL, PIP, LIL, LIP

setwd("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\differential_accessibility")
load("ATAC-seq.DAA.raw.LFC.RData")

setwd("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\ChIPseeker.out")

peaks_count <- read_xlsx("Senesce.geneset.annotated.ATAC-seq.consensus.peaks.normalized.counts.xlsx",
                         sheet = 1, guess_max = 21474836)
senescence_genes <- read.delim(file.path("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020", 
                                         "All.mouse.senescence.geneset.refmt.txt"), 
                               header = TRUE, check.names =  FALSE)
peaks_count <- peaks_count[peaks_count$geneId %in% rownames(senescence_genes), ]
peaks_count <- peaks_count[, c(1:3, 7, 13,16,18:19, 50:75)]
peaks_count <- peaks_count[!is.na(peaks_count$PIP_V1), ]
peaks_count <- as.data.frame(peaks_count)

output_DESeq <- function(i, dds, contrast.matrix, 
                         threshold, shrink)
{
  print(i)
  
  
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
  res <- as.data.frame(res)
  # res <- res[!is.na(res$padj) & res$padj < 0.05 
  #            & abs(res$log2FoldChange) >= log2(1.5), ]
  res <- merge(peaks, res, by.x = "Geneid", 
               by.y = "row.names", all.y = TRUE)
  
  res 
}

## apply raw FLC
DESeq_out_rawLFC <- lapply(1:nrow(contrast.matrix), 
                    output_DESeq, dds = dds, 
                    contrast.matrix = contrast.matrix, 
                    threshold = 1, shrink = FALSE)

DESeq_out_shrunkenLFC <- lapply(1:nrow(contrast.matrix), 
                           output_DESeq, dds = dds, 
                           contrast.matrix = contrast.matrix, 
                           threshold = 1, shrink = TRUE)
DESeq_out_rawLFC <- DESeq_out_shrunkenLFC

geneset_changes <- lapply(DESeq_out_rawLFC, function(.x){
   temp <- merge(.x[, c("Geneid", "baseMean", "log2FoldChange", "padj")], 
                  peaks_count,
         by = "Geneid")
    temp <- temp[order(temp$Geneid), ]
    rownames(temp) <- with(temp, 
                           paste0(Symbol,"|",seqnames, ":", start, "-", end))
    temp
})

geneset_logFC <- lapply(geneset_changes, function(.x){
  .x[, c(3,12:37)]
})
geneset_logFC_cmb <- cbind(geneset_logFC[[1]][, 1, drop = FALSE], 
                       geneset_logFC[[2]][, 1, drop = FALSE],
                       geneset_logFC[[3]][, 1, drop = FALSE],
                       geneset_logFC[[4]][, 1, drop = FALSE],
                       geneset_logFC[[5]][, 1, drop = FALSE],
                       geneset_logFC[[6]])
names(geneset_logFC_cmb)[1:6] <- rownames(contrast.matrix)
peaks_count <- geneset_logFC_cmb[, order(colnames(geneset_logFC_cmb)[1:6])]
senescence_genes <-geneset_logFC_cmb[, 7:32] 

## plot heatmap for each gene sets

col_annot <- data.frame(group = colnames(peaks_count))
rownames(col_annot) <- colnames(peaks_count)
col_annot$group <- as.factor(col_annot$group)
annotation_colors <- list(group = c("LIL-LIP" = "#911eb4", "PIL-LIL" = "#dcbeff",
                                    "PIL-LIP" = "#3cb44b", "PIP-LIL" = "#bfef45",
                                    "PIP-LIP" = "#4363d8", "PIP-PIL" = "#42d4f4"))

out <- "C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\Heatmap_shrunken.logFC" 
if(!dir.exists(out))
{
  dir.create(out, recursive = TRUE)
}


for (i in seq_along(senescence_genes))
{
  temp_count <- peaks_count[!is.na(senescence_genes[, i]), ]
  print(colnames(senescence_genes)[i])
  print(dim(temp_count))
  
  out_file <- file.path(out, 
                        paste0(colnames(senescence_genes)[i], 
                               ".pdf"))
  
  pdf(out_file,
      height = 6*nrow(temp_count)*0.0139 + 0.5 , width = 4)
  heat <- pheatmap(as.matrix(temp_count),
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   clustering_method = "ward.D2",
                   clustering_distance_rows = "correlation",
                   scale="none",
                   cellheight = 6,
                   border_color = NA,
                   treeheight_col = 50,
                   treeheight_row = 0,
                   cellwidth = 6,
                   fontsize_row = 5,
                   annotation_col = col_annot,
                   annotation_colors = annotation_colors,
                   color =  rev(colorRampPalette(brewer.pal(3, "RdBu"))(256)),
                   show_colnames = FALSE,
                   show_rownames = TRUE)
  print(heat)
  dev.off()
}


## plot log2FC changes as scatter plots between treatment and control
setwd("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\differential_accessibility\\31 samples using effect model")
load("ATAC-seq.DAA.raw.LFC.RData")
# load("ATAC-seq.DAA.shrunken.LFC.RData")

output_DESeq <- function(i, dds, contrast.matrix, 
                         threshold, shrink)
{
  print(i)
  
  
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
  res <- as.data.frame(res)
  # res <- res[!is.na(res$padj) & res$padj < 0.05 
  #            & abs(res$log2FoldChange) >= log2(1.5), ]
  res <- merge(peaks, res, by.x = "Geneid", 
               by.y = "row.names", all.y = TRUE)
  
  res 
}

## apply raw FLC
DESeq_out_rawLFC <- lapply(1:nrow(contrast.matrix), 
                           output_DESeq, dds = dds, 
                           contrast.matrix = contrast.matrix, 
                           threshold = 1, shrink = FALSE)
DESeq_out_rawLFC <- lapply(5:8, 
                           output_DESeq, dds = dds, 
                           contrast.matrix = contrast.matrix, 
                           threshold = 1, shrink = TRUE)

names(DESeq_out_rawLFC) <- rownames(contrast.matrix)[5:8]


library(ggplot2)
library(ggExtra)
library("gridExtra")

plist <- mapply(function(.x, .y){
  .x$log10BaseMean <- log10(.x$baseMean)
  .x$Status <- as.factor(ifelse(!is.na(.x$padj) & 
                                  .x$padj <= 0.05 & 
                               abs(.x$log2FoldChange) >= log2(1.5), 
                             "DAR", "nonDAR"))
  # classic plot :
  p <- ggplot(.x, aes(x= log10(baseMean), 
                      y=log2FoldChange, 
                      color = Status)) +
       geom_point(alpha = 0.2, shape = 20, size = 2) +
       geom_hline(yintercept = 0, linetype = "dashed") +
    theme(legend.position="bottom") + ggtitle(.y) + 
    theme(plot.title = element_text(hjust = 0.5, siz = 12))
  
  # Show only marginal plot for x axis
  p1 <- ggMarginal(p, margins = 'y', 
                   type = "density",
                   color="purple", 
                   size = 3)
  p1
}, DESeq_out_rawLFC, names(DESeq_out_rawLFC), SIMPLIFY = FALSE)

#grid.arrange(grobs = plist, ncol = 2) ## display plot
ggsave(file = file.path(out, "Pairwise shrunken log2FC scatterplot.tiff"),
       arrangeGrob(grobs = plist, ncol = 2), 
       height = 20, width = 16,
       units = "in", dpi = 600)


## shrunken LFC
## find genes altered after combo treatment in PIP and LIP compared to vehicle

DESeq_out_rawLFC <- lapply(DESeq_out_rawLFC, function(.x){
  .x <- .x[, c(1, 8, 9, 12)]
  .x
})

four_cond_DAR <- do.call(cbind, DESeq_out_rawLFC)

four_cond_DAR <- four_cond_DAR[, -c(5, 6, 9, 10, 13, 14)]
rownames(four_cond_DAR) <- four_cond_DAR[, 1]
four_cond_DAR <- four_cond_DAR[, -1]

four_cond_DAR_sig <- four_cond_DAR[((!is.na(four_cond_DAR[, 3]) & four_cond_DAR[, 3] <= 0.05) & 
                                 abs(four_cond_DAR[, 2]) >= log2(1.5)) |
                                ((!is.na(four_cond_DAR[, 5]) & four_cond_DAR[, 5] <= 0.05) & 
                                   abs(four_cond_DAR[, 4]) >= log2(1.5)) |
                                ((!is.na(four_cond_DAR[, 7]) & four_cond_DAR[, 7] <= 0.05) &
                                  abs(four_cond_DAR[, 6]) >= log2(1.5)) |
                                ((!is.na(four_cond_DAR[, 9]) & four_cond_DAR[, 9] <= 0.05) &
                                   abs(four_cond_DAR[, 8]) >= log2(1.5)), ]


four_cond_DAR_sig <- four_cond_DAR_sig[, c(1, 8,9,4,5,2,3,6,7)]

setwd("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\ChIPseeker.out")

peaks_count <- read_xlsx("Senesce.geneset.annotated.ATAC-seq.consensus.peaks.normalized.counts.xlsx",
                         sheet = 1, guess_max = 21474836)
senescence_genes <- read.delim(file.path("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020", 
                                         "All.mouse.senescence.geneset.refmt.txt"), 
                               header = TRUE, check.names =  FALSE)
peaks_count <- peaks_count[peaks_count$geneId %in% rownames(senescence_genes), ]
peaks_count <- peaks_count[, c(1:3, 7, 13,16,18:19, 50:75)]
peaks_count <- peaks_count[!is.na(peaks_count$PIP_V1), ]
peaks_count <- as.data.frame(peaks_count)


four_cond_DAR_sig_senescence <- merge(four_cond_DAR_sig, peaks_count, 
                                      by.x = "row.names", by.y = "Geneid")
rownames(four_cond_DAR_sig_senescence) <- with(four_cond_DAR_sig_senescence,
                                               paste0(Symbol,"|",seqnames, ":", start, "-", end))

peaks_count <- four_cond_DAR_sig_senescence[, c(3,5,7,9)]
colnames(peaks_count) <- gsub("\\.log2FoldChange", "", colnames(peaks_count), perl = TRUE)

senescence_genes <-four_cond_DAR_sig_senescence[, 18:43] 

## plot heatmap for each gene sets

col_annot <- data.frame(group = colnames(peaks_count))
rownames(col_annot) <- colnames(peaks_count)
col_annot$group <- as.factor(col_annot$group)
annotation_colors <- list(group = c("PIP_C - PIP_V" = "#911eb4", 
                                    "LIP_C - LIP_V" = "#dcbeff",
                                    "LIL_C - LIL_V" = "#3cb44b",
                                    "PIL_C - PIL_V" = "#bfef45"))

out <- "C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\Heatmap_shrunken.logFC.C.vs.V" 
if(!dir.exists(out))
{
  dir.create(out, recursive = TRUE)
}


for (i in seq_along(senescence_genes))
{
  temp_count <- peaks_count[!is.na(senescence_genes[, i]), ]
  print(colnames(senescence_genes)[i])
  print(dim(temp_count))
  
  out_file <- file.path(out, 
                        paste0(colnames(senescence_genes)[i], 
                               ".pdf"))
  
  pdf(out_file,
      height = 6*nrow(temp_count)*0.0139 + 0.5 , width = 4)
  heat <- pheatmap(as.matrix(temp_count),
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   clustering_method = "ward.D2",
                   clustering_distance_rows = "correlation",
                   scale="none",
                   cellheight = 6,
                   border_color = NA,
                   treeheight_col = 50,
                   treeheight_row = 0,
                   cellwidth = 6,
                   fontsize_row = 5,
                   annotation_col = col_annot,
                   annotation_colors = annotation_colors,
                   color =  rev(colorRampPalette(brewer.pal(3, "RdBu"))(256)),
                   show_colnames = FALSE,
                   show_rownames = TRUE)
  print(heat)
  dev.off()
}

##
setwd("C:\\Users\\liuh\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\ChIPseeker.out")

peaks_count <- read_xlsx("Senesce.geneset.annotated.ATAC-seq.consensus.peaks.normalized.counts.xlsx",
                         sheet = 1, guess_max = 21474836)
peaks_count <- peaks_count[, c(1:3, 7, 13,16,18:19, 50:75)]
peaks_count <- peaks_count[!is.na(peaks_count$PIP_V1), ]


four_cond_DAR_sig_o <- merge(four_cond_DAR_sig, peaks_count[, 1:7], by.x = "row.names",
                             by.y = "Geneid")
rownames(four_cond_DAR_sig_o) <- with(four_cond_DAR_sig_o,
                                      paste0(Symbol,"|",seqnames, ":", start, "-", end))
four_cond_DAR_sig_n <- four_cond_DAR_sig_o[, c(1, 3, 5, 7, 9)]
four_cond_DAR_sig_n <- four_cond_DAR_sig_n[, -1]

colnames(four_cond_DAR_sig_n) <- gsub("\\.log2FoldChange|\\s+", "",
                                      colnames(four_cond_DAR_sig_n), 
                                      perl = TRUE)
col_annot <- data.frame(group = colnames(four_cond_DAR_sig_n))
rownames(col_annot) <- colnames(peaks_count)
col_annot$group <- as.factor(col_annot$group)
annotation_colors <- list(group = c("PIP_C-PIP_V" = "#911eb4", 
                                    "LIP_C-LIP_V" = "#dcbeff",
                                    "LIL_C-LIL_V" = "#3cb44b",
                                    "PIL_C-PIL_V" = "#bfef45"))

out_file <- file.path(out, "All.DAR.pdf")
temp_count <- four_cond_DAR_sig_n


library(genefilter)
plot(density(rowVars(as.matrix(temp_count))), 
     xlim =c(0, 0.8))

sum(rowVars(as.matrix(temp_count)) < 0.4)
temp_count_o <- temp_count
temp_count <- temp_count_o
temp_count <- temp_count[rowVars(as.matrix(temp_count)) >= 0.4, ]


heat <- pheatmap(as.matrix(temp_count),
                 cluster_rows = TRUE,
                 cluster_cols = FALSE,
                 clustering_method = "ward.D2",
                 clustering_distance_rows = "correlation",
                 scale="none",
                 cellheight = 1,
                 border_color = NA,
                 treeheight_col = 50,
                 treeheight_row = 300,
                 cellwidth = 10,
                 fontsize_row = 0.8,
                 annotation_col = col_annot,
                 annotation_colors = annotation_colors,
                 color =  rev(colorRampPalette(brewer.pal(3, "RdBu"))(256)),
                 show_colnames = FALSE,
                 show_rownames = TRUE)


pdf("all.C.vs.V.var.ls.0.1.pdf",
    height = 1*nrow(temp_count)*0.0139 + 2, width = 12)
print(heat)
dev.off()

## find patterned expression
temp_count <- temp_count_o

pip_lip_dist <- abs(temp_count[,1] - temp_count[,2])
pip_lil_dist <- abs(temp_count[,1] - temp_count[,3])
pip_pil_dist <- abs(temp_count[,1] - temp_count[,4])
lip_lil_dist <- abs(temp_count[,2] - temp_count[,3])
lip_lil_dist <- abs(temp_count[,2] - temp_count[,4])
lil_pil_dist <- abs(temp_count[,3] - temp_count[,4])


temp_count <- temp_count[pip_lip_dist < pip_lil_dist &
                           pip_lip_dist <  pip_pil_dist &
                           pip_lip_dist < lip_lil_dist &
                           pip_lip_dist < lip_lil_dist &
                           lil_pil_dist < pip_lil_dist &
                           lil_pil_dist < pip_pil_dist &
                           lil_pil_dist < lip_lil_dist &
                           lil_pil_dist < lip_lil_dist, ]

heat <- pheatmap(as.matrix(temp_count),
                 cluster_rows = TRUE,
                 cluster_cols = FALSE,
                 clustering_method = "ward.D2",
                 clustering_distance_rows = "correlation",
                 scale="none",
                 cellheight = 2,
                 border_color = NA,
                 treeheight_col = 50,
                 treeheight_row = 300,
                 cellwidth = 10,
                 fontsize_row = 1.8,
                 annotation_col = col_annot,
                 annotation_colors = annotation_colors,
                 color =  rev(colorRampPalette(brewer.pal(3, "RdBu"))(256)),
                 show_colnames = FALSE,
                 show_rownames = TRUE)


pdf("all.C.vs.V.two.groups.pdf",
    height = 2*nrow(temp_count)*0.0139 + 2, width = 8)
print(heat)
dev.off()

save(heat, four_cond_DAR_sig_o, temp_count, annotation_colors, col_annot,
     file = "all.C.vs.V.two.groups.RData")


## 07/22/2021 generate figures for publication
setwd("C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\differential_accessibility\\31 samples using effect model")

library("readxl")
library("tidyverse")
xls <- dir(".", ".xlsx$")
comparison_name <- c("(LIL_C - LIL_V) - (LIP_C - LIP_V)",
                     "(PIP_C - PIP_V) - (LIL_C - LIL_V)",
                     "(PIP_C - PIP_V) - (PIL_C - PIL_V)",
                     "LIL_C - LIL_V",
                     "LIP_C - LIP_V",
                     "PIL_C - PIL_V",
                     "PIP_C - PIP_V")
names(xls) <- comparison_name

DAR_num <- lapply(xls, function(.x){
  dat <- read_xlsx(.x, sheet = 1,  
                   guess_max = 21474836)
  n.up <- dat %>% filter(log2FoldChange >=1) %>% summarise(n = n())
  n.down <- dat %>% filter(log2FoldChange <= -1) %>%
    summarise(n = n())
  df <- data.frame(up = n.up, 
             down = n.down)
  colnames(df) <- c("n.up", "n.down")
  df
})
DAR_num <- do.call("rbind" ,DAR_num)
rownames(DAR_num) <- comparison_name

DAR_num <- DAR_num[c(4:7, 1:3), ]
DAR_num$contrast <- rownames(DAR_num)

library("reshape2")
library("ggplot2")
DAR_num_0 <- melt(DAR_num, variable = "Status",
                value.name = "DAR")
DAR_num_0$contrast <- factor(DAR_num_0$contrast, 
                             levels = c(comparison_name[4:7],
                                        comparison_name[1:3]))
DAR_num_0$Status <- ifelse(DAR_num_0$Status=="n.up", "Gain", "Loss" )

pdf("Fig 0. Number of DARs.pdf", height = 3.5, width = 5)
ggplot(DAR_num_0, 
       aes(x = contrast, y = DAR, 
           group = Status, fill = Status)) +
  geom_bar(stat="identity", position = position_dodge(width = 0.4), width = 0.4) +
  theme(text = element_text(size=8))+
  geom_text(data = DAR_num_0, 
            mapping = aes(x = contrast, y = DAR, label = DAR),
            stat="identity", position = position_dodge(width = 0.4), 
            size = 2, vjust = "center", hjust = "inward") + 
  theme(axis.text = element_text(size = 8), legend.key.size = unit(0.25, 'cm')) +
  xlab("Contrast") +  ylab("DARs") + 
  coord_flip()
dev.off()

## upset

DAR_names<- lapply(xls, function(.x){
  #.x <- xls[1]
  dat <- read_xlsx(.x, sheet = 1,  
                   guess_max = 21474836)
  peaks.up <- dat %>% filter(log2FoldChange >=1) %>% select(Geneid) 
  peaks.down <- dat %>% filter(log2FoldChange <= -1) %>%select(Geneid)
  
  list(peaks.up = as.data.frame(peaks.up)[,1], 
       peaks.down = as.data.frame(peaks.down)[,1])
})

DAR_up_names <- lapply(DAR_names, function(.x){
  .x[[1]]
})

library("UpSetR")
pdf("Fig 1.1 Gain.DAR.intersection.between.conditions.pdf",
     height = 5, width =10)
upset(fromList(DAR_up_names), nsets= 7, order.by = "freq")

dev.off()
DAR_down_names <- lapply(DAR_names, function(.x){
  .x[[2]]
})
pdf("Fig 1.2 Loss.DAR.intersection.between.conditions.pdf",
    height = 5, width =10)
upset(fromList(DAR_down_names), nsets= 7, order.by = "freq")
dev.off()

up_down <- append(DAR_up_names, DAR_down_names)
names(up_down) <- paste(names(up_down), rep(c("Gain", "Loss"), each = 7),sep =":")
pdf("Fig 1.3 Gain and Loss.DAR.intersection.between.conditions.pdf",
    height = 6, width =14)
upset(fromList(up_down), nsets= 14, sets = names(up_down),
      keep.order = TRUE, order.by = "freq", 
      nintersects = 65)
dev.off()

## DAR peaks
DAR_names<- lapply(xls, function(.x){
  dat <- read_xlsx(.x, sheet = 1,  
                   guess_max = 21474836)
  peaks.up <- dat %>% filter(log2FoldChange >=1) %>% select(Geneid) 
  peaks.down <- dat %>% filter(log2FoldChange <= -1) %>%select(Geneid)
  
  list(peaks.up = as.data.frame(peaks.up)[,1], 
       peaks.down = as.data.frame(peaks.down)[,1])
})
DAR_up_names <- lapply(DAR_names, function(.x){
  .x[[1]]
})

DAR_down_names <- lapply(DAR_names, function(.x){
  .x[[2]]
})

output_DAR <- function(outdir, data){
  out <- outdir
  if (!dir.exists(out)){
    dir.create(out)
  }
  
  null <- mapply(function(.x, .y){
    write.table(.x, file= gsub(" ", "", 
                               file.path(out, paste0(.y, ".DAR.txt"))),
                sep = "\t", quote = FALSE, 
                row.names = FALSE,
                col.names = FALSE)
  }, data, names(data), SIMPLIFY = FALSE)
}

output_DAR(outdir = "DAR_UP", data = DAR_up_names)
output_DAR(outdir = "DAR_DOWN", data = DAR_down_names)

