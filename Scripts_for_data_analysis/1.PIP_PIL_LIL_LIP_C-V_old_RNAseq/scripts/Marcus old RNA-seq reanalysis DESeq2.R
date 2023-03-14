install.packages("BiocManager")

pkgs <- c("ggplot2", "sva","DESeq2", "edgeR", "limma",
          "pheatmap", "RColorBrewer", "dplyr", "ggfortify",
          "genefilter", "gplots", "doParallel")
pkg2install <- pkgs[!pkgs%in% rownames(installed.packages())]

BiocManager::install(pkg2install, dependencies = TRUE)

library("ggplot2")
library("sva")
library("DESeq2")
library("edgeR")
library("limma")
library("pheatmap")
library("RColorBrewer")
library("dplyr")
library("ggfortify")
library("genefilter")
library("gplots")
library("doParallel")
library("WriteXLS")
library("readxl")

setwd("C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus\\old.RNA-seq.analysis")

filename <- "featureCounts_Annotated 8-2-21.cnt.csv"

count_PE <- read.delim(filename, sep = ",",
                       header = TRUE, as.is =TRUE, check.names = FALSE)

colnames(count_PE) <- gsub(" ",
                           "_", colnames(count_PE))
rownames(count_PE) <- count_PE[, 1]
count_PE <- count_PE[, 2:31]

sample_name <- colnames(count_PE)
metadata <- data.frame(sample_id = sample_name,  
                       group = gsub("(.+?[CV]).+", "\\1", 
                                    sample_name, perl = TRUE),
                       replicates = gsub(".+?[CV](.+)", "\\1", 
                                         sample_name, perl = TRUE))
metadata$group <- factor(metadata$group)

count_PE <- count_PE[rowSums(count_PE) != 0, ]

## total reads assigned to gene features
total <- data.frame(total = colSums(count_PE))
total <- cbind(metadata, total)
total$sample_id <- paste0(total$group, total$replicates)

pdf("Fig. 1. Total number of reads assigned to gene features.pdf",
    width = 6, height =5)
ggplot(total, aes(sample_id, total, fill = group))+
  geom_col() + 
  theme(axis.text.x = 
          element_text(angle=90,
                       hjust = 0, vjust = 0.5))
dev.off()


## remove lowly expressed genes
cpms = cpm(count_PE)
keep = rowSums(cpms > 1) >= 2
count_PE = count_PE[keep, ]
rm("cpms", "keep")

dim(count_PE)  ## 14779
all(colnames(count_PE) == metadata$sample_id)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_PE),
                              colData = metadata,
                              design = ~0 + group)

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

## one sample normal:SRR1313203.1 is an outlier, so it should be 
## removed and reanalyze it.
pdf("Fig 2. Heatmap showing sample distances.pdf", 
    width = 4.5, height = 4)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = paste0(vsd$group, vsd$replicates))
dev.off()


## PCA using adjusted value
pc <- prcomp(t(assay(vsd)), scale = T, 
             center = T, retx = TRUE)
pdf(file = "Fig 3. PCA.plot.pdf", 
    width = 12.0, height = 10)
autoplot(pc, data = metadata, 
         colour = 'group', size = 4)
dev.off()


normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
write.table(normCounts, "Table 1. Normalized.count.table.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)


model <- model.matrix(~0 + group, data = metadata)
contrast.matrix <- matrix(c(1, -1, 0, 0, rep(0,4),       ## LIL_C-LIL_V
                            0,  0, 1, -1, rep(0,4),      ## LIP_C-LIP_V
                            rep(0, 4), 1, -1, rep(0,2),  ## PIL_C-PIL_V
                            rep(0, 6), 1, -1             ## PIP_C-PIP_V
                            
), nrow = 4, byrow = TRUE)
rownames(contrast.matrix) <- c("LIL_C-LIL_V", 
                               "LIP_C-LIP_V",
                               "PIL_C-PIL_V",
                               "PIP_C-PIP_V")

dds <- nbinomWaldTest(dds, modelMatrix = NULL,
                      betaPrior=FALSE,
                      maxit = 50000, 
                      useOptim = TRUE, 
                      quiet = FALSE, 
                      useT = FALSE, useQR = TRUE)

gene_name <- read.delim("Ensembl.Mus.GCRm38.101.gene.name.txt", 
                        header = FALSE, as.is = TRUE)
colnames(gene_name) <- c("Gene_id", "Symbol", "Biotype")

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
  counts=normCounts
  
  res <- merge(res, counts, by.x = "row.names", 
               by.y = "row.names", all.x = TRUE)
  colnames(res)[1] <- "Gene_id"
  res <- merge(res, gene_name, by= "Gene_id", all.x = TRUE)
  res 
}

## apply raw FLC
DESeq_out <- lapply(1:nrow(contrast.matrix), 
                    output_DESeq, dds = dds, 
                    contrast.matrix = contrast.matrix, 
                    threshold = 1, shrink = FALSE)

WriteXLS(x = DESeq_out, 
         ExcelFileName = "Table 3.2. DESeq.differential expressed genes.xlsx",
         row.names = FALSE, 
         SheetNames = rownames(contrast.matrix))

save.image(file = "old.RNA-seq.DESeq2.RData")


## prepare GSEA data
out <- "rank.files"
if (!dir.exists(out)){
  dir.create(out)
}

ranked_gene <- mapply(function(.x, .y){
  .x <- .x[!is.na(.x$Symbol) & !duplicated(.x$Symbol), ]
  gene_list <- .x[order(.x$log2FoldChange),
                  c("Symbol", "log2FoldChange")]
  write.table(gene_list, file = file.path(out, paste0(.y, ".rnk")),
              sep = "\t", quote = FALSE, 
              row.names =FALSE, col.names =FALSE)
}, DESeq_out, rownames(contrast.matrix), SIMPLIFY = FALSE)




library("biomaRt")
library("readxl")

setwd("C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus\\old.RNA-seq.analysis")

ezh2_gs <- readxl::read_xlsx("EZH2 genesets.xlsx", sheet =1)
convertHumanGeneList <- function(x){
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol",
                   values = x , mart = human, 
                   attributesL = c("mgi_symbol"), 
                   martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

ezh2_gs <- ezh2_gs[-1, ]
Ezh2_mouse_gs <- lapply(ezh2_gs, function(.x){
  gs <- .x[!is.na(.x)]
  gs <- convertHumanGeneList(gs)
})
load(file = "old.RNA-seq.DESeq2.RData")

Ezh2_mouse_gs_log2FC <- lapply(Ezh2_mouse_gs, function(.x){
  fc <- do.call("cbind", lapply(DESeq_out, function(.y){
    tmp_fc <- .y[.y$Symbol %in% .x,  "log2FoldChange"]
    names(tmp_fc) <- .y$Symbol[.y$Symbol %in% .x]
    tmp_fc
  }))
  colnames(fc) <- rownames(contrast.matrix)
  fc
})


plot_sasp_heatmap <- function(exp, 
                              filename_heatmap, height, width =4)
{
  annotation_col<- data.frame(Group = colnames(exp))
  rownames(annotation_col) <- colnames(exp)
  annotation_colors <- list(Group=c("LIL_C-LIL_V" = "#3cb44b",
                                    "LIP_C-LIP_V" = "#bfef45",
                                    "PIL_C-PIL_V" = "#f58231",
                                    "PIP_C-PIP_V" = "#ffe119"))
  
  pdf(filename_heatmap, height = height , width = width)
  heat <- pheatmap(as.matrix(exp), 
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   clustering_method ="ward.D2",
                   clustering_distance_rows = "correlation",
                   scale="none",
                   border_color= NA,
                   treeheight_row = 20,
                   treeheight_col = 10,
                   annotation_col = annotation_col,
                   annotation_colors =annotation_colors,
                   show_colnames= FALSE,
                   show_rownames= TRUE,
                   cellheight = 5,
                   fontsize= 6,
                   fontsize_row = 4)
  print(heat)
  dev.off()
}

null <- mapply(plot_sasp_heatmap, exp = Ezh2_mouse_gs_log2FC, 
               filename_heatmap = paste0(names(Ezh2_mouse_gs_log2FC), ".heatmap.pdf"),
               height = sapply(Ezh2_mouse_gs_log2FC, nrow)*5/72 +0.5) 




con <- file("Ezh2.target.genesets.gmt", open = "w")
null <- mapply(function(.x, .y){
  writeLines(c(.y, "na", .x[-length(.x)]), con = con, sep= "\t")
  writeLines(.x[length(.x)], con = con, sep= "\n")},
  Ezh2_mouse_gs, 
  names(Ezh2_mouse_gs),
  SIMPLIFY = TRUE)
close(con)

## new senescence genesets
senescence_gs <- readxl::read_xlsx("New senescence datasets August 2021.xlsx", sheet =1)
senescence_mouse_gs <- lapply(senescence_gs, function(.x){
  gs <- .x[!is.na(.x)]
  gs <- convertHumanGeneList(gs)
})

con <- file("New.SASP.genesets.gmt", open = "w")
null <- mapply(function(.x, .y){
  writeLines(c(.y, "na", .x[-length(.x)]), con = con, sep= "\t")
  writeLines(.x[length(.x)], con = con, sep= "\n")},
  Ezh2_mouse_gs, 
  names(Ezh2_mouse_gs),
  SIMPLIFY = TRUE)
close(con)

senescence_mouse_gs_log2FC <- lapply(senescence_mouse_gs, function(.x){
  fc <- do.call("cbind", lapply(DESeq_out, function(.y){
    tmp_fc <- .y[.y$Symbol %in% .x,  "log2FoldChange"]
    names(tmp_fc) <- .y$Symbol[.y$Symbol %in% .x]
    tmp_fc
  }))
  colnames(fc) <- rownames(contrast.matrix)
  fc
})

null <- mapply(plot_sasp_heatmap, exp = senescence_mouse_gs_log2FC, 
               filename_heatmap = paste0(names(senescence_mouse_gs_log2FC), ".heatmap.pdf"),
               height = sapply(senescence_mouse_gs_log2FC, nrow)*5/72 +0.5) 


## 12/14/2022 replot heatmap showing expression genesets of interest across replicates
library("ggplot2")
library("sva")
library("DESeq2")
library("edgeR")
library("limma")
library("pheatmap")
library("RColorBrewer")
library("dplyr")
library("ggfortify")
library("genefilter")
library("gplots")
library("doParallel")
library("WriteXLS")


## only significant DEGs
path <- "C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus\\old.RNA-seq.analysis/Table 3.2. DESeq.differential expressed genes_2.xlsx"
DESeq_out_sig <- lapply(1:10, function(i)
  {
     de <- read_xlsx(path,
               sheet = i)
     de <- de[abs(de$log2FoldChange) >= log2(2) & !is.na(de$padj) &
                de$padj <= 0.05, ]
     de
})

names(DESeq_out_sig) <- excel_sheets(path = path)

WriteXLS(x= DESeq_out_sig,
         ExcelFileName = "Table 3.3.old.RNAseq.sig.DESeq.differential.expressed.genes.xlsx",
         row.names = FALSE,
         SheetNames = names(DESeq_out_sig))




plot_heatmap <- function(exp,
                         filename_heatmap, 
                         height, width =4)
{
  annotation_col<- data.frame(Group = gsub("\\d+$", "", colnames(exp),
                                           perl = TRUE))
  rownames(annotation_col) <- colnames(exp)
  annotation_colors <- list(Group=c("LIL_V" = "#999999",
                                    "LIL_C" = "#E69F00",
                                    "PIL_V" = "#56B4E9",
                                    "PIL_C" = "#009E73",
                                    "LIP_V" = "#F0E442",
                                    "LIP_C" = "#0072B2",
                                    "PIP_V" = "#D55E00",
                                    "PIP_C" = "#CC79A7"))

  pdf(filename_heatmap, height = height , width = width)
  heat <- pheatmap(log2(as.matrix(exp+1)), 
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   clustering_method ="ward.D2",
                   clustering_distance_rows = "correlation",
                   scale="row",
                   border_color= NA,
                   treeheight_row = 20,
                   treeheight_col = 10,
                   annotation_col = annotation_col,
                   annotation_colors =annotation_colors,
                   show_colnames= FALSE,
                   show_rownames= TRUE,
                   cellheight = 5,
                   fontsize= 6,
                   fontsize_row = 4)
  print(heat)
  dev.off()
}

null <- mapply(plot_heatmap, exp = geneset_exp, 
               filename_heatmap = paste0(names(geneset_exp), ".heatmap.pdf"),
               height = sapply(geneset_exp, nrow)*5/72 +0.5)


## extract gene expression for genesets
library("readxl")
library("WriteXLS")

setwd("C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus\\old.RNA-seq.analysis")

load("old.RNA-seq.DESeq2.RData")


deg_file <- "Table 3.2. DESeq.differential expressed genes_2.xlsx"
sheets <- excel_sheets(deg_file)

degs <- lapply(1:4, function(x){
     deg <- read_xlsx(deg_file, sheet = x)
})

genesets <- read.delim("12142022.genesets.of.interest4heatmap.txt", 
                       as.is = TRUE, header = TRUE)

## remove "" entries
genesets <- lapply(genesets, function(.x)
{
    .x <- .x[.x != ""]
})

null <- mapply(function(.x, .x.name)
  {
      gs_exp <- lapply(degs, function(.y){
           .y <- .y[.y$Symbol %in% .x, ]
      })
      
      WriteXLS(x = gs_exp, 
               ExcelFileName = paste("Table.6 old.RNAseq", .x.name, "4heatmap.xlsx"),
               row.names = FALSE,
               SheetNames = sheets[1:4])
}, genesets, names(genesets), SIMPLIFY = FALSE)

