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

cur_dir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\UMass_work\Marcus\11-08-2021 RNA-seq\SVA-adjusted DEG analysis)"
setwd(cur_dir)

filename <- "SMA-TK_RNAseq.formatted.txt"

count_PE <- read.delim(filename, sep = "\t",
                       header = TRUE, as.is =TRUE, check.names = FALSE)

rownames(count_PE) <- count_PE[, 1]
count_PE <- count_PE[, -1]

sample_name <- colnames(count_PE)
metadata <- data.frame(sample_id = sample_name,  
                       group = gsub("(.+)\\d+", "\\1", 
                                    sample_name, perl = TRUE),
                       replicates = gsub(".+(\\d+)", "\\1", 
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
keep = rowSums(cpms > 1) >= 4
count_PE = count_PE[keep, ]
rm("cpms", "keep")

dim(count_PE)  ## 13927
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
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
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


## using SVA to adjust hidden variations in the PC3, consider pig responding effects
normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
write.table(normCounts, "Table 1. Normalized.count.table.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

# using sva to reduce hidden variations
get.sva <- function (expr.data=NULL, meta.data=NULL)
{
  mod <- model.matrix(~ group, data=meta.data)
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
  adj.exp <- v$E -fit$coefficients[, 5:7] %*% t(design[, 5:7])
  adj.exp
}

## get sva: 4 hidden variables
meta.sva <- get.sva(expr.data= normCounts, meta.data=metadata)
model <- model.matrix(~ 0 + group + sv1 + sv2 + sv3,
                      data=meta.sva)


## adjust for sva
adj.exp <- get.adj(expr.data=normCounts, design=model,  meta.data=meta.sva)
write.table(adj.exp, "Table 2. sva-adjusted.gene.expression.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

#### Heatmap showing sample distances after adjusting for hidden variation
sampleDists <- dist(t(adj.exp))
pdf("Fig 2.2 Heatmap showing SVA-adjusted sample distances.pdf",
    width = 6, height = 5)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = paste0(vsd$group, vsd$replicates))
dev.off()

## PCA using adjusted value
pc2 <- prcomp(t(adj.exp), scale = T,
              center = T, retx = TRUE)
pdf(file = "Fig 3.2. PCA.plot.SVA.adjusted.samples.pdf",
    width = 12.0, height = 10)
autoplot(pc2, data = metadata,
         colour = 'group', size = 4)
dev.off()

## SVA is much better
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_PE),
                              colData = meta.sva,
                              design = ~0 + group + sv1 + sv2 + sv3)

dds <- estimateSizeFactors(dds)

dds <- estimateDispersions(dds, fitType="parametric",
                           maxit=1000)
normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
write.table(normCounts, "Table 1.1 Normalized.count.table.SVA-adjusted.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

load("SMA-TK.RNA-seq.DESeq2.RData")

model <- model.matrix(~0 + group + sv1 + sv2 + sv3, data = meta.sva)
contrast.matrix <- matrix(c(1, 0, 0, -1, rep(0,3),      ## GCV - V
                            0,  1, 0, -1, rep(0,3),     ## TP - V
                            rep(0, 2), 1, -1, rep(0,3), ## TPGCV - V
                            1, 0, -1, 0, rep(0,3),      ## GCV - TPGCV
                            0, 1, -1, 0, rep(0,3),      ## TP - TPGCV
                            1, -1, 0, 0, rep(0,3)       ## GCV - TP
), nrow = 6, byrow = TRUE)
rownames(contrast.matrix) <- c("GCV-V", 
                               "TP-V",
                               "TPGCV-V",
                               "GCV-TPGCV",
                               "TP-TPGCV",
                               "GCV-TP")

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
  # i=1
  # threshold =1
  # shrink = FALSE
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

save.image(file = "SMA-TK.RNA-seq.DESeq2.RData")

##only output the significant DEGs
DESeq_out_sig <- lapply(DESeq_out, function(.x){
     .x <- .x[abs(.x$log2FoldChange) >= 1 & !is.na(.x$padj) &
                .x$padj <= 0.05, ]
})

WriteXLS(x = DESeq_out_sig, 
         ExcelFileName = "Table 3.2.sig.DESeq.differential.expressed.genes.xlsx",
         row.names = FALSE, 
         SheetNames = rownames(contrast.matrix))


## No SVA adjustment
# dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_PE),
#                               colData = meta.sva,
#                               design = ~0 + group)
# 
# dds <- estimateSizeFactors(dds)
# 
# dds <- estimateDispersions(dds, fitType="parametric",
#                            maxit=1000)
# normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
# write.table(normCounts, "Table 1.1 Normalized.count.table.no.SVA-adjusted.txt", sep = "\t",
#             quote = FALSE, row.names = TRUE)
# 
# model <- model.matrix(~0 + group, data = meta.sva)
# contrast.matrix <- matrix(c(1, 0, 0, -1,      ## GCV - V
#                             0,  1, 0, -1,     ## TP - V
#                             rep(0, 2), 1, -1,  ## TPGCV - V
#                             1, 0, -1, 0,       ## GCV - TPGCV
#                             0, 1, -1, 0      ## TP - TPGCV
# ), nrow = 5, byrow = TRUE)
# rownames(contrast.matrix) <- c("GCV - V", 
#                                "TP - V",
#                                "TPGCV - V",
#                                "GCV - TPGCV",
#                                "TP - TPGCV")
# 
# dds <- nbinomWaldTest(dds, modelMatrix = NULL,
#                       betaPrior=FALSE,
#                       maxit = 50000, 
#                       useOptim = TRUE, 
#                       quiet = FALSE, 
#                       useT = FALSE, useQR = TRUE)
# 
# gene_name <- read.delim("Ensembl.Mus.GCRm38.101.gene.name.txt", 
#                         header = FALSE, as.is = TRUE)
# colnames(gene_name) <- c("Gene_id", "Symbol", "Biotype")
# 
# output_DESeq <- function(i, dds, contrast.matrix, 
#                          threshold, shrink)
# {
#     # i=1
#     # threshold =1
#     # shrink = FALSE
#     res <- results(dds, alpha = 0.01,
#                    contrast=contrast.matrix[i,], 
#                    lfcThreshold=log2(threshold),
#                    format = "DataFrame",
#                    altHypothesis="greaterAbs") 
#     if (shrink)
#     {
#         res <- lfcShrink(dds,
#                          contrast = contrast.matrix[i,],
#                          res = res,
#                          format = "DataFrame",
#                          lfcThreshold = log2(threshold),
#                          type = "ashr")
#     }
#     res <- data.frame(res)
#     counts=normCounts
#     
#     res <- merge(res, counts, by.x = "row.names", 
#                  by.y = "row.names", all.x = TRUE)
#     colnames(res)[1] <- "Gene_id"
#     res <- merge(res, gene_name, by= "Gene_id", all.x = TRUE)
#     res 
# }
# 
# ## apply raw FLC
# DESeq_out <- lapply(1:nrow(contrast.matrix), 
#                     output_DESeq, dds = dds, 
#                     contrast.matrix = contrast.matrix, 
#                     threshold = 1, shrink = FALSE)
# 
# WriteXLS(x = DESeq_out, 
#          ExcelFileName = "Table 3.2. DESeq.differential expressed genes.xlsx",
#          row.names = FALSE, 
#          SheetNames = rownames(contrast.matrix))
# 
# save.image(file = "SMA-TK.RNA-seq.DESeq2.RData")


## prepare GSEA data
out <- "rank.files"
if (!dir.exists(out)){
  dir.create(out)
}

ranked_gene <- mapply(function(.x, .y){
  .x <- .x[!is.na(.x$Symbol) & !duplicated(.x$Symbol), ]
  gene_list <- .x[order(.x$log2FoldChange,  decreasing = TRUE),
                  c("Symbol", "log2FoldChange")]
  write.table(gene_list, file = file.path(out, paste0(.y, ".rnk")),
              sep = "\t", quote = FALSE, 
              row.names =FALSE, col.names =FALSE)
}, DESeq_out, rownames(contrast.matrix), SIMPLIFY = FALSE)


## 
de_genes_exp <- lapply(seq_along(DESeq_out), function(.x){
  expr <- DESeq_out[[.x]]
  expr <- expr[abs(expr$log2FoldChange) > log2(2) &
                 !is.na(expr$padj) & expr$padj < 0.05, ]
  adj.exp <- adj.exp[rownames(adj.exp) %in% expr$Gene_id, ]
  as.data.frame(adj.exp)
})

WriteXLS(x = de_genes_exp, 
         ExcelFileName = "Table 4.1. DESeq.differential.expressed.genes.expression.value4heatmap.xlsx",
         row.names = TRUE, 
         SheetNames = rownames(contrast.matrix))

plot_heatmap <- function(de_genes_exp_cmb, contrast.name, metadata, row_names=NA,
                         filename_heatmap, height, width =4 )
{
  
  if (contrast.name == "GCV - V")
  {
    annotation_col<- data.frame(group = metadata$group[c(1:5, 16:19)])
    rownames(annotation_col) <- colnames(de_genes_exp_cmb)
    annotation_colors <- list(group=c("GCV" = "#3cb44b",
                                      "V" = "#bfef45"))
  } else if (contrast.name == "TP - V") {
    annotation_col<- data.frame(group = metadata$group[c(6:10, 16:19)])
    rownames(annotation_col) <- colnames(de_genes_exp_cmb)
    annotation_colors <- list(group=c("TP" = "#f58231",
                                      "V" = "#ffe119"))
  } else if (contrast.name == "TPGCV - V") {
    annotation_col<- data.frame(group = metadata$group[c(11:15, 16:19)])
    rownames(annotation_col) <- colnames(de_genes_exp_cmb)
    annotation_colors <- list(group=c("TPGCV" = "#3cb44b",
                                      "V" = "#bfef45"))
  } else if (contrast.name == "GCV - TPGCV") {
    annotation_col<- data.frame(group = metadata$group[c(1:5, 11:15)])
    rownames(annotation_col) <- colnames(de_genes_exp_cmb)
    annotation_colors <- list(group=c("GCV" = "#3cb44b",
                                      "TPGCV" = "#bfef45"))
  } else if (contrast.name == "TP - TPGCV") {
    annotation_col<- data.frame(group = metadata$group[c(6:10, 11:15)])
    rownames(annotation_col) <- colnames(de_genes_exp_cmb)
    annotation_colors <- list(group=c("TP" = "#3cb44b",
                                      "TPGCV" = "#bfef45"))
  } else {
    annotation_col<- data.frame(group = metadata$group)
    rownames(annotation_col) <- colnames(de_genes_exp_cmb)
    annotation_colors <- list(group=c("GCV" = "#f58231",
                                      "TP" = "#3cb44b",
                                      "TPGCV" = "#bfef45",
                                      "V" = "#ffe119"))
  }
  show_rownames <- {if (all(!is.na(row_names))) {TRUE} else {FALSE}}
  library("pheatmap")
  pdf(filename_heatmap, height = height , width = width)
  heat <- pheatmap(as.matrix(de_genes_exp_cmb), 
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   clustering_method ="ward.D2",
                   clustering_distance_rows = "correlation",
                   scale="row",
                   border_color= NA,
                   treeheight_row = 20,
                   treeheight_col = 10,
                   annotation_col = annotation_col,
                   annotation_colors =annotation_colors,
                   show_colnames= FALSE,
                   show_rownames= show_rownames,
                   fontsize_row = 5,
                   labels_row = row_names )
  print(heat)
  dev.off()
}

plot_heatmap(de_genes_exp[[1]][, c(1:5, 16:19)], 
             contrast.name = rownames(contrast.matrix)[1],
             metadata = metadata,
             filename_heatmap = "Fig 4.1.GCV-V heatmap.pdf",
             height = 10, width =4)

plot_heatmap(de_genes_exp[[2]][, c(6:10, 16:19)], 
             contrast.name = rownames(contrast.matrix)[2],
             metadata = metadata,
             filename_heatmap = "Fig 4.2.TP-V heatmap.pdf",
             height = 4, width =4)

plot_heatmap(de_genes_exp[[3]][, c(11:15, 16:19)], 
             contrast.name = rownames(contrast.matrix)[3],
             metadata = metadata,
             filename_heatmap = "Fig 4.3.TPGCV -V heatmap.pdf",
             height = 12, width = 6)

plot_heatmap(de_genes_exp[[4]][, c(1:5, 11:15)], 
             contrast.name = rownames(contrast.matrix)[4],
             metadata = metadata,
             filename_heatmap = "Fig 4.4.GCV-TPGCV heatmap.pdf",
             height = 10, width = 6)
plot_heatmap(de_genes_exp[[5]][, c(6:10, 11:15)], 
             contrast.name = rownames(contrast.matrix)[5],
             metadata = metadata,
             filename_heatmap = "Fig 4.5.TP-TPGCV heatmap.pdf",
             height = 4, width = 6)


## two senescence gene sets
library("readxl")

gene_set_file <- "all.mouse.senescence.genelist.xlsx"
gene_set <- read_xlsx(gene_set_file)

## "Custom SASP_Chibaya et al_mouse" only 
custom.gs <- gene_set[, 28] 


adj.exp_annotated <- merge(gene_name, adj.exp, 
                           by.x = "Gene_id", by.y = "row.names", 
                           all.y = TRUE)
adj.exp_annotated <- adj.exp_annotated[adj.exp_annotated$Symbol %in% as.data.frame(custom.gs)[,1], ]
rownames(adj.exp_annotated) <- adj.exp_annotated$Symbol
asap.adj <- adj.exp_annotated[, -c(1:3)]

WriteXLS(x = asap.adj, 
         ExcelFileName = "Table 5. DESeq.differential.expressed.custom.senescence.genes.expression.value4heatmap.xlsx",
         row.names = TRUE)


plot_heatmap(asap.adj, 
             contrast.name = "",
             row_names = rownames(asap.adj),
             metadata = metadata,
             filename_heatmap = "Fig 5.6.custom.ASAP.heatmap.pdf",
             height = 6, width = 4)


## volcano plot
library("EnhancedVolcano")
pdf("Figure 6.V.TP.GCV.TPGCV.Volcano plots.pdf", 
    width = 15.5, height = 10) 
for (i in 1:nrow(contrast.matrix))
{
  p <- EnhancedVolcano(toptable = DESeq_out[[i]],
                       lab = DESeq_out[[i]]$Symbol,
                       x = 'log2FoldChange',
                       y = 'padj',
                       legendLabSize = 10,
                       legendPosition = "right",
                       title = rownames(contrast.matrix)[i],
                       ylab=  bquote(~-Log[10]~italic(adj.P)),
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       legendLabels = c("NS", 
                                        expression(Log[2]~FC), 
                                        "adj.p-value", 
                                        expression(adj.p-value~and~log[2]~FC)),
                       legendIconSize = 8,
                       drawConnectors = TRUE,
                       widthConnectors = 0.5,
                       pointSize = 1.0,
                       xlim = ifelse(i ==1, c(-5, 18), c(-6, 18)),
                       labSize = 4)
  print(p)
}
dev.off()


## Reactome pathway analysis
library("readxl")
library("ReactomePA")
library("clusterProfiler")
library("magrittr")
library("dplyr")
library("biomaRt")

ensembl = useEnsembl(biomart="ensembl",
                     dataset="mmusculus_gene_ensembl", 
                     version=101)

xl <- "Table 3.2. DESeq.differential expressed genes.xlsx"

all_genes <- lapply(1:5, function(.x)
{
  temp <- read_xlsx(xl, sheet =.x)
  mouse_entrez_id <- getBM(attributes=c('ensembl_gene_id',
                                        'entrezgene_id'),
                           filters = 'ensembl_gene_id',
                           values =  temp$Gene_id, mart = ensembl)
  temp <- merge(temp, mouse_entrez_id, 
                by.x = "Gene_id", 
                by.y = "ensembl_gene_id", all.x = TRUE)
})

library("WriteXLS")

samples <- rownames(contrast.matrix)

WriteXLS(x = all_genes, 
         ExcelFileName = "Table 3.3 DESeq.differential expressed genes with Entrez_geneid.xlsx",
         row.names = FALSE, 
         SheetNames = samples)

de_genes <- lapply(all_genes, function(.x){
  up <- .x[.x$log2FoldChange >1 & !is.na(.x$padj) &
             .x$padj <= 0.05 & !is.na(.x$entrezgene_id), ]
  down <- .x[.x$log2FoldChange < -1 & !is.na(.x$padj) &
               .x$padj <= 0.05 & !is.na(.x$entrezgene_id), ]
  list(up=up, down = down)
})

search_kegg_organism('mouse', by='common_name')

directions <- c("up", "down")
for (direction in 1:2){
  de_genes_direction <- lapply(de_genes, "[[", direction)
  
  for (i in seq_along(de_genes_direction))
  {   
    kk <- enrichKEGG(
      gene = as.character(unique(de_genes_direction[[i]]$entrezgene_id)),
      organism = "mmu",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      universe = as.character(all_genes[[i]]$entrezgene_id[!is.na(all_genes[[1]]$entrezgene_id)]),
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      use_internal_data = FALSE)
    write.table(kk, file = paste(samples[i], 
                                 paste0(directions[direction],
                                        "-regulated KEGG pathway overrepresentation analysis.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    tryCatch({pdf(paste(samples[i], 
                        paste0(directions[direction],
                               "-regulated KEGG pathway overrepresentation analysis.pdf")),
                  width = 6, height = 4)
      print(dotplot(kk))},
      error=function(cond) {return(NA)},
      finally= {dev.off()})
    
    reactome <- enrichPathway(gene = as.character(unique(de_genes_direction[[i]]$entrezgene_id)),
                              organism = "mouse",
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.2,
                              universe = as.character(all_genes[[i]]$entrezgene_id[!is.na(all_genes[[i]]$entrezgene_id)]),
                              minGSSize = 10,
                              maxGSSize = 500,
                              readable = TRUE)
    write.table(reactome, 
                file = paste(samples[i], 
                             paste0(directions[direction], 
                                    "-regulated Reactome pathway overrepresentation analysis.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    tryCatch({pdf(paste(samples[i], 
                        paste0(directions[direction],
                               "-regulated Reactome pathway overrepresentation analysis.pdf")),
                  width = 6, height = 4)
      print(dotplot(reactome))},
      error=function(cond) {return(NA)},
      finally = {dev.off()})
  }  
}
save.image(file = "03112022.Marcus.RNA-seq.RData")



## update figures
cur_dir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\UMass_work\Marcus\11-08-2021 RNA-seq\SVA-adjusted DEG analysis)"
setwd(cur_dir)
load("03112022.Marcus.RNA-seq.RData")

##
install.packages("enrichR")

library(enrichR)
library(cowplot)
library(patchwork)
listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()

#[1] "TRRUST_Transcription_Factors_2019"        
#[2] "ChEA_2016"                                
#[3] "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"
#[4] "ENCODE_TF_ChIP-seq_2015"                  
#[5] "Genome_Browser_PWMs"                      
#[6] "TRANSFAC_and_JASPAR_PWMs"
#[7] "TF-LOF_Expression_from_GEO"

dbs <- dbs$libraryName[c(147, 104, 96, 50, 1, 2, 12)]

up_dow_genes <- mapply(function(.x, .y){
  up <- .x[.x$log2FoldChange >1 & !is.na(.x$padj) &
             .x$padj <= 0.05, "Symbol"]
  down <- .x[.x$log2FoldChange < -1 & !is.na(.x$padj) &
               .x$padj <= 0.05, "Symbol"]
  all <- .x[abs(.x$log2FoldChange) > 1 & !is.na(.x$padj) &
              .x$padj <= 0.05, "Symbol"]
  gene_list <- list(up = up, down = down, all= all)
  null <- lapply(seq_along(gene_list), function(k) {
    cat(names(gene_list)[k], ":\t", sep ="",  file = paste(.y, "DEGs.txt"), append = TRUE)
    cat(gene_list[[k]], sep = "\t", file = paste(.y, "DEGs.txt"), append = TRUE)
    cat("\n", file = paste(.y, "DEGs.txt"), append = TRUE)
  })
  
  null <- mapply(function(x, y) {
    
    enriched <- enrichr(x, dbs)
    
    out <- paste0(.y, y, "-genes.TF.target.enrichment.out")
    if (!dir.exists(out)) {
      dir.create(out)
    }
    null <-  mapply(function(i, j){
      write.table(i, file = file.path(out, paste0(j, "_table.txt")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }, enriched, names(enriched), SIMPLIFY = FALSE)
    
  }, gene_list, names(gene_list), SIMPLIFY = FALSE)
}, all_genes, samples, SIMPLIFY = FALSE)


## visualize the EnrichR result
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\UMass_work\Marcus\11-08-2021 RNA-seq\SVA-adjusted DEG analysis\EnrichR)")
folders <- dir(".", "target.enrichment.out$", full.names = TRUE)

out <- "visualization_top10"
if (!dir.exists(out)){
  dir.create(out)
}

for (f in folders){
  
  files <- dir(f, full.names = TRUE)
  out_names <- gsub(".+/(.+?)-genes.TF.target.enrichment.out", "\\1", f, perl = TRUE)
  set_names <- gsub(".+/(.+?)_table.txt", "\\1", files, perl = TRUE)
  tf <- mapply(function(.x, .y){
    tf_enrichment <- read.delim(.x, header = TRUE, as.is = TRUE)[, c(1, 4, 8)]
    tf_enrichment <- tf_enrichment[order(tf_enrichment$Combined.Score, 
                                         decreasing = TRUE), ]
    tf_enrichment <- tf_enrichment[1:10, ]
    if (nrow(tf_enrichment) > 0) {
      tf_enrichment$set <- .y
    }
    
    tf_enrichment <- tf_enrichment[order(tf_enrichment$Combined.Score), ]
    tf_enrichment$Term <- factor(tf_enrichment$Term, 
                                 levels = tf_enrichment$Term)
    tf_enrichment
  }, files, set_names, SIMPLIFY = FALSE, USE.NAMES = TRUE)
  
  pdf(file.path(out, paste0(out_names, "enriched.transcription.factors.pdf")), 
      height = 8, width = 11)
  p <- lapply(tf, function(x) {
    ggplot(x, aes(Term, Combined.Score, fill = set)) + 
      geom_bar(stat = "identity",  width= 0.5, 
               position = position_dodge(width=0.1)) + 
      coord_flip() + 
      ylab("Combined Score") +
      xlab("Term") + theme_cowplot(font_size = 8) +
      theme(legend.position="top", 
            axis.text = element_text(size = 6)) 
  })
  
  print(wrap_plots(p, nrow = 4, ncol = 2))
  dev.off()
}
