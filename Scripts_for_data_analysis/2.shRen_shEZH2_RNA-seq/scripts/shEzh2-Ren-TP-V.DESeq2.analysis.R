install.packages("BiocManager")
pkgs <- c("ggplot2", "sva","DESeq2", "edgeR", "limma",
          "pheatmap", "RColorBrewer", "dplyr", "ggfortify",
          "genefilter", "gplots", "doParallel", "readxl",
          "WriteXLS", "EnhancedVolcano", "msigdbr", 
          "ReactomePA", "clusterProfiler", "magrittr",
          "biomaRt")
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
library("readxl")
library("WriteXLS")
library("EnhancedVolcano")
library("msigdbr")
library("ReactomePA")
library("clusterProfiler")
library("magrittr")
library("biomaRt")

setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\UMass_work\Marcus\RNA-seq 202103)")

filename <- "Marcus_RNAseq.PE.formatted.txt"
count_PE <- read.delim(filename, header = TRUE, 
                       as.is =TRUE, check.names = FALSE)

sample_name <- colnames(count_PE)[-1]
metadata <- data.frame(sample_id = sample_name,  
                       group = gsub("(.+?[CV]).+", "\\1", 
                                    sample_name, perl = TRUE),
                       replicates = gsub(".+?[CV](.+)", "\\1", 
                                         sample_name, perl = TRUE))
metadata$group <- factor(metadata$group, 
                         levels = unique(metadata$group),
                         labels = c("shEzh2_TP","shEzh2_Veh",
                                    "shRen_TP","shRen_Veh"))
rownames(count_PE) <- count_PE[, 1]
count_PE <- count_PE[, -1]
count_PE <- count_PE[rowSums(count_PE) != 0, ]

## total reads assigned to gene features
total <- data.frame(total = colSums(count_PE))
total <- cbind(metadata, total)
total$sample_id <- paste0(total$group, total$replicates)
pdf("Fig.1. Total number of reads assigned to gene features.pdf",
    width = 6, height =5)
ggplot(total, aes(sample_id, total, fill = group))+
  geom_col() + 
  theme(axis.text.x = 
          element_text(angle=90, hjust = 0, vjust = 0.5))
dev.off()

## remove lowly expressed genes
cpms = cpm(count_PE)
keep = rowSums(cpms > 1) >= 5
count_PE = count_PE[keep, ]
rm("cpms", "keep")

## sanity check
all(colnames(count_PE) == metadata$sample_id)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_PE),
                              colData = metadata,
                              design = ~0 + group)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType="parametric", 
                           maxit=1000)
## exploratory analysis
vsd <- vst(dds, blind = TRUE)
sampleDists <- dist(t(assay(vsd)))

## Heatmap showing sample distances
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

pdf("Fig 2. Heatmap showing sample distances.pdf", 
    width = 4.5, height = 4)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = paste0(vsd$group,
                                                      vsd$replicates))
dev.off()

## PCA plot
pc <- prcomp(t(assay(vsd)), scale = T, center = T, retx = TRUE)
pdf(file = "Fig 3. PCA.plot.pdf", width = 12.0, height = 10)
autoplot(pc, data = metadata, colour = 'group', size = 4)
dev.off()

normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
write.table(normCounts, "Table 1. Normalized.count.table.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

## using SVA to adjust for hidden variations
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
## get sva: 4 hidden variables
meta.sva <- get.sva(expr.data= normCounts, meta.data=metadata)
model <- model.matrix(~ 0 + group + sv1 + sv2 + sv3 + sv4 + sv5 + sv6,
                      data=meta.sva)

## adjust expression for hidden variations for EDA plots
get.adj <- function(expr.data=NULL, design=NULL,  meta.data=NULL)
{
  ## adjusted expression using voom()
  v <- voom(expr.data, design=design)
  fit <- lmFit(v, design=design)
  adj.exp <- v$E -fit$coefficients[, 5:9] %*% t(design[, 5:9])
  adj.exp
}

adj.exp <- get.adj(expr.data=normCounts, design=model, meta.data=meta.sva)
write.table(adj.exp, "Table 2. sva-adjusted.gene.expression.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

## Heatmap showing sample distances after adjusting for hidden variation
sampleDists <- dist(t(adj.exp))
pdf("Fig 2.2 Heatmap showing SVA-adjusted sample distances.pdf", 
    width = 6, height = 5)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists, 
                                 sampleNames = paste0(vsd$group, vsd$replicates))
dev.off()

## PCA using adjusted value
pc2 <- prcomp(t(adj.exp), scale = T, center = T, retx = TRUE)
pdf(file = "Fig 3.2. PCA.plot.SVA.adjusted.samples.pdf", 
    width = 12.0, height = 10)
autoplot(pc2, data = metadata, colour = 'group', size = 4)
dev.off()

## SVA is much better
dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_PE),
                              colData = meta.sva,
                              design = ~0 + group + sv1 +
                                sv2 + sv3 + sv4 + sv5 + sv6)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType="parametric", 
                           maxit=1000)
normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
write.table(normCounts, 
            "Table 1.1 Normalized.count.table.SVA-adjusted.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

load("04052021.Marcus.RNA-seq.RData")
model <- model.matrix(~0 + group + sv1 + sv2 + sv3 + sv4 + sv5 + sv6, data = meta.sva)

contrast.matrix <- matrix(c(1, -1, 0,  0, rep(0,6), ## shEzh2_TP - shEzh2_Veh
                            0, 0,  1, -1, rep(0,6), ## shRen_TP - shRen_Veh
                            1, -1, -1, 1, rep(0,6), ## interaction
                            1, 0,  -1, 0, rep(0,6), ## shEzh2_TP - shRen_TP
                            0, 1,  0, -1, rep(0,6), ## shEzh2_Veh - shRen_Veh
                            1, 0,  0, -1, rep(0,6), ## Ezh2_TP-Ren_V
                            0, 1,  -1, 0, rep(0,6)  ## EZh2_V-Ren_TP
), nrow = 7, byrow = TRUE)
rownames(contrast.matrix) <- c("Ezh2_TP-Ezh2_V", 
                               "Ren_TP-Ren_V",
                               "(Ezh2_TP-Ezh2_V)-(Ren_TP-Ren_V)",
                               "Ezh2_TP-Ren_TP", 
                               "Ezh2_V-Ren_V",
                               "Ezh2_TP-Ren_V",
                               "EZh2_V-Ren_TP")

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
  
  if (i ==1){
    counts=normCounts[, c(1:6, 12:16)]
  } else if (i==2) {
    counts=normCounts[, c(7:11, 17:21)]
  } else {
    counts=normCounts
  }
  res <- merge(res, counts, by.x = "row.names", 
               by.y = "row.names", all.x = TRUE)
  colnames(res)[1] <- "Gene_id"
  res <- merge(res, gene_name, by= "Gene_id", all.x = TRUE)
  res 
}

## raw FLC
DESeq_out <- lapply(1:nrow(contrast.matrix), 
                    output_DESeq, dds = dds, 
                    contrast.matrix = contrast.matrix, 
                    threshold = 1, shrink = FALSE)

WriteXLS(x = DESeq_out, 
         ExcelFileName = "Table 3.2. DESeq.differential expressed genes.xlsx",
         row.names = FALSE, 
         SheetNames = rownames(contrast.matrix))

##only output the significant DEGs
DESeq_out_sig <- lapply(DESeq_out, function(.x){
  .x <- .x[abs(.x$log2FoldChange) >= 1 & !is.na(.x$padj) &
             .x$padj <= 0.05, ]
})

WriteXLS(x = DESeq_out_sig, 
         ExcelFileName = "Table 3.2.sig.DESeq.differential.expressed.genes.xlsx",
         row.names = FALSE, 
         SheetNames = rownames(contrast.matrix))

save.image(file = "04052021.Marcus.RNA-seq.RData")


de_genes_exp <- lapply(seq_along(DESeq_out), function(.x){
  expr <- DESeq_out[[.x]]
  expr <- expr[abs(expr$log2FoldChange) > log2(2) &
                 !is.na(expr$padj) & expr$padj < 0.05, ]
  adj.exp <- adj.exp[rownames(adj.exp) %in% expr$Gene_id, ]
  as.data.frame(adj.exp)
})
WriteXLS(x =de_genes_exp, 
         ExcelFileName = 
           "Table 4.1. DESeq.differential.expressed.genes.expression.value4heatmap.xlsx",
         row.names = TRUE, 
         SheetNames = rownames(contrast.matrix))
plot_heatmap <- function(de_genes_exp_cmb, metadata, row_names=NA,
                         filename_heatmap, height, width =4 )
{
  if (ncol(de_genes_exp_cmb) == 11)
  {
    colnames(de_genes_exp_cmb) <- paste0(metadata$group[c(1:6, 12:16)], 
                                         metadata$replicates[c(1:6, 12:16)])
    
    annotation_col<- data.frame(group = metadata$group[c(1:6, 12:16)])
    rownames(annotation_col) <- colnames(de_genes_exp_cmb)
    annotation_colors <- list(group=c("shEzh2_TP" = "#3cb44b",
                                      "shRen_TP" = "#bfef45"))
  } else if (ncol(de_genes_exp_cmb) == 10) {
    colnames(de_genes_exp_cmb) <- 
      paste0(metadata$group[c(7:11, 17:21)], 
             metadata$replicates[c(7:11, 17:21)])
    
    annotation_col<- data.frame(group = metadata$group[c(7:11, 17:21)])
    rownames(annotation_col) <- colnames(de_genes_exp_cmb)
    annotation_colors <- list(group=c("shEzh2_Veh" = "#f58231",
                                      "shRen_Veh" = "#ffe119"))
  } else {
    colnames(de_genes_exp_cmb) <- paste0(metadata$group[1:21], 
                                         metadata$replicates[1:21])
    
    annotation_col<- data.frame(group = metadata$group[1:21])
    rownames(annotation_col) <- colnames(de_genes_exp_cmb)
    annotation_colors <- list(group=c("shEzh2_TP" = "#3cb44b",
                                      "shEzh2_Veh" = "#bfef45",
                                      "shRen_TP" = "#f58231",
                                      "shRen_Veh" = "#ffe119"))
  }                               
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
                   show_rownames= FALSE,
                   fontsize_row = 4,
                   labels_row = row_names )
  print(heat)
  dev.off()
}

plot_heatmap(de_genes_exp[[1]][, 1:11], metadata = metadata,
             filename_heatmap = "Fig 4.1.shEzh2_TP-shEzh2_Veh heatmap.pdf",
             height = 10, width =4)

plot_heatmap(de_genes_exp[[2]][, 12:21], metadata = metadata,
             filename_heatmap = "Fig 4.2.shRen_TP-shRen_Veh heatmap.pdf",
             height = 6, width =4)

plot_heatmap(de_genes_exp[[3]], metadata = metadata,
             filename_heatmap = "Fig 4.3.shEzh2_TP-shEzh2_Veh -(shRen_TP-shRen_Veh) heatmap.pdf",
             height = 6, width = 6)

plot_heatmap(de_genes_exp[[1]][, c(1:6, 12:16)], metadata = metadata,
             filename_heatmap = "Fig 4.4.shEzh2_TP-shRen_TP heatmap.pdf",
             height = 10, width = 6)
plot_heatmap(de_genes_exp[[2]][, c(7:11, 17:21)], metadata = metadata,
             filename_heatmap = "Fig 4.5.shEzh2_Veh-shRen_Veh heatmap.pdf",
             height = 6, width = 6)


## two senescence gene sets
gene_set_file <- "all.mouse.senescence.genesets.xlsx"
gene_set <- read_xlsx(gene_set_file)

## Fridman and Tasdemir senescence gene sets
gs <- gene_set[, c(17,18,26)] 
Fridman <- as.vector(as.matrix(gs[-1, 1:2]))
Fridman <- Fridman[!is.na(Fridman)]
Tasdemir <- as.vector(as.matrix(gs[-1, 3]))
Tasdemir <- Tasdemir[!is.na(Tasdemir)]

de_genes_exp <- lapply(seq_along(DESeq_out), function(.x){
  expr <- DESeq_out[[.x]]
  expr <- expr[abs(expr$log2FoldChange) > log2(2) &
                 !is.na(expr$padj) & expr$padj < 0.05, ]
  adj.exp <- adj.exp[rownames(adj.exp) %in% expr$Gene_id, ]
  
  adj.exp_annotated <- merge(gene_name, adj.exp, 
                             by.x = "Gene_id", by.y = "row.names", 
                             all.y = TRUE)
  adj.exp
})

de_genes_exp_senescence <- lapply(de_genes_exp, function(x){
  x <- x[x$Symbol %in% c(Fridman, Tasdemir), ]
  rownames(x) <- x$Symbol
  x <- x[, -c(1:3)]
})

WriteXLS(x =de_genes_exp_senescence, 
         ExcelFileName = 
           paste0("Table 5.DESeq.differential.expressed.Fridman.and.Tasdemir.",
                  "senescence.genes.expression.value4heatmap.xlsx"),
         row.names = TRUE, 
         SheetNames = c("shEzh2_TP-shEzh2_Veh", 
                        "shRen_TP-shRen_Veh"))


Fridman_gs_expr <- 
  de_genes_exp_senescence[[1]][rownames(de_genes_exp_senescence[[1]]) %in%
                                 Fridman, 1:11]
plot_heatmap(Fridman_gs_expr, 
             metadata = metadata, row_names = rownames(Fridman_gs_expr),
             filename_heatmap = "Fig 5.1.shEzh2_TP-shEzh2_Veh.Fridman.et.al.heatmap.pdf",
             height = 2, width = 4)

Tasdemir_gs_expr <- 
  de_genes_exp_senescence[[1]][rownames(de_genes_exp_senescence[[1]]) %in%
                                 Tasdemir, 1:11]
plot_heatmap(Tasdemir_gs_expr, metadata = metadata,
             row_names = rownames(Tasdemir_gs_expr),
             filename_heatmap = "Fig 5.1.shEzh2_TP-shEzh2_Veh.Tasdemir.et.al.heatmap.pdf",
             height = 4, width =4 )

WriteXLS(x = list(Fridman_Senescence = Fridman_gs_expr, 
                  Tasdemir_Senescence = Tasdemir_gs_expr), 
         ExcelFileName = paste0("Table 6.1. shEzh2_TP-shEzh2_Veh.DE.Fridman.and.",
                                "Tasdemir.senescence.genes.expression.value4heatmap.xlsx"),
         row.names = TRUE, 
         SheetNames = c("Fridman_Senescence", "Tasdemir_Senescence"))


Fridman_gs_expr <- 
  de_genes_exp_senescence[[2]][rownames(de_genes_exp_senescence[[2]]) %in%
                                 Fridman, 12:21]
plot_heatmap(Fridman_gs_expr, 
             metadata = metadata, row_names = rownames(Fridman_gs_expr),
             filename_heatmap = "Fig 5.2.shRen_TP-shRen_Veh.Fridman.et.al.heatmap.pdf",
             height = 1.5, width = 4)

Tasdemir_gs_expr <- 
  de_genes_exp_senescence[[2]][rownames(de_genes_exp_senescence[[2]]) %in%
                                 Tasdemir, 12:21]
plot_heatmap(Tasdemir_gs_expr, metadata = metadata,
             row_names = rownames(Tasdemir_gs_expr),
             filename_heatmap = "Fig 5.2.shRen_TP-shRen_Veh.Tasdemir.et.al.heatmap.pdf",
             height = 4, width =4 )

WriteXLS(x = list(Fridman_Senescence = Fridman_gs_expr,
                  Tasdemir_Senescence = Tasdemir_gs_expr), 
         ExcelFileName = 
           paste0("Table 6.2. shRen_TP-shRen_Veh.DE.Fridman.and.Tasdemir",
                  ".senescence.genes.expression.value4heatmap.xlsx"),
         row.names = TRUE, 
         SheetNames = c("Fridman_Senescence", "Tasdemir_Senescence"))



## volcano plot
pdf("Figure 5.shEzh2_TP-shRen_TP and shEzh2_V-shRen_V Volcano plots.pdf", 
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

all_gene_sets = msigdbr(species = "Mus musculus")
msigdbr_collections()

gs_catogery <- unique(all_gene_sets$gs_cat)

for (group in seq_along(gs_catogery))
{
  geneset_category <- all_gene_sets %>%
    dplyr::filter(gs_cat == gs_catogery[group]) %>% 
    select(gs_subcat, gs_name, gene_symbol)
  
  geneset_category <- split(as.data.frame(geneset_category),
                            f = geneset_category$gs_subcat)
  if (names(geneset_category)[1] == "") {
    names(geneset_category)[1] <- gs_catogery[group]
  }
  for (subgroup in names(geneset_category))
  {
    fn <- paste0("Mus.geneset ", gs_catogery[group], 
                 gsub(":", "-", subgroup, perl = TRUE), 
                 ".gmt")
    sink(file = fn,
         append = TRUE)
    temp <- geneset_category[[subgroup]]
    temp_gs <- split(temp,
                     f = temp$gs_name)
    for (gs_n in names(temp_gs))
    {
      cat(gs_n, "na", 
          temp_gs[[gs_n]]$gene_symbol, 
          sep = "\t")
      cat("\n")
    }
    
    sink()
  }
}

## Reactome pathway analysis
ensembl = useEnsembl(biomart="ensembl",
                     dataset="mmusculus_gene_ensembl", 
                     version=101)

xl <- "Table 3.2. DESeq.differential expressed genes.xlsx"

all_genes <- lapply(1:2, function(.x)
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

samples <- c("shEzh2-TP_shEzh2-veh", "shRen-TP_shRen-veh",
             "(Ezh2-TP_Ezh2-veh)-(Ren-TP_Ren-veh)", 
             "shEzh2-TP_shRen-TP", "shEzh2-V_shRen-V")

WriteXLS(x = all_genes, 
         ExcelFileName = paste0("Table 3. DESeq.differential expressed genes",
                                " with Entrez_geneid.xlsx"),
         row.names = FALSE, 
         SheetNames = samples)

de_genes <- lapply(all_genes, function(.x){
  up <- .x[.x$log2FoldChange >1 & !is.na(.x$padj) &
             .x$padj <= 0.5 & !is.na(.x$entrezgene_id), ]
  down <- .x[.x$log2FoldChange < -1 & !is.na(.x$padj) &
               .x$padj <= 0.5 & !is.na(.x$entrezgene_id), ]
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
    write.table(kk, file = 
                  paste(samples[i],
                        paste0(directions[direction], 
                               "-regulated KEGG pathway overrepresentation analysis.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    pdf(paste(samples[i], paste0(directions[direction],
                                 "-regulated KEGG pathway overrepresentation analysis.pdf")),
        width = 8, height = 6)
    print(dotplot(kk))
    dev.off()
    
    reactome <- enrichPathway(gene = as.character(unique(de_genes_direction[[i]]$entrezgene_id)),
                              organism = "mouse",
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.2,
                              universe = as.character(all_genes[[i]]$entrezgene_id[!is.na(all_genes[[i]]$entrezgene_id)]),
                              minGSSize = 10,
                              maxGSSize = 500,
                              readable = TRUE)
    write.table(reactome, file = paste(samples[i], paste0(directions[direction], 
                                                          "-regulated Reactome pathway overrepresentation analysis.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    pdf(paste(samples[i], paste0(directions[direction],
                                 "-regulated Reactome pathway overrepresentation analysis.pdf")),
        width = 14, height = 8)
    print(dotplot(reactome))
    dev.off()
  }  
}

## senescence genes
gene_set_file <- "all.mouse.senescence.genesets.xlsx"
gene_set <- read_xlsx(gene_set_file)
pooled_genesets <- unique(as.vector(as.matrix(as.data.frame(gene_set[-c(1), ]))))

de_genes <- lapply(c(1,2), function(.x) {
  temp <- read_xlsx("Table 3. DESeq.differential expressed genes with Entrez_geneid.xlsx",
                    sheet = .x)
  #temp <- temp[abs(temp$log2FoldChange) >= log2(2) & !is.na(temp$padj) & temp$padj < 0.05, c(1, 3, 7,19)]
  temp <- temp[temp$Symbol %in% pooled_genesets, c(1, 3, 7,19)]
})

## volcano plots
contrasts_name <- c("shEzh2_TP-shEzh2_Veh", "shRen_TP-shRen_Veh")
pdf("Fig 6. Volcano plots for ASAP genes.pdf", width = 15, height = 10) 
for (i in 1:2)
{
  p<- EnhancedVolcano(toptable = de_genes[[i]],
                      lab = de_genes[[i]]$Symbol,
                      x = 'log2FoldChange',
                      y = 'padj',
                      legendLabSize = 10,
                      legendPosition = "right",
                      title = contrasts_name[i],
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
                      labSize = 4)
  print(p)
}
dev.off()

## heatmap
adj_expr <- read.delim("Table 2. sva-adjusted.gene.expression.txt", 
                       header = TRUE)

de_genes_exp <- lapply(seq_along(de_genes), function(.x){
  expr <- de_genes[[.x]]
  expr <- expr[abs(expr$log2FoldChange) > log2(2) &
                 !is.na(expr$padj) & expr$padj < 0.05, ]
  adj_expr <- adj_expr[rownames(adj_expr) %in% expr$Gene_id, ]
  adj_expr <- merge(as.data.frame(adj_expr), expr[, c("Gene_id", "Symbol")],
                    by.x = "row.names", by.y = "Gene_id", all.x = TRUE)
  adj_expr <- unique(adj_expr)
  rownames(adj_expr) <- adj_expr[, 1]
  adj_expr <- adj_expr[, -1]
})

plot_heatmap(de_genes_exp[[1]][, 1:11], metadata = metadata, 
             row_names = de_genes_exp[[1]][, 22],
             filename_heatmap = "Fig 4.1 ASAP genes shEzh2_TP-shEzh2_Veh heatmap.pdf",
             height = 20, width =4)

plot_heatmap(de_genes_exp[[2]][, 12:21], 
             metadata = metadata, 
             row_names = de_genes_exp[[1]][, 22],
             filename_heatmap = 
               "Fig 4.2 ASAP genes shRen_TP-shRen_Veh heatmap.pdf",
             height = 12, width =4)

## differentially expressed ASAP genes
WriteXLS::WriteXLS(x= de_genes, 
                   ExcelFileName = "Table 5. DESeq.differential expressed ASAP.genes.xlsx",
                   row.names = TRUE, 
                   SheetNames = contrasts_name)

## Adding TF annotation
mouse_TF <- read.delim("Mus_musculus_TF.cofactors.animalTFDB3.0.txt", 
                       header = TRUE)

de_genes <- lapply(c(1,2), function(.x) {
  temp <- read_xlsx("Table 3. DESeq.differential expressed genes with Entrez_geneid.xlsx",
                    sheet = .x)
  temp <- merge(temp, mouse_TF[, c(3,4,6)], 
                by.x = "Gene_id", by.y = "Ensembl", 
                all.x = TRUE)
})

WriteXLS::WriteXLS(x= de_genes, 
                   ExcelFileName = "Table 6. DESeq.differential expressed gene.TF.annotated.xlsx",
                   row.names = FALSE, 
                   SheetNames = contrasts_name)

## DE ASAP Gene list 
de_genes_exp <- lapply(seq_along(de_genes), function(.x){
  expr <- de_genes[[.x]]
  expr <- expr[abs(expr$log2FoldChange) > log2(2) &
                 !is.na(expr$padj) & expr$padj < 0.05, ]
  
  expr$regulation <- ifelse(expr$log2FoldChange > 0, "up", "down")
  expr
})

WriteXLS::WriteXLS(x= de_genes_exp, 
                   ExcelFileName = "Table 7. DESeq.differential expressed ASAP genelist.xlsx",
                   row.names = FALSE, 
                   SheetNames = contrasts_name)

de_genes <- lapply(c(1,2), function(.x) {
  expr <- read_xlsx("Table 3. DESeq.differential expressed genes with Entrez_geneid.xlsx",
                    sheet = .x)
  expr <- expr[abs(expr$log2FoldChange) > log2(2) &
                 !is.na(expr$padj) & expr$padj < 0.05, ]
  expr$regulation <- ifelse(expr$log2FoldChange > 0, "up", "down")
  expr
})

WriteXLS::WriteXLS(x= de_genes, 
                   ExcelFileName = "Table 8. DESeq.differential expressed genelist.xlsx",
                   row.names = FALSE, 
                   SheetNames = contrasts_name)

## replot heatmap for interaction with row labels
f <- "Table 4. DESeq.differential.expressed.genes.expression.value4heatmap.xlsx"
inter <- as.data.frame(read_xlsx(f, sheet =3))
rownames(inter) <- inter$Symbol 
inter <- inter[, -c(1:3)]

load("04052021.Marcus.RNA-seq.RData")
plot_heatmap(inter, metadata = metadata, row_names = rownames(inter),
             filename_heatmap = "Fig 4.3.shEzh2_TP-shEzh2_Veh -(shRen_TP-shRen_Veh) heatmap with gene names.pdf",
             height = 25, width = 6)

## plot heatmap for SASP genes across samples
load("04052021.Marcus.RNA-seq.RData")
sasp <- read_xlsx("all.mouse.senescence.genelist.xlsx", sheet = 1)
sasp <- sasp[, 27]
extra_genes <- c("CXCL9", "CCL20", "ICAM1", "VEGFb", 
                 "PDGFa", "PDGFb", "IL12a", "IL12b", 
                 "MMP2")
extra_genes <- paste0(gsub("(^.).+", "\\1", extra_genes, perl = TRUE),
                      tolower(gsub("^.(.+)", "\\1", extra_genes, perl = TRUE)))
target_genes <- c(extra_genes, as.vector(as.matrix(sasp)))
target_genes <- unique(target_genes[!is.na(target_genes)])
target_genes <- list(target_genes = target_genes)

adj.exp <- merge(adj.exp, gene_name, by.x = "row.names",
                 by.y = "Gene_id", all.x = TRUE)

sasp_expr <- lapply(target_genes, function(.x){
  dat <- adj.exp[adj.exp$Symbol %in% .x, ]
  dat <- dat[, c(2:23)]
  rownames(dat) <- dat[, 22]
  dat <- dat[, -22]
})

plot_sasp_heatmap <- function(exp, metadata,
                              filename_heatmap, height, width =4)
{
  colnames(exp) <- paste0(metadata$group, 
                          metadata$replicates)
  
  annotation_col<- data.frame(Group = metadata$group)
  rownames(annotation_col) <- colnames(exp)
  annotation_colors <- list(Group=c("shEzh2_TP" = "#3cb44b",
                                    "shEzh2_Veh" = "#bfef45",
                                    "shRen_TP" = "#f58231",
                                    "shRen_Veh" = "#ffe119"))
  
  pdf(filename_heatmap, height = height , width = width)
  heat <- pheatmap(as.matrix(exp), 
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
                   show_rownames= TRUE,
                   cellheight = 5,
                   fontsize= 6,
                   fontsize_row = 4)
  print(heat)
  dev.off()
}

null <- mapply(plot_sasp_heatmap, exp=sasp_expr, 
               filename_heatmap = paste0(names(sasp_expr), ".heatmap.2.pdf"),
               height = sapply(sasp_expr, nrow)*5/72 +0.5, 
               MoreArgs = list(metadata = metadata)) 