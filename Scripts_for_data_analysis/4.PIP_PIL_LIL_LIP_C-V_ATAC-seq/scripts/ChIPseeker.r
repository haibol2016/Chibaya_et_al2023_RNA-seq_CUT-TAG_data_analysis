library("ChIPseeker")
library("GenomicFeatures")
library("ChIPpeakAnno")
library("ggplot2")

## make TxDb from GTF:prepare their own TxDb object by retrieving information from UCSC Genome Bioinformatics and BioMart data resources by R function makeTxDbFromBiomart and  makeTxDbFromUCSC. 
# hsa_TxDb <- makeTxDbFromGFF(file = "~/work/mccb/genome/Homo_sapiens/human_gencode_v34/gencode.v34.primary_assembly.annotation.gtf", format="gtf", dataSource="Gencode", organism="Sus scrofa", taxonomyId= NA, circ_seqs="chrM", chrominfo = NULL, miRBaseBuild=NA)
# # 
# saveDb(hsa_TxDb, file="~/work/mccb/genome/Homo_sapiens/human_gencode_v34/hsa.gencode.v34.sqlite")

mmu_TxDb <- loadDb("/home/hl84w/work/mccb/genome/Mus_musculus_GRCm38.p5/Mus_musculus.GRCm38.p6.v101.sqlite")
# peak_files <- dir("./", "Peak$")
peak_files <- dir(".", "q005.consensus.peaks.col6.bed$", full.names = TRUE)

macsOutput <- lapply(peak_files, function(.x){
    toGRanges(.x, format = "BED")
}) 

names(macsOutput) <- c("consensusPeaks.q005")

## genome-wide coverage plot
pdf("Figure1.Genome-wide ATAC-seq peak distribution.pdf", width =12, height = 10)
null <- lapply(seq_along(macsOutput), function(.x){
    print(covplot(macsOutput[[.x]], weightCol = "score",
                  title = names(macsOutput)[[.x]],
                  chrs = paste0(c(as.character(1:19), "X", "Y"))))
})
dev.off()

## Profile of ChIP peaks binding to TSS regions
pdf("Figure2.consensus.peak.tagHeatmap.2.TSS.pdf", width = 2, height =6)

    promoter <- getPromoters(TxDb=mmu_TxDb, upstream=3000, downstream=3000)
    tagMatrixList <- lapply(macsOutput, getTagMatrix, windows=promoter)
    names(tagMatrixList) <- names(macsOutput)
    tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color="blue")
    
dev.off()

pdf("Figure3.consensus.peak.AvgProfile.2.TSS.pdf", width = 6, height = 3)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), facet="none") + theme(legend.position="bottom")
dev.off()

# ## Average Profile of ChIP peaks binding to TSS region
# 
# #### one-step profile plot
# plotAvgProf2(files[[4]], TxDb=txdb, upstream=3000, downstream=3000,
#              xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
# 
# ## adding confidence interval
# plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)


## peak annotation: The output of annotatePeak is csAnno instance. ChIPseeker provides as.GRanges to convert csAnno to GRanges instance, and as.data.frame to convert csAnno to  data.frame which can be exported to file by write.table.ChIPseeker adopted the following priority in genomic annotation.

# Promoter
# 5' UTR
# 3' UTR
# Exon
# Intron
# Downstream
# Intergenic

# ChIPseeker also provides parameter genomicAnnotationPriority for user to prioritize this hierachy.

peakAnnoList <- lapply(seq_along(macsOutput), function(.x){
    peakAnno <- annotatePeak(macsOutput[[.x]], tssRegion=c(-3000, 1000),
                             TxDb=mmu_TxDb, 
             genomicAnnotationPriority = c("Promoter", "5UTR", 
                                           "Intergenic", "Exon", "Intron",
                                           "3UTR", "Downstream")) 
    pdf(paste0(names(macsOutput)[.x],
     "-Pie charts showing distribution of peaks among genomic features.pdf"), 
        width = 10, height = 5)
    plotAnnoPie(peakAnno)
    vennpie(peakAnno)
    ## some annotation overlap
    upsetplot(peakAnno)
    dev.off()
    peakAnno
})
## Visualize Genomic Annotation

names(peakAnnoList) <- names(macsOutput)

## output

peakAnnoList_df <- lapply(peakAnnoList, as.data.frame)

geneID <- unique(read.delim("/home/hl84w/work/mccb/genome/Mus_musculus_GRCm38.p5/Mus_musculus.GRCm38.101.geneid.name.mapping.txt", header = TRUE, sep = "\t", as.is = TRUE))
colnames(geneID) <- c("geneId", "Symbol", "Biotype")

peakAnnoList_df_namesadded <- lapply(peakAnnoList_df, function(.x){
    .y <- merge(.x, geneID, by = "geneId", all.x = TRUE)
    .y <- .y[, c(colnames(.x), "Symbol", "Biotype")]
    .y
})


library("WriteXLS")
WriteXLS(peakAnnoList_df_namesadded, 
         ExcelFileName = "Annotated.ATAC-seq.consensus.peaks.xlsx", 
         row.names = TRUE, SheetNames = names(peakAnnoList))

pdf("Figure7.Barplots.showing.distribution.of.consensus.peaks.among.genomic.features.pdf", width = 6, height = 2)

plotAnnoBar(peakAnnoList) +
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.title=element_text(size=8),
          text = element_text(size=10),
          legend.key.size = unit(0.3, "cm"),
          legend.key.width = unit(0.3,"cm")) +
    guides(fill=guide_legend(nrow=4,byrow=TRUE))
plotDistToTSS(peakAnnoList,
              title="Distribution consensus peaks\nrelative to TSS") +
    theme(legend.position="bottom", 
          legend.title=element_text(size=8),
          plot.title = element_text(hjust = 0.5, size = 10),
          text = element_text(size=10), 
          legend.key.size = unit(0.3, "cm"),
          legend.key.width = unit(0.3,"cm")) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
dev.off()

save.image(file = "q005.consensus.peaks.ChIPseeker.RData")

### peak overlapping

pdf("Figure8.union.merged.peak.overlap.pdf", width = 3.5, height = 3.5)
p <- makeVennDiagram(macsOutput, totalTest=1e+5,
                         fill= c("#e6194B", "#42d4f4", "#f032e6", "#3cb44b"), # circle fill color
                         col=c("black"), #circle border color
                         cat.col=c("#e6194B", "#42d4f4", "#f032e6", "#3cb44b"))
    
dev.off()


## convert human gene symbol to mouse gene symbol

# Basic function to convert human to mouse gene names
library("biomaRt")
library("readxl")

setwd("C:\\Users\\Haibo\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020")

senescence <- as.data.frame(read_excel("Senescence.assocsiate.genesets.xlsx", 
                         sheet = 1))

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

GeneList_hsp2mmu <- function(x, uniq = TRUE, mart, martL){
    x <- x[!is.na(x)]
    genesV2 <- getLDS(attributes = c("hgnc_symbol"), 
                      filters = "hgnc_symbol", values = x , 
                      mart = mart, attributesL = c("mgi_symbol"),
                      martL = martL, uniqueRows = uniq)
    genesV2
}

# unique_gene_mapping <- lapply(senescence[-2, ], GeneList_hsp2mmu, uniq = TRUE)
# names(unique_gene_mapping) <- colnames(senescence)

one2many_gene_mapping <- lapply(senescence[-1, ], GeneList_hsp2mmu, uniq = FALSE, mart = human, martL = mouse)
names(one2many_gene_mapping) <- colnames(senescence)


## some genes is not mappable
one2many_gene_human <- lapply(one2many_gene_mapping, function(.x){
    .x[, 1]
})

umapped_genes <- mapply(function(.x, .y){
    x_use <- .x[!is.na(.x)]
    x_use[!x_use %in% .y]
}, senescence[-1, ],  one2many_gene_human, SIMPLIFY = TRUE)

save(one2many_gene_mapping, file = "human2mouse.gene.mapping.RData")

one2many_gene_mouse <- lapply(one2many_gene_mapping, function(.x){
    c(.x[, 2], rep(NA, 236 - length(.x[, 2])))
})


## all unique input genes

uniq_human <- unique(as.vector(as.matrix(senescence[-1, ])))
write.table(uniq_human, file = "unique.human.genes.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

## 
install.packages("msigdbr")
library("msigdbr")
msigdbr_show_species()

# Retrieve human genes for all gene sets in the database.
m_df = msigdbr(species = "Mus musculus")
head(m_df)

# Check the available collections and sub-collections.

collections <- m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)


## hallmark
m_df_h = msigdbr(species = "Mus musculus", category = "H")
head(m_df_h)

m_df_h <- as.data.frame(m_df_h)
colnames(senescence)

hallmark <- colnames(senescence)[grepl("^HALLMARK_", colnames(senescence))]

hallmark_mouse <- lapply(hallmark, function(.x){
    mouse_gene <- data.frame(gene_symbol = m_df_h[m_df$gs_name == .x, ]$gene_symbol)
})
names(hallmark_mouse) <- hallmark


## Biocarta
m_df_C2 = msigdbr(species = "Mus musculus", category = "C2",
                 subcategory = "CP:BIOCARTA")
head(m_df_C2)

m_df_C2 <- as.data.frame(m_df_C2)
colnames(senescence)

hallmark <- colnames(senescence)[grepl("^BIOCARTA_", colnames(senescence))]

C2_mouse <- lapply(hallmark, function(.x){
    mouse_gene <- data.frame(gene_symbol = m_df_C2[m_df_C2$gs_name == .x, ]$gene_symbol)
})
names(C2_mouse) <- hallmark

## PID
hallmark <- "PID_IL12_2PATHWAY"

m_df_C2_PID = msigdbr(species = "Mus musculus", category = "C2",
                  subcategory = "CP:PID")
m_df_C2_PID <- as.data.frame(m_df_C2_PID)
C2_PID_mouse <- lapply(hallmark, function(.x){
    mouse_gene <- data.frame(gene_symbol = m_df_C2_PID[m_df_C2_PID$gs_name == .x, ]$gene_symbol)
})
names(C2_PID_mouse) <- hallmark

## REACTOME
hallmark <- colnames(senescence)[grepl("^REACTOME_", colnames(senescence))]

m_df_C2_REACTOME = msigdbr(species = "Mus musculus", category = "C2",
                      subcategory = "CP:REACTOME")
m_df_C2_REACTOME <- as.data.frame(m_df_C2_REACTOME)
C2_REACTOME_mouse <- lapply(hallmark, function(.x){
    mouse_gene <- data.frame(gene_symbol = m_df_C2_REACTOME[m_df_C2_REACTOME$gs_name == .x, ]$gene_symbol)
})
names(C2_REACTOME_mouse) <- hallmark


## CGP
hallmark <- c("FRIDMAN_SENESCENCE_DN", "FRIDMAN_SENESCENCE_UP",
"CHICAS_RB1_TARGETS_SENESCENT")

m_df_C2_CGP = msigdbr(species = "Mus musculus", 
                           category = "C2",
                           subcategory = "CGP")
m_df_C2_CGP <- as.data.frame(m_df_C2_CGP)
C2_CGP_mouse <- lapply(hallmark, function(.x){
    mouse_gene <- data.frame(gene_symbol = m_df_C2_CGP[m_df_C2_CGP$gs_name == .x, ]$gene_symbol)
})
names(C2_CGP_mouse) <- hallmark

test <- c(hallmark_mouse, C2_mouse, C2_PID_mouse, C2_REACTOME_mouse, C2_CGP_mouse)

library("WriteXLS")
WriteXLS(test, 
         ExcelFileName = "MsigDB.mouse.senescence.genesets.xlsx", 
         row.names = FALSE, SheetNames = substr(gsub("_", "", names(test)),1,30))
names(test)


## adding count values to the annotation table

library("readxl")

annotated_peaks <- as.data.frame(read_xlsx("C:\\Users\\Haibo\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\ChIPseeker.out/Annotated.ATAC-seq.consensus.peaks.xlsx"))

count_dir <- "C:\\Users\\Haibo\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\differential_accessibility"
peak_rawcount <- read.delim(file.path(count_dir,"q005.all.samples.merged.peaks.final.refmt.count.table"), as.is = TRUE)[, 3:6]
peak_normalCount <- read.delim(file.path(count_dir, "31 samples using effect model", "Table 1. Qsmooth-Normalized.count.table.txt"), as.is = TRUE)

head(peak_rawcount)
head(peak_normalCount)

peak_normalCount <- merge(peak_rawcount, peak_normalCount, by.x = "Geneid", by.y = "row.names", all.x = TRUE)

annotated_peaks <- merge(annotated_peaks, peak_normalCount, by.x =c("seqnames",	"start",	"end") , by.y = c("Chr", "Start","End" ), all.x = TRUE)

annotated_peaks <- annotated_peaks[, -4]
WriteXLS(annotated_peaks, 
         ExcelFileName = file.path("C:\\Users\\Haibo\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\ChIPseeker.out", "Annotated.ATAC-seq.consensus.peaks.normalized.counts.xlsx"), row.names = FALSE, SheetNames = "")





