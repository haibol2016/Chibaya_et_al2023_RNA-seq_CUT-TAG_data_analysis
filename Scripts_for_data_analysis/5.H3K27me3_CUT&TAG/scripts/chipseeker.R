
# conda activate seurat4
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
peak_files <- dir(".", ".6col.bed$", full.names = TRUE)

macsOutput <- lapply(peak_files, function(.x){
  toGRanges(.x, format = "BED")
}) 

names(macsOutput) <- c("consensusPeaks.q05")

## genome-wide coverage plot
pdf("Figure1.Genome-wide H3K27me3 peak distribution.pdf", width =12, height = 10)
null <- lapply(seq_along(macsOutput), function(.x){
  print(covplot(macsOutput[[.x]], weightCol = "score",
                title = names(macsOutput)[[.x]],
                chrs = paste0(c(as.character(1:19), "X", "Y"))))
})
dev.off()

## Profile of ChIP peaks binding to TSS regions
pdf("Figure2.consensus.peak.tagHeatmap.2.TSS.pdf", width = 8, height =6)

promoter <- getPromoters(TxDb=mmu_TxDb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(macsOutput, getTagMatrix, windows=promoter)
names(tagMatrixList) <- names(macsOutput)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color="blue")

dev.off()

pdf("Figure3.consensus.peak.AvgProfile.2.TSS.pdf", width = 6, height = 3)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), facet="none") + theme(legend.position="bottom")
dev.off()


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
         ExcelFileName = "Annotated.H3K27me3.q.0.05.consensus.peaks.xlsx", 
         row.names = TRUE, SheetNames = names(peakAnnoList))

write.table(peakAnnoList_df_namesadded[[1]], 
         file = "Annotated.H3K27me3.q.0.05.consensus.peaks.txt",
         sep = "\t",
         row.names = FALSE, col.names = TRUE, quote = TRUE)

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

save.image(file = "q05.consensus.peaks.ChIPseeker.RData")


