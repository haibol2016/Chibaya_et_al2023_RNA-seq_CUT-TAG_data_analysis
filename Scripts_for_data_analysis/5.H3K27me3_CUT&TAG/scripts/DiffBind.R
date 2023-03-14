BiocManager::install("ggfortify")
library(ggfortify)
library(DiffBind)


out_dir <- "results/010.DiffBind.out"

if (!dir.exists(out_dir))
{
  dir.create(out_dir)
}

sample_sheet="/home/hl84w/marcus_ruscetti/haibo/cuttag_10272022/docs/H3K27me3.sample.sheet.csv"
AVM <- dba(sampleSheet= sample_sheet)

setwd(out_dir)
pdf("Figure.1.raw.peakSet.correlation.pdf", height = 7, width =7)
p <- plot(AVM, attributes= DBA_ID, 
          ColAttributes = c(DBA_CONDITION),
          colScheme = "YlOrBr",
          cexRow = 1,
          cexCol = 1,
          keysize = 1.5,
          key.title = NULL)
print(p)
dev.off()

AVM <- dba.count(AVM, minOverlap=1, 
                 score = DBA_SCORE_RPKM,
                 fragmentSize = 0,
                 bParallel = TRUE,
                 summits = TRUE)

## correlation heatmap based on the affinity scores  ( 109323 peaks)
#
pdf("Figure2.Genrich.peaks.RPKM.correlation.pdf", height = 6, width = 6)
plot(AVM, attributes= DBA_ID, 
     ColAttributes = c(DBA_CONDITION),
     colScheme = "YlOrBr",
     score = DBA_SCORE_RPKM,
     cexRow = 1,
     cexCol = 1,
     keysize = 1.5,
     key.title = NULL)
dev.off()
#
pdf("Figure3.Genrich.peaks.PCA.pdf", height = 5, width =5)
dba.plotPCA(AVM, attributes= DBA_CONDITION, 
            score = DBA_SCORE_RPKM,
            cor = TRUE,
            components = 1:2)
dev.off()
#
save(AVM, file = "Genrich.peaks.RPKM.RData")
#
Score_RPKM <- do.call("cbind", lapply(AVM$peaks, function(.x){
  .x$Score
}))
#
colnames(Score_RPKM) <- AVM$sample$SampleID

rownames(Score_RPKM) <- paste(AVM$peak[[1]][,1], AVM$peaks[[1]][,2], AVM$peaks[[1]][,3], sep=":")

write.table(Score_RPKM, file = "ATACseq.RPKM.table.txt", sep = "\t", quote = FALSE)
pc <- prcomp(t(Score_RPKM), scale = T, center = T, retx = TRUE)
#
save(pc, file ="Genrich.peaks.RPKM.PCA.RData" )
#
#
samples <- read.csv(sample_sheet)
samples$group <- samples$Condition
#
#
#
pdf(file = "Genrich.peaks.RPKM.PCA.plot.showing.ATAC-seq.samples.relationship.pdf", width = 7, height = 5)
autoplot(pc, data = samples, colour = 'group', size = 2)
dev.off()
#
AVM <- dba(sampleSheet= sample_sheet)
AVM <- dba.count(AVM, minOverlap=1, 
                 score = DBA_SCORE_READS,
                 fragmentSize = 0,
                 bParallel = TRUE,
                 summits = TRUE)
save(AVM, file = "Genrich.peaks.counts.RData")

Score_count <- do.call("cbind", lapply(AVM$peaks, function(.x){
  .x$Reads
}))
#
colnames(Score_count) <- AVM$sample$SampleID

rownames(Score_count) <- paste(AVM$peak[[1]][,1], AVM$peaks[[1]][,2], AVM$peaks[[1]][,3], sep=":")

write.table(Score_count, file = "10272022.CutTag.count.table.txt", sep = "\t", quote = FALSE)
