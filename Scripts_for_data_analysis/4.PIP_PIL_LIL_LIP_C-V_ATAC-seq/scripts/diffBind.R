#BiocManager::install("ggfortify")
library(ggfortify)
library(DiffBind)


out_dir <- "results/029.DiffBind.out"

if (!dir.exists(out_dir))
{
   dir.create(out_dir)
}

sample_sheet="/home/hl84w/michael_lodato/Narendra/docs/ATAC-seq.sample.sheet.csv"

AVM <- dba(sampleSheet= sample_sheet)

setwd(out_dir)
pdf("Figure.1.raw.peakSet.correlation.pdf", height = 7, width =7)
p <- plot(AVM, attributes= DBA_ID, 
     ColAttributes = c(DBA_TISSUE, DBA_CONDITION),
     colScheme = "YlOrBr",
     cexRow = 0.5,
     cexCol = 0.5,
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
pdf("Figure2.MACS2.peaks.RPKM.correlation.pdf", height = 6, width = 6)
plot(AVM, attributes= DBA_ID, 
     ColAttributes = c(DBA_TISSUE, DBA_CONDITION),
     colScheme = "YlOrBr",
     score = DBA_SCORE_RPKM,
     cexRow = 0.5,
     cexCol = 0.5,
     keysize = 1.5,
     key.title = NULL)
dev.off()
#
pdf("Figure3.MCS2.peaks.PCA.pdf", height = 5, width =5)
dba.plotPCA(AVM, attributes= DBA_FACTOR, 
            score = DBA_SCORE_RPKM,
            cor = TRUE,
            components = 1:2)
dev.off()
#
save.image(file = "MACS2.peaks.RPKM.RData")
#
Score_RPKM <- do.call("cbind", lapply(AVM$peaks, function(.x){
    .x$Score
}))
#
colnames(Score_RPKM) <- paste(AVM$sample$Tissue,AVM$sample$SampleID, sep = "_")

rownames(Score_RPKM) <- paste(AVM$peak[[1]][,1], AVM$peaks[[1]][,2], AVM$peaks[[1]][,3], sep=":")

write.table(Score_RPKM, file = "ATACseq.RPKM.table.txt", sep = "\t", quote = FALSE)
pc <- prcomp(t(Score_RPKM), scale = T, center = T, retx = TRUE)
#
save(pc, file ="MACS2.peaks.RPKM.PCA.RData" )
#
#
samples <- read.csv(sample_sheet)
samples$group <- as.factor(paste(samples$Tissue, samples$Condition, sep = "_"))
#
#
#
pdf(file = "MACS2.peaks.RPKM.PCA.plot.showing.ATAC-seq.samples.relationship.pdf", width = 7, height = 5)
autoplot(pc, data = samples, colour = 'group', size = 2)
dev.off()
#
AVM <- dba.count(AVM, minOverlap=1, 
                 score = DBA_SCORE_READS,
                 fragmentSize = 0,
                 bParallel = TRUE,
                 summits = TRUE)
save(AVM, file = "MACS2.peaks.counts.RData")

Score_count <- do.call("cbind", lapply(AVM$peaks, function(.x){
    .x$Reads
}))
#
colnames(Score_count) <- paste(AVM$sample$Tissue,AVM$sample$SampleID, sep = "_")

rownames(Score_count) <- paste(AVM$peak[[1]][,1], AVM$peaks[[1]][,2], AVM$peaks[[1]][,3], sep=":")

write.table(Score_count, file = "ATACseq.count.table.txt", sep = "\t", quote = FALSE)

