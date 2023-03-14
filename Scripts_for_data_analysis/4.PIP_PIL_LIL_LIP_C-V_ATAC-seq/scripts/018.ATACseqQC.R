#!/usr/bin/env Rscript
## This script accept one preprocessed BAM file for QC analysis. An index file per BAM file must be included in the same directory.
## Loading all required packages
library(motifStack)
library(ATACseqQC)

## load the BSgenome package
library("BSgenome.Mmusculum.ENSEMBL.mm10")
## load TxDb package/database
library("GenomicFeatures")
library("DescTools")
library(MotifDb)
library(ChIPpeakAnno)
## Getting the BAM file name and sample ID
args <- commandArgs(trailingOnly = TRUE)
bamfile <- args[1]
bamfile.sample.ID <- gsub(".bam", "", basename(bamfile))
outPath <- file.path("results", "018.ATACseqQC.out",
                     paste0(bamfile.sample.ID, ".splited.bam"))
if (!dir.exists(outPath)) {
    dir.create(outPath, recursive = TRUE)
}
# Plotting size distribution of ATAC-seq library insert (Figure 1G)
#pdf(file.path(
#    outPath,
#    paste0(bamfile.sample.ID, ".fragment.size.distribution.pdf")),
#width = 10,
#height = 8)
#fragSize <- fragSizeDist(bamFiles = bamfile, 
#                         bamFiles.labels = bamfile.sample.ID)
#dev.off()

## BAM file tags to be included when read in by the Rsamtools
tags <- c("AS", "NM", "MD")

## Build GRangs for the genome excluding unplaced scaffoldsand chrY
genome <- read.delim("docs/mouse.auto.x.chrom.size.txt", header = F)
seqlev <- as.character(genome[, 1])
gr <- GRanges(as.character(genome[, 1]),
              IRanges(genome[, 2], genome[, 3]))
## For QC, read alginments from chromosomes 1 and 2 are representatitive enough

# which <- gr[seqnames(gr) %in% c("chr1", "chr2")]
which <- gr[seqnames(gr) %in% c("1", "2")]
## Reading in paired end read alignments
gal <- readBamFile(bamfile,
                   tag = tags,
                   which = which,
                   asMates = TRUE)
## Shifting the coordinates of 5' ends of the aligned reads in the bam file, +4 for reads mapping
## to the positive strand, -5 for reads mapping to the negative strand
gal1 <- shiftGAlignmentsList(gal)
shiftedBamfile <- file.path(outPath, paste0(bamfile.sample.ID, ".shifted.bam"))
## Outputting the shifted BAM file
export(gal1, shiftedBamfile)
## Getting information of known transcripts from UCSC genome browser database
## or build TxDb from GTF/GFF
# txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb  <- loadDb("/home/hl84w/work/mccb/genome/Mus_musculus_GRCm38.p5/Mus_musculus.GRCm38.p6.v101.sqlite")
txs <- transcripts(txdb)

## TSSenrichment

tsse <- TSSEscore(gal1, txs)

pdf(file = file.path(outPath, paste0(bamfile.sample.ID,
                                      ".TSSenrichment.pdf")),
   width = 4, height = 4)

x <- 100*(-9:10-.5)
y <- tsse$values
plot(x, y, type="b", 
     xlab="Distance to TSS",
     ylab="Aggregated TSS score")
lines(loess(y~x), col="red")
dev.off()

## Classifying reads into nucleosome-free, mono-, di- and tri-nucleosome
## bins based on their fragment sizes. Don't use the machine learning option. TOO SLOW.
# genome <- Hsapiens

genome <- Mmusculum

objs <- splitGAlignmentsByCut(gal1, txs = txs)
## output split BAM files
null <- writeListOfGAlignments(objs, outPath)
## Heatmap and coverage curve for nucleosome-free amd mono-, di- and tri-nucleosome occupied regions
bamfiles <- file.path(
    outPath,
    c(
        "NucleosomeFree.bam",
        "mononucleosome.bam",
        "dinucleosome.bam",
        "trinucleosome.bam"))
## Extracting TSSs coordinates
TSS <- promoters(txs, upstream = 0, downstream = 1)
TSS <- unique(TSS)
## Estimating the library size for normalization
librarySize <- estLibSize(bamfiles)
## Calculating the signals around TSSs.
NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(
    bamfiles,
    TSS = TSS,
    librarySize = librarySize,
    seqlev = seqlev,
    TSS.filter = 0.5,
    n.tile = NTILE,
    upstream = ups,
    downstream = dws
)
## log2 transformed signals
names(sigs) <- gsub(".bam", "", basename(names(sigs)))
sigs.log2 <- lapply(sigs, function(.ele)
    log2(.ele + 1))
## Plotting heatmap showing signals  for nucleosome-free and oligonucleosome-bound
## regions around TSSs.  (Figure 1 H and 1I)
pdf(file = file.path(outPath, paste0(bamfile.sample.ID, 
                                      ".heatmap and averaged coverage.pdf")))
featureAlignedHeatmap(sigs.log2,
                      reCenterPeaks(TSS, width = ups + dws),
                      zeroAt = .5,
                      n.tile = NTILE)
out <- featureAlignedDistribution(
    sigs,
    reCenterPeaks(TSS, width = ups + dws),
    zeroAt = .5,
    n.tile = NTILE,
    type = "l",
    ylab = "Averaged coverage"
)

## rescale the nucleosome-free and nucleosome signals to 0~1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
matplot(out, type="l", xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()
## Plotting CTCF footprints. (Figure 2C)
CTCF <- query(MotifDb, c("CTCF"))
CTCF <- as.list(CTCF)
print(CTCF[[1]], digits = 2)
# seqlev <- c("chr1", "chr2")

seqlev <- c("1", "2")

nucleosome_free_bamfile <- file.path(outPath, "NucleosomeFree.bam")

pdf(file = file.path(outPath, paste0(bamfile.sample.ID, ".CTCF.footprint.pdf")))
factorFootprints(
    nucleosome_free_bamfile,
    pfm = CTCF[[1]],
    genome = genome,
    min.score = "95%",
    seqlev = seqlev,
    upstream = 100,
    downstream = 100
)
dev.off()
save.image(file = file.path(outPath, paste0(bamfile.sample.ID, ".RData")))

