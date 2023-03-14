BiocManager::install("tximport")

library(tximport)
dir <- "results/006.salmon.quant.out"
sample <- read.delim("docs/metadata.txt",as.is = TRUE)
colnames(sample) <- c("file_name", "run")
files <- file.path(dir, "salmon", sample$run, "quant.sf")
names(files) <- sample$run

tx2gene <- read.delim("docs/tx2gene.txt",  as.is = TRUE)
colnames(tx2gene) <- c("TxName", "GeneID")
tx2gene$GeneId <- gsub("\\.\\d+", "", tx2gene$GeneId, perl = TRUE)

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)

gene2symbol <- read.delim("~/work/mccb/genome/Mus_musculus_GRCm38.p5/Mus_musculus.GRCm38.101.geneid.name.mapping.txt", as.is = TRUE)
