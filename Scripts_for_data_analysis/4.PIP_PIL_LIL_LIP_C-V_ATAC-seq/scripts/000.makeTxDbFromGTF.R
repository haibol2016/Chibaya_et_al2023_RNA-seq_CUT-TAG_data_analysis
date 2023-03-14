
library("GenomicFeatures")
library("GenomeInfoDb")

gtf <- "/home/hl84w/work/mccb/genome/Mus_musculus_GRCm38.p5/Mus_musculus.GRCm38.101.gtf"

mouse_genome <- read.delim("docs/mouse.auto.x.chrom.size.txt", header = F)
seq_info <- Seqinfo(seqnames = mouse_genome[, 1], 
                    seqlengths= mouse_genome[, 3],
                    isCircular= ifelse(mouse_genome[, 1] == "MT", TRUE, FALSE),
                    genome= "GRCm38.p6")

mmu_TxDb <-  makeTxDbFromGFF(file = gtf,
         format= "gtf",
         dataSource = "Ensembl",
         organism= "Mus musculus",
         taxonomyId= 10090,
         circ_seqs= "MT",
         chrominfo=seq_info,
         miRBaseBuild=NA)

saveDb(mmu_TxDb, file="/home/hl84w/work/mccb/genome/Mus_musculus_GRCm38.p5/Mus_musculus.GRCm38.p6.v101.sqlite")
