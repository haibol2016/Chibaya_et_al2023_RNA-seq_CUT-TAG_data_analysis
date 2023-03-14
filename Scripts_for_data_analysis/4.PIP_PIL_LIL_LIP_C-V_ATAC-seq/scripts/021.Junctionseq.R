
#Load JunctionSeq:
library("JunctionSeq");
#The sample decoder:
decoder <- read.table(
"docs/QoRTs.decoder.txt",
header=T,stringsAsFactors=F);
#The count files:
countFiles <- paste0(
"results/003.QoRTs.out/",
decoder$sample.ID,
"/QC.spliceJunctionAndExonCounts.forJunctionSeq.txt"
);

#Run the analysis:
jscs <- runJunctionSeqAnalyses(
sample.files = countFiles,
sample.names = decoder$sample.ID,
condition = decoder$group.ID,
flat.gff.file = "results/003.QoRTs.out/withNovel.forJunctionSeq.gff.gz",
nCores = 1,
analysis.type = "junctionsAndExons",
verbose=TRUE,
debug.mode = TRUE
);

if (!dir.exists("results/021.JunctionSeq.out/results")
{
   dir.create("results/021.JunctionSeq.out/results", recursive = TRUE)
}

writeCompleteResults(jscs,
outfile.prefix="results/021.JunctionSeq.out/",
save.jscs = TRUE
);


buildAllPlots(
jscs=jscs,
FDR.threshold = 0.01,
outfile.prefix = "results/021.JunctionSeq.out/results",
variance.plot = TRUE,
ma.plot = TRUE,
rawCounts.plot=TRUEverbose = TRUE);
verbose = TRUE);
