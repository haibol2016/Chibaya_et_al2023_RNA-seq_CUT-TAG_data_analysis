
library("DESeq2")
library("QoRTs");
#Read in the QC data:
res <- read.qc.results.data("results/003.QoRTs.out/",
                         decoder.files = "docs/QoRTs.decoder.txt",
                         calc.DESeq2 = TRUE, calc.edgeR = FALSE);
out <- "results/003.QoRTs.out"

if (!dir.exists(out))
{
  dir.create(out, recursive = TRUE)
}

makeMultiPlot.all(res,
	outfile.dir = out,
	plot.device.name = "pdf");
