install.packages("tidyverse")
install.packages("WriteXLS")

library(readxl)
library(WriteXLS)

dir_geneset <- "C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020"

gencode <- read.delim(file.path(dir_geneset, "Gencode.M25.gene.mapping.txt"),
                      header = FALSE, as.is = TRUE)
senescence_gene <- read_excel(file.path(dir_geneset, 
                                    "all.mouse.senescence.genelist.xlsx"), 
                          sheet = 1, guess_max = 21474836)
head(gencode)
head(senescence_gene)

## all senescence genes are in mouse genecode naming
lapply(senescence_gene, function(.x){
    sum(!is.na(.x)) == sum(!is.na(.x) & .x %in% gencode$V4)
})

## how many uniq genes in the lists
uniq_senescence <-  as.vector(as.matrix(senescence_gene))
head(uniq_senescence)
uniq_senescence <- unique(uniq_senescence[!is.na(uniq_senescence)]) 
length(uniq_senescence)  ## 2174

gencode_senescence <- gencode[gencode[, 4] %in% uniq_senescence, ]
check_1 <- gencode_senescence[gencode_senescence[, 4] %in% gencode_senescence[, 4][duplicated(gencode_senescence[, 4])], ]
length(unique(check_1[,4]))  ## 27 gene symbol with two Ensembl IDs

## generate a data frame for the union of 2201 senescence genes
dumpy_df <- as.data.frame(matrix("", nrow = 2201, 
                                  ncol = ncol(senescence_gene)))
rownames(dumpy_df) <- gencode_senescence[,1]
colnames(dumpy_df) <- colnames(senescence_gene)
dumpy_df$GeneSymbol <- gencode_senescence[,4]

lapply(seq_along(dumpy_df[, 1:26]), function(.x)
    {
        dumpy_df[,.x] <<- ifelse(dumpy_df$GeneSymbol %in% 
                                    senescence_gene[, .x][!is.na(senescence_gene[, .x])], 
                                 dumpy_df$GeneSymbol, "")
    })


## adding to the peak annotation
dir_peaks <- paste0("C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus",
                    "\\ATAC-seq-09112020\\ChIPseeker.out")

peaks <- read_excel(file.path(dir_peaks, 
                              "Annotated.ATAC-seq.consensus.peaks.normalized.counts.xlsx"),
                    sheet = 1, guess_max = 21474836)

peaks <- merge(peaks, dumpy_df, by.x = "geneId", 
               by.y = "row.names", sort = FALSE, all.x = TRUE)

peaks <- peaks[, c(2:13, 1, 14:75)]
peaks$geneChr <- peaks$seqnames

## number of genes within each set
gene_openCR <- lapply(peaks[,50:75], function(.x)
    {
        unique(.x[!is.na(.x) & .x !=""])
})


WriteXLS(peaks, ExcelFileName = file.path(dir_peaks, 
     "Senesce.geneset.annotated.ATAC-seq.consensus.peaks.normalized.counts.xlsx"))

write.table(dumpy_df, file.path(dir_geneset,
                                "All.mouse.senescence.geneset.refmt.txt"),
            row.names = TRUE, sep = "\t", quote = TRUE)

## Annotate DAR regions
dir_dar <- paste0("C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work",
                 "\\Marcus\\ATAC-seq-09112020", 
                 "\\differential_accessibility\\31 samples using effect model")
xlsx_files <- dir(dir_dar, "xlsx$")

out_dir <- paste0("C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work", 
                  "\\Marcus\\ATAC-seq-09112020", 
                  "\\differential_accessibility/Annotate.DAR.out")
if (!dir.exists(out_dir))
{
    dir.create(out_dir)
}

null <- lapply(xlsx_files, function(.x) {
    temp_xlsx <- read_excel(file.path(dir_dar,.x), sheet = 1, guess_max = 21474836)
    temp_xlsx <- merge(temp_xlsx, peaks, 
                       by.x= c("Chr","Start", "End"), 
                       by.y = c("seqnames", "start", "end"),
                       all.x = TRUE, sort = FALSE)
    WriteXLS(temp_xlsx, ExcelFileName = file.path(out_dir, .x))
})

## Annotate DAR regions of simple effect
dir_dar <- paste0("C:\\Users\\liuh\\Desktop\\UMass_work", "\\Marcus\\ATAC-seq-09112020", 
                  "\\differential_accessibility\\31 samples using effect model")

xlsx_files <- dir(dir_dar, "xlsx")

null <- lapply(xlsx_files, function(.x) {
    temp_xlsx <- read_excel(file.path(dir_dar, .x), sheet = 1, guess_max = 21474836)
    temp_xlsx <- merge(temp_xlsx, peaks, 
                       by.x= c("Chr","Start", "End"), 
                       by.y = c("seqnames", "start", "end"),
                       all.x = TRUE, sort = FALSE)
    temp_xlsx$geneChr <- temp_xlsx$Chr
    temp_xlsx$Geneid.y <- NULL
    WriteXLS(temp_xlsx, ExcelFileName = file.path(out_dir, .x))
})


## over-representation analysis
all_peak_genes <- unique(peaks$Symbol)

## 1967 of 2201 genes involved
dumpy_df_bkg <- dumpy_df[rownames(dumpy_df) %in% all_peak_genes, ]

dar_peaks_dir <- paste0("C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work",
                        "\\Marcus\\ATAC-seq-09112020", 
                        "\\differential_accessibility\\Annotate.DAR.out")
xlsx_files_in <- dir(dar_peaks_dir, ".xlsx$")

library("tidyverse")
library("broom")
null <- lapply(xlsx_files_in, function(.x)
    {
        temp <- as.data.frame(read_excel(file.path(dar_peaks_dir, .x), 
                                         sheet =1, guess_max = 21474836))
        temp_sig <- temp[temp$log2FoldChange <= -log2(2), 
                         c(colnames(dumpy_df)[1:26], "Symbol")]
        dar_gene <- unique(temp_sig$Symbol[!is.na(temp_sig$Symbol)])
        
        ft_result <- lapply(seq_along(colnames(dumpy_df)[1:26]), function(.x){
            gene_set <- dumpy_df[, .x][dumpy_df[, .x]!=""]
            dar_geneset <- unique(temp_sig[, .x][!is.na(temp_sig[,.x])])
            
            test_mtx <- matrix(c(length(dar_geneset), 
                                 sum(!gene_set %in% dar_geneset),
                                 sum(!dar_gene %in% dar_geneset),
                                 sum(!all_peak_genes %in% 
                                         c(dar_gene, gene_set[!gene_set %in% dar_geneset]))),
                               nrow = 2, byrow = TRUE)
            #print(.x)
            #print(test_mtx)
            ft <- fisher.test(test_mtx, alternative = "two.sided")
            ft <- as.data.frame(tidy(ft))
            ft$GenesetName <- colnames(dumpy_df)[.x]
            ft$Genes <- paste(dar_geneset, collapse = ",")
            #ft$EnrichmentFold <- ifelse(length(dar_geneset) == 0, 0, length(dar_gene)/length(gene_set))
            ft
        })
        results <- do.call("rbind", ft_result)[, c(7,8,1,2)]
        results$adj.p <- p.adjust(results$p.value, method = "BH")
        results
})

sheetNames <-gsub("shrunken", "shrnk", gsub(" |\\.", "", 
                                            gsub(".LFC.xlsx|.DAR","", 
                                                 xlsx_files_in)))



out_dir <- paste0("C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work", 
                  "\\Marcus\\ATAC-seq-09112020", 
                  "\\differential_accessibility")
if (!dir.exists(out_dir))
{
    dir.create(out_dir)
}
WriteXLS(null, 
         ExcelFileName = file.path(out_dir,
                                   "0.ORA.down.DAR.enriched.genesets.two-sided.xlsx"),
         SheetNames = sheetNames)


## generate genomic regions for HOMER motif analysis
library("readxl")
dar_peaks_dir <- "/home/hl84w/marcus_ruscetti/results/013.DAR.out"
files <- dir(dar_peaks_dir, ".xlsx$", full.names = TRUE)

out_names <- gsub("\\s+", "", gsub(".+/(.+?).xlsx", "\\1.bed", 
                                   files, perl = TRUE), perl = TRUE)
out_dir <- "/home/hl84w/marcus_ruscetti/results/014.DAR.BED.out"
if (!dir.exists(out_dir))
{
    dir.create(out_dir)
}

for (i in seq_along(files))
{
    temp <- as.data.frame(read_excel(files[i], sheet = 1, 
                                     guess_max = 21474836))
    
    ## only log2(FC) > log2(1.5)
    temp <- unique(temp [abs(temp$log2FoldChange) >= log2(1.5), c(4:6, 1)])
    temp$score <- "."
    temp$strand <- "+"
    write.table(temp, file.path(out_dir, out_names[i]), 
                sep = "\t", quote = FALSE, row.names =FALSE,
                col.names = FALSE)
}






