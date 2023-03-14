BiocManager::install("clusterProfiler")
BiocManager::install("ReactomePA")
BiocManager::install("enrichplot")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("DOSE")

## Disease Ontology (DO) Semantic and Enrichment analysis.
library("DOSE")
library("clusterProfiler")
library("ReactomePA")
library("enrichplot")
library("org.Mm.eg.db")
library("readxl")

## universe 

universe_dir <- "C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\ChIPseeker.out"

universe <- read_excel(file.path(universe_dir, "Senesce.geneset.annotated.ATAC-seq.consensus.peaks.normalized.counts.xlsx"),
                       col_names = TRUE, guess_max = 150000)
universe <- as.data.frame(universe[, "geneId"])
universe <- unique(universe[,1])

columns(org.Mm.eg.db)
keytypes(org.Mm.eg.db)
keys(org.Mm.eg.db, keytype = "ENSEMBL")
mapIds(x, keys, column, keytype, ..., multiVals)

egids <- select(org.Mm.eg.db, keys = unique(universe), 
       columns = "ENTREZID", keytype ="ENSEMBL")

universe <- egids[,2][!is.na(egids[,2])]

## test gene lists

test_dir <- "C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\differential_accessibility\\Annotate.DAR.out"

setwd(test_dir)
files <- dir(test_dir, "shrunken.LFC.xlsx$", full.names = TRUE)
names(files) <- gsub(" ", "", gsub(".+/(.+?).DAR.shrunken.LFC.xlsx",
                                   "\\1", files, perl = TRUE))
out_dir <- "C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\differential_accessibility/GO.pathway.analysis"

if(!dir.exists(out_dir)){
    dir.create(out_dir)
}

load("C:\\Users\\liuh\\OneDrive - University Of Massachusetts Medical School\\Desktop\\UMass_work\\Marcus\\ATAC-seq-09112020\\differential_accessibility/GO.pathway.analysis/GO.KEGG.pathway.RData")
null <- mapply(function(.x, .y){
    #.x <- files[1]
    temp <- read_excel(.x, col_names =  TRUE, guess_max = 150000)
    temp <- as.data.frame(unique(temp[temp$log2FoldChange <= -1, 
                                         "geneId"]))
    temp_egids <- select(org.Mm.eg.db, keys = unique(temp[,1]), 
                         columns = "ENTREZID", keytype ="ENSEMBL")
    gene <- temp_egids[, 2][!is.na(temp_egids[,2])]
    
    ## Gene ontology terms
    for (GO_ont in c("BP", "MF", "CC")){
        ego <- enrichGO(gene  = gene,
                        keyType = "ENTREZID",
                        universe      = universe,
                        OrgDb         = org.Mm.eg.db,
                        ont           = GO_ont,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.1,
                        minGSSize = 5,
                        maxGSSize = 500,
                        readable      = TRUE)
        df <- as.data.frame(ego)
        if (nrow(df) > 0) {
            write.table(df, file=file.path(out_dir, 
                                           paste("Table 1",.y, GO_ont, "down_DAR.txt")),
                        sep = "\t", quote = FALSE)
            pdf(file=file.path(out_dir, 
                               paste("Fig 1.", .y, GO_ont, "down_DAR.pdf")),
                width = 8, height = 3.5 + nrow(df)*5/72)
            print(dotplot(ego, showCategory=30))
            dev.off()
        } else {
            print("No enrichment")
        }
    }
    
    ## Enrichment Reactome pathways
    reactome_pathway <- enrichPathway(gene,   # entrez gene id
                                      organism = "mouse", 
                                      pvalueCutoff = 0.05, 
                                      pAdjustMethod = "BH", 
                                      qvalueCutoff = 0.1, 
                                      universe = "", 
                                      minGSSize = 5, 
                                      maxGSSize = 500, 
                                      readable =  TRUE)
    df <- as.data.frame(reactome_pathway)
    if (nrow(df) > 0) {
        write.table(df, file=file.path(out_dir, 
                                       paste("Table 2.",.y, ".Reactome.down_DAR.txt")),
                    sep = "\t", quote = FALSE)
        pdf(file=file.path(out_dir, 
                           paste("Fig 2.", .y, ".Reactome.down_DAR.pdf")),
            width = 8, height = 3.5 + nrow(df)*5/72)
        print(dotplot(reactome_pathway, showCategory=30))
        dev.off()
    }
    
    ## KEGG pathway
    kegg_pathway <- enrichKEGG(gene = gene, #entrez gene id
                               organism     = 'mmu',
                               keyType = "kegg",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               universe = universe,
                               minGSSize = 5,
                               maxGSSize = 500,
                               qvalueCutoff = 0.1,
                               use_internal_data = FALSE)
    df <- as.data.frame(kegg_pathway)
    if (nrow(df) > 0){
        write.table(df, file=file.path(out_dir, 
                                       paste("Table 3.", .y, ".KEGG.down_DAR.txt")),
                    sep = "\t", quote = FALSE)
        pdf(file=file.path(out_dir, 
                           paste("Fig 3.", .y, ".KEGG.down_DAR.pdf")),
            width = 8, height = 3.5 + nrow(df)*5/72)
        print(dotplot(kegg_pathway, showCategory=30))
        dev.off()
    }
    
}, files, names(files), SIMPLIFY = FALSE)
