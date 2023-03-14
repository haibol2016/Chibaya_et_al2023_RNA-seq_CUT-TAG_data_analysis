library("ggplot2")

setwd("C:\\Users\\Haibo\\Desktop\\UMASS work\\Marcus\\Enrichr.out")
folders <- dir(full.names = TRUE)

out <- "visualization"
if (!dir.exists(out)){
    dir.create(out)
}

for (f in folders){
    files <- dir(f, full.names = TRUE)
    out_names <- gsub(".+/(.+)", "\\1", f, perl = TRUE)
    set_names <- gsub(".+/(.+?)_table.txt", "\\1", files, perl = TRUE)
    tf <- mapply(function(.x, .y){
        tf_enrichment <- read.delim(.x, header = TRUE, as.is = TRUE)[, c(1, 4, 8)]
        tf_enrichment <- tf_enrichment[tf_enrichment$Adjusted.P.value < 0.05, ]
        if (nrow(tf_enrichment) > 0) {
            tf_enrichment$set <- .y
        }
        tf_enrichment
    }, files, set_names, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    tf <- do.call("rbind", tf)
    tf <- tf[order(tf$set, tf$Combined.Score), ]
    tf$Term <- factor(tf$Term, levels = tf$Term)
    
    height = nrow(tf) * 0.15
    if (height < 3){
        height = 2.5
    }
    pdf(file.path(out, paste0(out_names, "enriched.transcription.factors.pdf")), height = height, width = 10)
    p <- ggplot(tf, aes(Term, Combined.Score, fill = set)) + 
        geom_bar(stat = "identity",  width= 0.2, 
                 position = position_dodge(width=0.2)) + 
        coord_flip() + 
        ylab("Combined Score") + theme_bw() +
        xlab("Term") + 
        guides(fill=guide_legend(title="Data Source", 
                                 nrow = length(unique(tf$set)))) +
        theme(legend.position="top", 
              axis.text = element_text(size = 8))
        
    print(p)
    dev.off()
}