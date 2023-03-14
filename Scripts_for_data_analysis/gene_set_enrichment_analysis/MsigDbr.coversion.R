library(msigdbr)
library("dplyr")

out_dir <- "."
setwd(out_dir)

species <-  msigdbr_species()

species <- "Mus musculus"
null <- lapply(species, function(sp){
    
    # sub directory for each species
    out_subdir <- file.path(out_dir,gsub("\\s+", "_", sp, perl = TRUE))
    if (!dir.exists(out_subdir)){
        dir.create(out_subdir)
    }
    
    all_gene_sets <- msigdbr(species = sp) %>% 
        select(gs_cat, gs_subcat, 
               gs_name, entrez_gene,
               gene_symbol)
    
    all_gene_sets$gs_subcat <- paste(all_gene_sets$gs_cat,
                                     all_gene_sets$gs_subcat,
                                     sep = ":")
    all_gene_sets$gs_subcat <- gsub(":$", "", all_gene_sets$gs_subcat, 
                                    perl = TRUE)
    all_gene_sets <- as.data.frame(all_gene_sets)
    genesets_list <- split(all_gene_sets,
                           f = all_gene_sets$gs_subcat)
    for (i in seq_along(genesets_list))
    {
        id_types <- c("entrez_gene", "gene_symbol")
        output_subdir <- file.path(out_subdir, id_types)
        
        null <- lapply(output_subdir, function(.x){
            if (!dir.exists(.x)){
                dir.create(.x)
            }
        })
        
        for (j in seq_along(id_types))
        {
            fn <- file.path(output_subdir[j],
                            paste0(gsub(":", "-", names(genesets_list)[i], 
                                        perl = TRUE), ".gmt"))
            ## remove existing files
            if (file.exists(fn))
            {
                unlink(fn)
            }
            
            temp_subcat <- genesets_list[[i]][, c("gs_name", id_types[j])]
            temp_subcat_list <- split(temp_subcat, f= temp_subcat$gs_name)
            for (k in seq_along(temp_subcat_list))
            {
                cat(names(temp_subcat_list)[k], "na", 
                    temp_subcat_list[[k]][, id_types[j]],
                    file = fn, sep = "\t",
                    append = TRUE)
                cat("\n",  file = fn, sep = "\t", append = TRUE)
            }
        }
    }
})

