library(dplyr)
library(maftools)
library(stringr)
library(tidyr)
library(tibble)



build_primary_maf_df_hmf <- function() {
    base_path <- "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/mafs/purple_mafs/primary/mafs/"
    primary_df <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/filtered_cleaned_df.csv")
    
    primary_mafs <- list()
    for(i in 1:nrow(primary_df)){
        row <- primary_df[i,]
        sample_id = row[["Sample"]]  
        tumor_id = row[["DNA_primary_cegat"]]
        primary_maf_file <-  paste(base_path, tumor_id, ".purple.maf" , sep="")

        df <- read.csv(primary_maf_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
        
        df <- df %>% filter(
            FILTER =="PASS", 
            (t_alt_count / t_depth) >= 0.05, 
        )

        df$tumor_id <- tumor_id
        df$sample <- sample_id
        df$origin <- "primary"
        primary_mafs[[i]] <- df
        
    }
        
    full_primary_mafs <- bind_rows(primary_mafs)
    return(full_primary_mafs)
}   


