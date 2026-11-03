library(dplyr)
library(stringr)
library(tidyr)
library(tibble)


build_ln_met_maf_df_hmf <- function() {

    base_path <- "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/mafs/purple_mafs/ln_met/mafs/"
    ln_met_df = read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/filtered_cleaned_df.csv")
    ln_mafs <- list()
    for(i in 1:nrow(ln_met_df)){
        row <- ln_met_df[i,]
        tumor_id = row[["DNA_ln_cegat"]]
        sample_id = row[["Sample"]]  

        ln_met_maf_file <-  paste(base_path, tumor_id, ".purple.maf" , sep="")

        df <- read.csv(ln_met_maf_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#")
        
        df <- df %>% filter(
            FILTER =="PASS", 
            (t_alt_count / t_depth) >= 0.05, 
        )

        df$tumor_id <- tumor_id
        df$sample <- sample_id
        df$origin <- "ln_met"
                ln_mafs[[i]] <- df
        
    }
        
    full_ln_met_mafs <- bind_rows(ln_mafs)
    return(full_ln_met_mafs)
}   

