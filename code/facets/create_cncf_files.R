library(dplyr)


files <- list.files("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/facets_redux", pattern = "*_hisens\\.rds", recursive=TRUE, full.names=TRUE)

for(file in files) {
    df <- readRDS(file)
    
    target_dir <- dirname(file)
    output_path <- file.path(target_dir, "cncf.csv")
    write.csv(df$segs, file = output_path, row.names = FALSE)
}
