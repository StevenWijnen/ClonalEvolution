library(dplyr)

files <- list.files("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/facets_redux", pattern = "*_hisens\\.rds", recursive=TRUE, full.names=TRUE)
clinical_data <- read.csv()

for(file in files) {
    df <- readRDS(file)
    str_split
    output_path <- file.path(target_dir, "cncf.csv")
}
