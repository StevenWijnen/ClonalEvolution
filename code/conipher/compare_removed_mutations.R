# Load tidyr for reshaping the dataframe at the end
# install.packages("tidyr") # run this once if you don't have it installed
library(tidyr)

base_path <- "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/"
options <- c(
    "3" = "CONIPHER_w_BAM_adjusted_and_tcn_adjusted_3muts_standard",
    "4" = "CONIPHER_w_BAM_adjusted_and_tcn_adjusted_4muts",
    "5" = "CONIPHER_w_BAM_adjusted_and_tcn_adjusted_standard_settings",
    "7" = "CONIPHER_w_BAM_adjusted_and_tcn_adjusted_0.05_7muts"
)

full_paths <- setNames(paste0(base_path, options), names(options))

# Initialize an empty list to store row data dynamically
results_list <- list()

for (opt_name in names(full_paths)) {
    dir <- full_paths[[opt_name]]
    print(paste("Processing directory:", dir))
    
    files <- list.files(
        path = dir, 
        pattern = "SCoutput\\.DIRTY\\.tsv$", 
        full.names = TRUE, 
        recursive = TRUE
    )
    
    if (length(files) > 0) {
        for (file_path in files) {
            # Read file and count rows
            data <- read.delim(file_path, header = TRUE)
            n_rows <- nrow(data)
            
            # Extract the sample name from the filename
            # E.g., turns "PatientA.SCoutput.DIRTY.tsv" into "PatientA"
            filename <- basename(file_path)
            sample_name <- gsub("\\.SCoutput\\.DIRTY\\.tsv$", "", filename)
            
            # Create a 1-row dataframe for this specific file
            temp_df <- data.frame(
                Sample = sample_name,
                Setting = opt_name,
                Removed_Muts = n_rows,
                stringsAsFactors = FALSE
            )
            
            # Append it to our results list
            results_list[[length(results_list) + 1]] <- temp_df
        }
    } else {
        print(paste("-> Warning: No matching files found in setting", opt_name))
    }
}

# Combine all the individual rows into one massive long dataframe
long_df <- do.call(rbind, results_list)

# Reshape the dataframe so each sample is a row and settings are columns
final_df <- pivot_wider(
    data = long_df, 
    names_from = Setting, 
    values_from = Removed_Muts,
    names_prefix = "Setting_" # Makes columns "Setting_3", "Setting_4", etc.
)

# View the final dataframe
print(final_df)