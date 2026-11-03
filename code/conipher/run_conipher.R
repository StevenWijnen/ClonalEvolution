library(dplyr)
library(CONIPHER)
library(parallel)


files <- list.files("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/conipher_tables_full_oncoanalyser")
base_path <- "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/conipher_tables_full_oncoanalyser"
out_dir <- "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/CONIPHER_oncoanalyser_no_cleaning"

run_conipher <- function(filename) {
  case_id <- sub("\\.tsv$", "", filename)
  input <- file.path(base_path, filename)
  case_out_dir <- file.path(out_dir, case_id)
  
  message("Starting: ", case_id)
  
  # Clustering
  clustering_ok <- tryCatch({
    conipher_clustering(
      case_id = case_id,
      out_dir = case_out_dir,
      input_tsv_loc = input,
      min_cluster_size = 4,
      clonal_cutOff = 0.9,
      # clean_clusters = TRUE,
      nProcs = 4         
    )
    TRUE
  }, error = function(e) {
    message("❌ Clustering failed for ", case_id, ": ", conditionMessage(e))
    FALSE
  })
  
  if (!clustering_ok) return(list(case_id = case_id, status = "clustering_failed"))
  
  # Tree building
  tree_building_input <- file.path(out_dir, case_id, paste0(case_id, ".SCoutput.CLEAN.tsv"))
  
  tryCatch({
    conipher_treebuilding(
      input_tsv_loc = tree_building_input,
      out_dir = case_out_dir,
      prefix = case_id,
      pval_cutoff = 0.05, #Higher p values cause clusters to merge easier (p value for testing, are these the same clusters as H0)
      min_cluster_size = 4,
      ccf_buffer = 15
    )
    list(case_id = case_id, status = "success")
  }, error = function(e) {
    message("❌ Treebuilding failed for ", case_id, ": ", conditionMessage(e))
    list(case_id = case_id, status = "treebuilding_failed")
  })
}

n_parallel <- 32

results <- mclapply(files, run_conipher, mc.cores = n_parallel)

# Summarise results
results_df <- bind_rows(results)
failed <- results_df %>% filter(status != "success")

print("Done. Failed cases:")
print(failed)


