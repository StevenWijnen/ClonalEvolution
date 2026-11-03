library(dplyr)
library(facets)
library(tidyr)
library(facetsSuite)
library(GenomicRanges)
library(readr)

source("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/code/tree_analyses/load_tree_struct.R")

base_path = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/facets_redux/"
facets_primary_path <- paste0(base_path, "primary/results/")
facets_ln_met_path <- paste0(base_path, "ln_met/results/")

clinical_data <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/filtered_cleaned_df.csv")
#Subset to the 26 samples which we used for tree analyses

good_samples <- get_topologies_oncoanalyser() %>% pull(case_id)

clinical_data <- clinical_data %>% filter(Sample %in% good_samples)
    
process_facets_file <- function(path, case_id, id, site) {
    
    facets_file <- readRDS(path)
    ssegs <- facets_file$segs
    ssegs <- ssegs %>% filter(nhet >= 10, !is.na(mafR))
    
    facets_file$segs <- ssegs
    primary_cna_df <- facetsSuite::gene_level_changes(facets_file, genome = "hg38")
    primary_cna_df$ploidy <- facets_file$ploidy
    primary_cna_df$purity <- facets_file$purity
    primary_cna_df$Tumor_Sample_Barcode <- id
    primary_cna_df$Sample <- case_id
    
    filtered_facets_df <- primary_cna_df %>% 
    filter(filter %in% c("PASS", "suppress_segment_too_large")) %>% 
    filter(cn_state != "DIPLOID") %>%
    mutate(status = case_when(
            median_cnlr_seg > 0.1 & tcn.em - ploidy >= 0.3 ~ "Amp",
            median_cnlr_seg < -0.1 & tcn.em - ploidy <= -0.3 ~ "Del",
            TRUE ~ "Neutral"
    )) %>% 
    filter(chrom != 23) %>% 
    dplyr::select(
        gene,
        Tumor_Sample_Barcode,
        Sample,
        status, 
        chrom,
        median_cnlr_seg, 
        ploidy,
        purity, 
        cn_state, 
        filter,
        tcn.em,
        lcn.em,
        cf.em
        # Start_Position = gene_start,
        # End_Position = gene_end,
    )
    
    return(filtered_facets_df)
}


process_patient <- function(p_id, normal_id, primary_tumor_id, ln_tumor_id) {

    primary_rds_path <- paste(base_path, "primary/results", primary_tumor_id, 
                            paste0(primary_tumor_id, "_", normal_id, "_hisens.rds"), sep="/")
    ln_rds_path <- paste(base_path, "ln_met/results", ln_tumor_id, 
                       paste0(ln_tumor_id, "_", normal_id, "_hisens.rds"), sep="/")

    primary_df <- process_facets_file(primary_rds_path, p_id, primary_tumor_id, "primary")
    ln_df      <- process_facets_file(ln_rds_path, p_id, ln_tumor_id, "ln_met")

    ## full join on gene — keep every gene seen as altered in EITHER site
    comparison <- full_join(
        primary_df %>% select(gene, chrom, status_primary = status, 
                               tcn_primary = tcn.em, cf_primary = cf.em,
                               ploidy_primary = ploidy, purity_primary = purity, cn_state_primary = cn_state),
        ln_df %>% select(gene, chrom, status_ln = status, 
                          tcn_ln = tcn.em, cf_ln = cf.em,
                          ploidy_ln = ploidy, purity_ln = purity, cn_state_ln = cn_state),
        by = "gene"
    ) %>%
    mutate(
        status_primary = replace_na(status_primary, "Neutral"),
        status_ln      = replace_na(status_ln, "Neutral"),
        transition = case_when(
            status_primary == "Neutral" & status_ln == "Amp"  ~ "Acquired_Amp",
            status_primary == "Neutral" & status_ln == "Del"  ~ "Acquired_Del",
            status_primary == "Amp"     & status_ln == "Neutral" ~ "Lost_Amp",
            status_primary == "Del"     & status_ln == "Neutral" ~ "Lost_Del",
            status_primary == "Amp"     & status_ln == "Del"   ~ "Flip_AmpToDel",
            status_primary == "Del"     & status_ln == "Amp"   ~ "Flip_DelToAmp",
            status_primary == "Amp"     & status_ln == "Amp"   ~ "Persistent_Amp",
            status_primary == "Del"     & status_ln == "Del"   ~ "Persistent_Del",
            TRUE ~ "Neutral_both"
        ),
        patient = p_id
    )

    return(comparison)
}

make_gistic_seg <- function(ssegs, dipLogR, sample_id) {
  ssegs %>%
    transmute(
      Sample = sample_id,
      Chromosome = chrom,
      Start_Position = start,
      End_Position = end,
      Num_Markers = num.mark,
      Seg_CN = cnlr.median - dipLogR
    )
}


primary_gistic2_ls <- list()
ln_met_gistic2_ls <- list()

for(i in 1:nrow(clinical_data)) {
    row <- clinical_data[i, ]
    sample <-  row[["Sample"]]
    primary_tumor_id <- row[["DNA_primary_cegat"]]
    ln_tumor_id <- row[["DNA_ln_cegat"]]

    normal_id <- row[["DNA_normal_cegat"]]
    primary_rds_path <- paste(base_path, "primary/results", primary_tumor_id, 
                                paste0(primary_tumor_id, "_", normal_id, "_hisens.rds"), sep="/")


    facets_file <- readRDS(primary_rds_path)
    ssegs <-  facets_file$segs   
    ssegs <- ssegs %>% filter(nhet >= 10, !is.na(mafR))

    primary_gistic2_ls[[sample]] <-  make_gistic_seg(ssegs, facets_file$dipLogR, sample)

    ln_rds_path <- paste(base_path, "ln_met/results", ln_tumor_id, 
                       paste0(ln_tumor_id, "_", normal_id, "_hisens.rds"), sep="/")

    facets_file <- readRDS(ln_rds_path)
    ssegs <-  facets_file$segs   
    ssegs <- ssegs %>% filter(nhet >= 10, !is.na(mafR))

    ln_met_gistic2_ls[[sample]] <-  make_gistic_seg(ssegs, facets_file$dipLogR, sample)

}



gistic_df_primary <- do.call(rbind, primary_gistic2_ls)
write.table(gistic_df_primary, "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/GISTIC2/input/full_primary.seg", sep="\t", row.names=FALSE, quote=FALSE)

gistic_df_ln_met <- do.call(rbind, ln_met_gistic2_ls)
write.table(gistic_df_ln_met, "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/GISTIC2/input/full_ln_met.seg", sep="\t", row.names=FALSE, quote=FALSE)
