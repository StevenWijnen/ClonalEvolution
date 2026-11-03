library(GenVisR)
library(readr)
library(dplyr)
library(ggplot2)


patient_data <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/filtered_cleaned_df.csv")
good_samples <- c("DONIMI_D14","DONIMI_D24","DONIMI_D4","OpACINneo_105","OpACINneo_16","OpACINneo_21","OpACINneo_28","OpACINneo_41","OpACINneo_54","OpACINneo_9","OpACINneo_98","PRADO_107","PRADO_111","PRADO_112","PRADO_117","PRADO_119","PRADO_120","PRADO_121","PRADO_134","PRADO_137","PRADO_152","PRADO_157","PRADO_164","PRADO_182","PRADO_186","PRADO_193")
base_path <- "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/facets_redux"
patient_data <- patient_data %>% filter(Sample %in% good_samples)

process_facets_ <- function(file_path, origin) {
        facets_file <- readRDS(file_path)
        sample_df <- facets_file$segs
        
        sample_df <- sample_df %>% filter(nhet >= 10)
        sample_ploidy <- facets_file$ploidy
        

        sample_df <- sample_df %>% dplyr::select(chrom, start, end, tcn.em, lcn.em, cnlr.median.clust, cnlr.median)
        sample_df$relative_change <- sample_df$tcn.em - sample_ploidy

        sample_df$gain_or_loss <- sign(sample_df$relative_change)
        sample_df$length <- sample_df$end - sample_df$start

        # sample_df <- sample_df %>% filter(length >= 50000)

    return(sample_df)
}

process_patient <- function(p_id, normal_id, primary_tumor_id, ln_tumor_id) {

    primary_rds_path <- paste(base_path, "primary/results", primary_tumor_id, 
                            paste0(primary_tumor_id, "_", normal_id, "_hisens.rds"), sep="/")
  
    ln_rds_path <- paste(base_path, "ln_met/results", ln_tumor_id, 
                       paste0(ln_tumor_id, "_", normal_id, "_hisens.rds"), sep="/")

    primary_df <- process_facets_(primary_rds_path, "primary")
    primary_df$patient_id <- p_id
    primary_df$sample <- primary_tumor_id
    primary_df$origin <- "primary"

    ln_df <- process_facets_(ln_rds_path, "ln_met")
    ln_df$patient_id <- p_id
    ln_df$sample <- ln_tumor_id
    ln_df$origin <- "ln_met"
    
    
    return(rbind(primary_df, ln_df))
}


results_list <- list()
for(i in 1:nrow(patient_data)) {
  row <- patient_data[i, ]
  case_id    <- row[["sample"]]
  normal_id  <- row[["DNA_normal_cegat"]]
  primary_id <- row[["DNA_primary_cegat"]]
  ln_tumor_id <- row[["DNA_ln_cegat"]]

    results_list[[i]] <- process_patient(case_id, normal_id, primary_id, ln_tumor_id)
    results_list[[i]]$group <- row[["IFN"]]
    results_list[[i]]$response <- row[["MPR"]]

}


final_master_df <- do.call(rbind, results_list)
chr_list <- paste0("chr", seq(1:22))

GenVis_df_primary <- final_master_df %>% filter(origin == "primary") %>% dplyr::select(chromosome = chrom, start, end, segmean = relative_change, sample = sample)
GenVis_df_ln_met <- final_master_df %>% filter(origin == "ln_met") %>% dplyr::select(chromosome = chrom, start, end, segmean = relative_change, sample = sample)


pub_theme <- theme(
  plot.title    = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10)),
  axis.title.x  = element_text(size = 8, face = "bold", margin = margin(t = 10)),
  axis.title.y  = element_text(size = 10, face = "bold", margin = margin(r = 10)),
  axis.text.y   = element_text(size = 8),
  legend.title  = element_text(size = 13, face = "bold"),
  legend.text   = element_text(size = 11),
  legend.position = "bottom",

  panel.spacing = unit(0.5, "lines"),
  plot.margin   = margin(t = 15, r = 20, b = 25, l = 20)  # slightly increased b
)


plot_primary <- cnFreq(
  GenVis_df_primary,
  genome        = "hg38",
  CN_low_cutoff  = -0.5,
  CN_high_cutoff =  0.5,
  plot_title    = "",
  plotChr       = chr_list,
  plotType = "proportion"
)  + pub_theme +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  theme(strip.text.x = element_text(size = 6, angle = 0))

plot_primary


pdf("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/primary_cna.pdf",
      width  = 7.25, 
        height = 3.6,  
    useDingbats = FALSE)
plot_primary
dev.off()



plot_ln_met <- cnFreq(
  GenVis_df_ln_met,
  genome        = "hg38",
  CN_low_cutoff  = -0.5,
  CN_high_cutoff =  0.5,
  plot_title    = "",
  plotChr       = chr_list,
  plotType = "proportion"

) + pub_theme +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +
  theme(strip.text.x = element_text(size = 6, angle = 0))

plot_ln_met


pdf("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/ln_met_cna.pdf",
      width  = 7.25, 
        height = 3.6,  
    useDingbats = FALSE)
plot_ln_met
dev.off()

plot_group <- function(df, target_group = "primlow_LNhigh") {
    #Low to high primary vs lnmet
    GenVis_df_primary_ <- df %>% 
        filter(origin == "primary", group == target_group) %>% 
        dplyr::select(chromosome = chrom, start, end, segmean = relative_change, sample = sample) #4 samples

    GenVis_df_ln_met_ <- df %>% filter(origin == "ln_met", group == target_group) %>% dplyr::select(chromosome = chrom, start, end, segmean = relative_change, sample = sample)
    n <- length(unique(GenVis_df_ln_met_$sample))
    plot_primary_ <- cnFreq(
        GenVis_df_primary_, 
        genome="hg38", 
        CN_low_cutoff = -0.5, CN_high_cutoff = 0.5, 
        plot_title= paste("primary population", target_group, "(n = ", n, ")"),
        plotChr =chr_list
        
    )  + pub_theme +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%"))+
    theme(strip.text.x = element_text(size = 6, angle = 0))


    plot_ln_met_ <- cnFreq(
            GenVis_df_ln_met_, 
            genome="hg38", 
            CN_low_cutoff = -0.5, CN_high_cutoff = 0.5, 
            plotChr =chr_list, 
            plot_title    = paste("lymph node population", target_group, "(n = ", n, ")"),
            plotType = "proportion"
        ) + pub_theme +
    scale_y_continuous(labels = function(x) paste0(x * 100, "%"))+
    theme(strip.text.x = element_text(size = 6, angle = 0))

    pdf(paste0("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/prim_", target_group, ".pdf"), 
        width  = 7.25, 
        height = 3.6,  
        useDingbats = FALSE)

        print(plot_primary_   )
    dev.off()

    pdf(paste0("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/ln_met_", target_group, ".pdf"), 
      width  = 7.25, 
        height = 3.6,  
        useDingbats = FALSE)

        print(plot_ln_met_   )
    dev.off()
}


plot_group(final_master_df, "primlow_LNlow")
plot_group(final_master_df, "primhigh_LNhigh")

plot_group(final_master_df, "primlow_LNhigh")
plot_group(final_master_df, "primhigh_LNlow")
