library(dplyr)
library(readxl)
library(readr)
library(tidyr)
library(survminer)
library(survival)
library(ggpubr)
library(DESeq2)

source("~/packages/R_Helpers/gene_name_mapper.R")
source("~/packages/R_Helpers/enrichment.R")

patient_annotation <- read.csv("/data_pool_1/NADINA/processed/patient_data_files/260914_summarised_nadina.csv")
patient_annotation$TMB_per_mb <- (as.numeric(patient_annotation$TMBmutmb) / 36715240) * 1e6



broad_data <- read_tsv("/data_pool_1/NADINA/processed/dna/oncoanalyser/neoadjuvant/output/GISTIC2/output_filtered/broad_values_by_arm.txt") 


amp_binarized <- broad_data %>%
  filter(`Chromosome Arm` == "1q") %>%
  mutate(across(-`Chromosome Arm`, ~ ifelse(. > 0, 1, 0)))

samples_1q_gain <- amp_binarized %>%
  pivot_longer(-`Chromosome Arm`, names_to = "Sample", values_to = "value") %>%
  filter(value == 1) %>%
  pull(Sample)


patient_annotation <- patient_annotation %>% 
    filter(ptid %in% colnames(broad_data)) %>%
    mutate(arm1q_gain = ptid %in% samples_1q_gain)

patient_annotation$arm1q_gain <- as.factor(patient_annotation$arm1q_gain)

#SANITY checks 
tmb_nadina <- ggplot(patient_annotation, aes(x = arm1q_gain, y = TMB_per_mb, fill = arm1q_gain)) +
 geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    linewidth = 0.5,
    color = "black",
    alpha = 0.75
  ) +
  geom_jitter(
    aes(fill = arm1q_gain),
    shape = 21,
    size = 1.8,
    alpha = 0.55,
    width = 0.15,
    show.legend = FALSE
  ) +   
  scale_fill_manual(values = c("#4A90D9", "#E05C5C")) +
  scale_x_discrete(labels = c("FALSE" = "1q neutal", "TRUE" = "1q gain")) +
  labs(
    x = "Chromosome 1q gain", 
    y = "Non-Synonymous TMB per MB"
  ) +
  geom_hline(
    yintercept = 10, 
    linetype = "dashed", 
    color = "gray30",    
    linewidth = 0.6      
  ) +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_blank(),
    axis.title.y     = element_text(size = 7, face = "bold", margin = margin(r = 8)),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(size = 9),
    axis.text.x      = element_text(size = 9),
    legend.position  = "none",
    legend.key.size  = unit(0.5, "cm"),
    panel.spacing    = unit(0, "lines"),
    plot.margin      = margin(10, 15, 10, 10)
  ) + stat_compare_means()

tmb_nadina




ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/CNA/arm1q_NADINA_tmb.pdf", 
  plot = tmb_nadina,
  width = 13.5,   
  height = 9.5,   
  units = "cm"
)




#Structural variants? 
base_path <- "/data_pool_1/NADINA/processed/dna/oncoanalyser/neoadjuvant/output/"

df_sv_tmb <- data.frame(
  ptid = patient_annotation$ptid,
  svTMB = NA_real_,
  ploidy = NA_real_,
  wgd = NA,
  arm1q_purple = NA_real_,
  arm1q_purple_call = NA,

  stringsAsFactors = FALSE
)

for (i in seq_along(patient_annotation$ptid)) {

  id <- patient_annotation$ptid[i]

  path <- file.path(base_path, id, "purple")
  file <- list.files(
    path,
    pattern = "\\.purple\\.purity\\.tsv$",
    full.names = TRUE
  )[1]

  if (is.na(file)) next

  purple_summary <- read_tsv(file, show_col_types = FALSE)
  ploidy <- purple_summary$ploidy[1]
  df_sv_tmb$svTMB[i] <- purple_summary$svTumorMutationalBurden[1]
  df_sv_tmb$ploidy[i] <- ploidy
  df_sv_tmb$wgd[i] <- purple_summary$wholeGenomeDuplication[1]

}


patient_annotation <- left_join(patient_annotation, df_sv_tmb)

#structural variants

svtmb <- ggplot(patient_annotation, aes(x = arm1q_gain, y = svTMB, fill = arm1q_gain)) +
 geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    linewidth = 0.5,
    color = "black",
    alpha = 0.75
  ) +
  geom_jitter(
    aes(fill = arm1q_gain),
    shape = 21,
    size = 1.8,
    alpha = 0.55,
    width = 0.15,
    show.legend = FALSE
  ) +   
  scale_fill_manual(values = c("#4A90D9", "#E05C5C")) +
  scale_x_discrete(labels = c("FALSE" = "1q neutal", "TRUE" = "1q gain")) +
  labs(
    x = "Chromosome 1q gain", 
    y = "Continous Ploidy estimates"
  ) +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_blank(),
    axis.title.y     = element_text(size = 7, face = "bold", margin = margin(r = 8)),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(size = 9),
    axis.text.x      = element_text(size = 9),
    legend.position  = "none",
    legend.key.size  = unit(0.5, "cm"),
    panel.spacing    = unit(0, "lines"),
    plot.margin      = margin(10, 15, 10, 10)
  ) + stat_compare_means()

svtmb

ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/CNA/arm1q_NADINA_ploidy.pdf", 
  plot = svtmb,
  width = 13.5,   
  height = 9.5,   
  units = "cm"
)




#Not more Structural variants!


surv_obj <- Surv(time = patient_annotation$Time_EFS, event = patient_annotation$EFS_event)
fit <- survfit(surv_obj ~ arm1q_gain, data = patient_annotation)



gain_vs_wt <- ggsurvplot(
  fit,
  data = patient_annotation,
  pval = TRUE,
  pval.size = 6,                 # Increase p-value font size
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.fontsize = 8,       # Risk table text size
  legend.title = "",
  # legend.labs = c("1q neutral | IFNy low", "1q neutral | IFNy High", "1q gain | IFNy low", "1q gain | IFNy high"), 
  xlab = "Time (Months)",
  ylab = "Event-Free Survival Probability",
  
  # palette = c("#5d666b", "#D55E00", "#0072B2", "#009E73")  ,
  
  # Increase overall font sizes
  font.main = c(20, "bold"),
  font.x = c(32),
  font.y = c(32),
  font.tickslab = c(28),
  font.legend = c(20),
  
  ggtheme = theme_classic(base_size = 16)
)
gain_vs_wt


pdf("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/CNA/NADINA_EFSS_1q.pdf", 
    width = 15,
    height = 10)

gain_vs_wt 
dev.off()



ifn_filtered_df <- patient_annotation %>% filter(IFNcat != ".") #%>% filter(IFNcat == "IFNy_low") #%>% filter(IFNcat != ".")

fisher.test(table(ifn_filtered_df$arm1q_gain, ifn_filtered_df$IFNcat))
cor.test(as.numeric(ifn_filtered_df$arm1q_gain), as.numeric(ifn_filtered_df$IFNscore), method = "spearman")

surv_obj <- Surv(time = ifn_filtered_df$Time_EFS, event = ifn_filtered_df$EFS_event)
fit <- survfit(surv_obj ~ arm1q_gain + IFNcat, data = ifn_filtered_df)


chr1qSurvPlot <- ggsurvplot(
  fit, 
  pval = TRUE, 
  risk.table = TRUE,
  risk.table.col = "strata", 
  pval.coord  = c(0.3,0.15),
  legend = "none",      
  conf.int = FALSE,
  legend.labs = c("1q neutral | IFNy high", "1q neutral | IFNy low", "1q gain | IFNy high", "1q gain | IFNy low"), 
  size = 1.5, censor.shape = 124,
  palette = c("#D55E00", "#5d666b",  "#009E73", "#0072B2"),
  xlim = c(0, 40), break.x.by = 12, surv.scale = "percent", break.y.by = 0.2
)
chr1qSurvPlot
surv_summarypat_36 <- summary(fit, times = 35)
surv_36allo <- list()

for (stratum in 1:4) {
  surv_36allo[[stratum]] <- surv_summarypat_36$surv[stratum]
}

for (i in 1:4) {
  chr1qSurvPlot$plot <- chr1qSurvPlot$plot + 
    geom_segment(aes(x = 35, y = 0, xend = 35, yend = surv_36allo[[i]]), 
                 linetype = "dashed", color = "black") + 
    annotate("text", x = 35, y = surv_36allo[[i]], 
             label = paste0(round(surv_36allo[[i]] * 100, 0), "%"), 
             vjust = -1, hjust = -0.1, color = "black")
}
chr1qSurvPlot$plot <- chr1qSurvPlot$plot + labs(x = "Months since registration", y = "Event-free survival (%)")
print(chr1qSurvPlot)






pdf(
  "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/CNA/NADINA_EFSS_1q.pdf",
  width = 9.13,
  height = 6.22
)

chr1qSurvPlot

dev.off()



testing_efs_low_ifn <- patient_annotation %>% filter(IFNcat == "IFNy_low")
surv_obj <- Surv(time = testing_efs_low_ifn$Time_EFS, event = testing_efs_low_ifn$EFS_event)
fit <- survfit(surv_obj ~ arm1q_gain, data = testing_efs_low_ifn)



low_gain_vs_wt <- ggsurvplot(
  fit,
  data = testing_efs_low_ifn,
  pval = TRUE,
  pval.size = 6,                 # Increase p-value font size
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.fontsize = 8,       # Risk table text size
  legend.title = "",
  # legend.labs = c("1q neutral | IFNy low", "1q neutral | IFNy High", "1q gain | IFNy low", "1q gain | IFNy high"), 
  xlab = "Time (Months)",
  ylab = "Event-Free Survival Probability",
  
  # palette = c("#5d666b", "#D55E00", "#0072B2", "#009E73")  ,
  
  # Increase overall font sizes
  font.main = c(20, "bold"),
  font.x = c(32),
  font.y = c(32),
  font.tickslab = c(28),
  font.legend = c(20),
  
  ggtheme = theme_classic(base_size = 16)
)
low_gain_vs_wt
