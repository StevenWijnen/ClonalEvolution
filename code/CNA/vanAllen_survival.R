library(dplyr)
library(readxl)
library(readr)
library(tidyr)
library(survminer)
library(survival)
library(ggpubr)



patient_annotation <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/GISTIC2/vanAllen/skcm_dfci_2015_clinical_data.tsv", sep="\t")

patient_annotation <- patient_annotation %>% filter(Primary.Tumor.Site != "mucosal")
arm_changes <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/GISTIC2/vanAllen/output/broad_values_by_arm.txt", sep="\t")
arm1q_patients <- colnames(arm_changes)[-1][arm_changes[2, -1] > 0]

all(patient_annotation$Patient.ID %in% colnames(arm_changes)[-1])

df <- patient_annotation %>% mutate(arm1q_gain = ifelse(Patient.ID %in% arm1q_patients,TRUE, FALSE))


#patient_annotation$Overall.Survival..Months.
df$Overall.Survival.Status <- as.numeric(sapply(strsplit(df$Overall.Survival.Status, ":"), `[[`, 1))
df$Disease.Free.status <- as.numeric(sapply(strsplit(df$Disease.Free.status, ":"), `[[`, 1))


cutoff <- 36
df$time36   <- pmin(df$Overall.Survival..Months., cutoff)
df$status36 <- ifelse(df$Overall.Survival..Months. > cutoff, 0, df$Overall.Survival.Status)


surv_obj <- Surv(time = df$time36, event = df$status36)
fit <- survfit(Surv(time36, status36) ~ arm1q_gain, data = df)



p <- ggsurvplot(
  fit,
  data = df,
  pval = TRUE,
  pval.size = 6,
  pval.coord = c(1, 0.15),          # p-value position
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.fontsize = 8,
  legend.title = "",
  xlab = "Time (Months)",
  ylab = "Overall Survival Probability",
  xlim = c(0, cutoff),
  break.time.by = 6,
  font.main = c(20, "bold"),
  font.x = c(32),
  font.y = c(32),
  font.tickslab = c(28),
  font.legend = c(20),
  ggtheme = theme_classic(base_size = 16)
)





surv_obj <- Surv(time = df$time36, event = df$status36)
fit <- survfit(surv_obj ~ arm1q_gain, data = df)

cox <- coxph(Surv(time36, status36) ~ arm1q_gain, data = df)
cox_sum <- summary(cox)
hr    <- cox_sum$conf.int[1, "exp(coef)"]
lower <- cox_sum$conf.int[1, "lower .95"]
upper <- cox_sum$conf.int[1, "upper .95"]

hr_label <- sprintf("HR = %.2f (95%% CI %.2f–%.2f)", hr, lower, upper)

chr1qSurvPlot <- ggsurvplot(
  fit, 
  pval = TRUE, 
  risk.table = TRUE,
  risk.table.col = "strata", 
  pval.coord  = c(0.3,0.15),
  legend = "none",      
  conf.int = FALSE,
  legend.labs = c("1q neutral", "1q gain"), 
  size = 1.5, censor.shape = 124,
  palette = c("#5d666b", "#0072B2"),
  xlim = c(0, 36), break.x.by = 12, surv.scale = "percent", break.y.by = 0.2
)
chr1qSurvPlot
surv_summarypat_36 <- summary(fit, times = 50)
surv_36allo <- list()

for (stratum in 1:2) {
  surv_36allo[[stratum]] <- surv_summarypat_36$surv[stratum]
}


chr1qSurvPlot$plot <- chr1qSurvPlot$plot +
  annotate("text", x = 0, y = 0.05, label = hr_label,
           hjust = 0, size = 6)

chr1qSurvPlot$plot <- chr1qSurvPlot$plot + labs(x = "Time (Months)", y = "Overall survival (%)")
print(chr1qSurvPlot)



pdf(
  "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/CNA/VanAllenchr1q.pdf",
  width = 9.13,
  height = 6.22
)

chr1qSurvPlot

dev.off()







tmb_vanAllen<-ggplot(df, aes(x = arm1q_gain, y = TMB..nonsynonymous., fill = arm1q_gain)) +
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
  ) +   scale_fill_manual(values = c("#4A90D9", "#E05C5C")) +
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

tmb_vanAllen


ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/CNA/arm1q_vanAllen_tmb.pdf", 
  plot = tmb_vanAllen,
  width = 13.5,    # Adjusted for a quarter page layout with margins
  height = 9.5,    # Keeps a clean, slightly wider landscape aspect ratio
  units = "cm"
)









#RNA STUFF
patient_annotation_rna_filtered <- df %>% filter(!is.na(CBSET.T.cells.CD8))
median_cd8 <- median(patient_annotation_rna_filtered$CBSET.T.cells.CD8)

patient_annotation_rna_filtered$high_cd8 <- patient_annotation_rna_filtered$CBSET.T.cells.CD8 > median_cd8


# cd8_cutoffs <- quantile(
#   patient_annotation_rna_filtered$CBSET.T.cells.CD8,
#   probs = c(1/3, 2/3),
#   na.rm = TRUE
# )

# patient_annotation_rna_filtered$high_cd8 <- patient_annotation_rna_filtered$CBSET.T.cells.CD8 > 0.07677361

# Assign groups
# patient_annotation_rna_filtered$cd8_group <- cut(
#   patient_annotation_rna_filtered$CBSET.T.cells.CD8,
#   breaks = c(-Inf, cd8_cutoffs[1], cd8_cutoffs[2], Inf),
#   labels = c("low", "medium", "high"),
#   include.lowest = TRUE
# )

# surv_obj <- Surv(time = patient_annotation_rna_filtered$Overall.Survival..Months., event = patient_annotation_rna_filtered$Overall.Survival.Status)
surv_obj <- Surv(time = patient_annotation_rna_filtered$Disease.Free.Survival..Months., event = patient_annotation_rna_filtered$Disease.Free.status)
fit <- survfit(surv_obj ~ high_cd8  , data = patient_annotation_rna_filtered)


ggsurvplot(
  fit,
  data = patient_annotation_rna_filtered,
  pval = TRUE,
  pval.size = 6,                 # Increase p-value font size
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.fontsize = 8,       # Risk table text size
  legend.title = "",

  xlab = "Time (Months)",
  ylab = "Event-Free Survival Probability",
  
  
  # Increase overall font sizes
  font.main = c(20, "bold"),
  font.x = c(32),
  font.y = c(32),
  font.tickslab = c(28),
  font.legend = c(20),
  
  ggtheme = theme_classic(base_size = 16)
)

