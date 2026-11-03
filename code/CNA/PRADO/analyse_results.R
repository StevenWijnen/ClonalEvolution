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

patient_annotation <- read_xlsx("/DATA/s.wijnen/Projects/PRADO_genomic_analyses/DATA/20250612_IFNstageIII - database - PRADO - 5year - Steven_clean.xlsx") %>% 
    dplyr::select(ptid = `Patient study code`, BRAF, MPR1, MPR2, Event, Time_EFS)
patient_ids <- read_xlsx("/DATA/s.wijnen/Projects/PRADO_genomic_analyses/DATA/PRADO_annotation_DNA_RNA_complete_20221124.xlsx") %>% 
    dplyr::select(ptid=studyID, CagetID_DNA_tumor, CagetID_RNA, IFNy_HL)
tmb <- read.csv("/DATA/s.wijnen/Projects/PRADO_genomic_analyses/output/FINAL_tmb_from_vep_vcf.csv")

patient_data <- left_join(patient_annotation, patient_ids)
patient_data <- left_join(patient_data, tmb, by=c("CagetID_DNA_tumor" ="ID"))

broad_data <- read_tsv("/data_pool_1/PRADO/processed/dna/SAREK_pipeline/output/variant_calling/GISTIC2/output/broad_values_by_arm.txt")

# Get 1p and 1q side by side per sample
arm_values <- broad_data %>%
  filter(`Chromosome Arm` %in% c("1p", "1q")) %>%
  pivot_longer(-`Chromosome Arm`, names_to = "Sample", values_to = "value") %>%
  pivot_wider(names_from = `Chromosome Arm`, values_from = value)

# Two binarized criteria
amp_binarized <- arm_values %>%
  mutate(
    q_gt_0 = ifelse(`1q` > 0, 1, 0),
    q_gt_p = ifelse(`1q` > `1p`, 1, 0)
  )



amp_binarized <- left_join(amp_binarized, patient_data, by=c("Sample" = "CagetID_DNA_tumor"))
amp_binarized$q_gt_0 <- as.factor(amp_binarized$q_gt_0)

tmb_prado<-ggplot(amp_binarized, aes(x = q_gt_0, y = nonSyn_TMB_per_mb, fill = q_gt_0)) +
 geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    linewidth = 0.5,
    color = "black",
    alpha = 0.75
  ) +
  geom_jitter(
    aes(fill = q_gt_0),
    shape = 21,
    size = 1.8,
    alpha = 0.55,
    width = 0.15,
    show.legend = FALSE
  ) +   scale_fill_manual(values = c("#4A90D9", "#E05C5C")) +
  scale_x_discrete(labels = c("0" = "1q neutal", "1" = "1q gain")) +
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

tmb_prado



# ggsave(
#   filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/CNA/arm1q_prado_tmb.pdf", 
#   plot = tmb_prado,
#   width = 13.5,    # Adjusted for a quarter page layout with margins
#   height = 9.5,    # Keeps a clean, slightly wider landscape aspect ratio
#   units = "cm"
# )

surv_obj <- Surv(time = amp_binarized$Time_EFS, event = amp_binarized$Event)
fit <- survfit(surv_obj ~ q_gt_0 + IFNy_HL, data = amp_binarized)

coxph(Surv(Time_EFS, Event) ~ q_gt_0 + IFNy_HL, data = amp_binarized)


chr1qSurvPlot <- ggsurvplot(
  fit, 
  pval = TRUE, 
  risk.table = TRUE,
  risk.table.col = "strata", 
  pval.coord  = c(0.3,0.15),
  legend = "none",      
  conf.int = FALSE,
  legend.labs = c("1q neutral | IFNy low", "1q neutral | IFNy High", "1q gain | IFNy low", "1q gain | IFNy high"), 
  size = 1.5, censor.shape = 124,
  palette = c("#5d666b", "#D55E00", "#0072B2", "#009E73")  ,
  xlim = c(0, 72), break.x.by = 12, surv.scale = "percent", break.y.by = 0.2
)

surv_summarypat_36 <- summary(fit, times = 60)
surv_36allo <- list()

for (stratum in 1:4) {
  surv_36allo[[stratum]] <- surv_summarypat_36$surv[stratum]
}

for (i in 1:4) {
  chr1qSurvPlot$plot <- chr1qSurvPlot$plot + 
    geom_segment(aes(x = 60, y = 0, xend = 60, yend = surv_36allo[[i]]), 
                 linetype = "dashed", color = "black") + 
    annotate("text", x = 60, y = surv_36allo[[i]], 
             label = paste0(round(surv_36allo[[i]] * 100, 0), "%"), 
             vjust = -1, hjust = -0.1, color = "black")
}
chr1qSurvPlot$plot <- chr1qSurvPlot$plot + labs(x = "Months since registration", y = "Event-free survival (%)")
print(chr1qSurvPlot)


testing_ifn_low <- amp_binarized %>% filter(IFNy_HL == 0)
surv_obj <- Surv(time = testing_ifn_low$Time_EFS, event = testing_ifn_low$Event)
fit <- survfit(surv_obj ~ q_gt_0 + IFNy_HL, data = testing_ifn_low)

coxph(Surv(Time_EFS, Event) ~ q_gt_0 + IFNy_HL, data = testing_ifn_low)


chr1qSurvPlot <- ggsurvplot(
  fit, 
  pval = TRUE, 
  risk.table = TRUE,
  risk.table.col = "strata", 
  pval.coord  = c(0.3,0.15),
  legend = "none",      
  conf.int = FALSE,
  # legend.labs = c("1q neutral | IFNy low", "1q neutral | IFNy High", "1q gain | IFNy low", "1q gain | IFNy high"), 
  size = 1.5, censor.shape = 124,
  # palette = c("#5d666b", "#D55E00", "#0072B2", "#009E73")  ,
  xlim = c(0, 72), break.x.by = 12, surv.scale = "percent", break.y.by = 0.2
)


chr1qSurvPlot


# pdf(
#   "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/CNA/PRADO_EFSS_1q.pdf",
#   width = 9.13,
#   height = 6.22
# )

# chr1qSurvPlot

# dev.off()







count_matrix <- read.csv("/data_pool_1/PRADO/processed/rna/Count_matrix_HTseq_PRADO.txt", sep="\t")
dim(count_matrix)
count_matrix <- count_matrix[startsWith(rownames(count_matrix), "ENSG") ,]
dim(count_matrix)

rna_filtered_amp_binarized <- amp_binarized %>% filter(CagetID_RNA!= "NA")
count_matrix <- count_matrix[,rna_filtered_amp_binarized$CagetID_RNA]
rna_filtered_amp_binarized$q_gt_0 <- as.factor(rna_filtered_amp_binarized$q_gt_0)

rna_filtered_amp_binarized$IFNy_HL <- as.factor(rna_filtered_amp_binarized$IFNy_HL)

dds <- DESeqDataSetFromMatrix(colData=rna_filtered_amp_binarized, 
                              countData=count_matrix, 
                              design = ~ q_gt_0 + IFNy_HL
)

keep_counts <- rowSums(counts(dds)) >= 10
dds <- dds[keep_counts, ]

keep_samples <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep_samples, ]

dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

plotPCA(vsd, intgroup = "q_gt_0", returnData = FALSE)



res_mut <- (results(dds, contrast = c("q_gt_0", "1", "0")))
res_mut_df <- as.data.frame(res_mut)

named_res <- ensembl_to_hgnc(res_mut_df)


write.csv(as.data.frame(named_res) %>% rownames_to_column(var="Gene"), "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/arm1qgain/prado_deg_winfy.csv", quote=FALSE, row.names=FALSE)

