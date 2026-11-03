library(dplyr)
library(purrr)
library(ggplot2)
library(tibble)
library(tidyr)

clinical_data <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/filtered_cleaned_df.csv")
clinical_data$IFN_ln <- ifelse(clinical_data$IFN %in% c("primhigh_LNlow", "primlow_LNlow"), "low", "high")
broad_significant_hits <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/GISTIC2/output/ln_met_full/broad_significance_results.txt", sep="\t")

significant_amps <- broad_significant_hits %>% 
  filter(Amp.q.value < 0.1) %>% 
  pull(Arm)

significant_dels <- broad_significant_hits %>% 
  filter(Del.q.value < 0.1) %>% 
  pull(Arm) 


per_sample_broad_value <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/GISTIC2/output/ln_met_full/broad_values_by_arm.txt", sep="\t")


per_sample_broad_value <- read.csv(
  "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/GISTIC2/output/ln_met_full/broad_values_by_arm.txt",
  sep = "\t"
)

amp_binarized <- per_sample_broad_value %>%
  filter(Chromosome.Arm %in% significant_amps) %>%
  mutate(across(-Chromosome.Arm, ~ ifelse(. > 0, 1, 0)))

del_binarized <- per_sample_broad_value %>%
  filter(Chromosome.Arm %in% significant_dels) %>%
  mutate(across(-Chromosome.Arm, ~ ifelse(. < 0, -1, 0)))


amp_binarized_t <- amp_binarized %>%
  pivot_longer(-Chromosome.Arm, names_to = "Sample", values_to = "value") %>%
  pivot_wider(names_from = Chromosome.Arm, values_from = value)
#Entire arm not specific to 

amp_binarized_t <- amp_binarized_t %>%
  mutate(`1q` = ifelse(Sample == "PRADO_193", 0, `1q`))


del_binarized_t <- del_binarized %>%
  pivot_longer(-Chromosome.Arm, names_to = "Sample", values_to = "value") %>%
  pivot_wider(names_from = Chromosome.Arm, values_from = value)




amp_df <- amp_binarized_t %>% left_join(
    clinical_data %>% dplyr::select(Sample, IFN), by="Sample"
)
# amp_df[amp_df$Sample == "OpACINneo_105", "IFN"] <- "primlow_LNhigh"

arm_vars <- setdiff(colnames(amp_df), c("Sample", "IFN"))


amp_fisher_results <- map_dfr(arm_vars, function(arm) {
  
  tab <- table(amp_df[[arm]], amp_df$IFN)
  
  # Fisher's test; simulate p-value if the table is too large for exact computation
  test <- tryCatch(
    fisher.test(tab),
    error = function(e) fisher.test(tab)
  )
  
  tibble(
    Arm = arm,
    p.value = test$p.value,
    table = list(tab)
  )
}) 

# View sorted results
amp_fisher_results %>% 
  dplyr::select(Arm, p.value) %>% 
  arrange(p.value)

amp_fisher_results <- amp_fisher_results %>%
  mutate(p.adjust = p.adjust(p.value, method = "hommel"))

amp_fisher_results

del_df <- del_binarized_t %>% 
  left_join(clinical_data %>% dplyr::select(Sample, IFN), by = "Sample")

arm_vars_del <- setdiff(colnames(del_df), c("Sample", "IFN"))

del_fisher_results <- map_dfr(arm_vars_del, function(arm) {
  
  tab <- table(del_df[[arm]], del_df$IFN)
  
  test <- tryCatch(
    fisher.test(tab),
    error = function(e) fisher.test(tab, simulate.p.value = TRUE, B = 1e5)
  )
  
  tibble(
    Arm = arm,
    p.value = test$p.value,
    table = list(tab)
  )
})

del_fisher_results %>% 
  dplyr::select(Arm, p.value) %>% 
  arrange(p.value)

test <- rbind(del_fisher_results, amp_fisher_results) %>%
  mutate(p.adjust = p.adjust(p.value, method = "hommel"))






data_long_1q <- amp_df %>%
  dplyr::select(IFN, `1q`) %>%
  count(IFN, `1q`) %>%
  group_by(IFN) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(`1q` = factor(`1q`, levels = c(0, 1)))

# n= labels for x-axis
n_labels_1q <- amp_df %>%
  count(IFN) %>%
  mutate(label = paste0(IFN, "\n(n=", n, ")")) %>%
  { setNames(.$label, .$IFN) }

arm_1q_plot <- ggplot(data_long_1q, aes(x = IFN, y = percentage, fill = `1q`)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 2, color = "white", fontface = "bold") +
  scale_fill_manual(
    values = c("0" = "#4A90D9", "1" = "#E05C5C"),
    labels = c("0" = "No 1q gain", "1" = "1q gain")
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_discrete(labels = n_labels_1q) +
  labs(
    title = "",
    x     = "IFN Group",
    y     = "Percentage",
    fill  = "1q status"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold"),
    legend.position    = "right",
    panel.grid.major.x = element_blank()
  )

arm_1q_plot

ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/arm1q_plot.pdf", 
  plot = arm_1q_plot,
  width = 13.5,    # Adjusted for a quarter page layout with margins
  height = 9.5,    # Keeps a clean, slightly wider landscape aspect ratio
  units = "cm"
)



#PEAK significance
sign_peaks  <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/GISTIC2/output/ln_met_full/all_lesions.conf_75.txt", sep="\t")[1:7,]
peaks_sub <- sign_peaks[, c(2, 11:(ncol(sign_peaks) - 1))]

peaks_sub$Descriptor <- make.unique(as.character(peaks_sub$Descriptor))

mat <- as.data.frame(t(peaks_sub[, -1]))
colnames(mat) <- peaks_sub$Descriptor
mat <- rownames_to_column(mat, "Sample")

mat <- mat %>% mutate(across(-Sample, ~ as.numeric(sign(.x))))

test_df <- mat %>% left_join(
    clinical_data %>% dplyr::select(Sample, IFN, IFN_ln), by="Sample"
)

fisher.test(table(test_df$`1p21.3`, test_df$IFN))
table(test_df$`1p21.3`, test_df$IFN)
