setwd("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/code/nanostring_analyses")
library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(ggpubr)
source("./load_data.R")

condition_colors <- c(
  "primhigh_LNhigh" = "#8b1812",
  "primhigh_LNlow"  = "#79a8da",
  "primlow_LNhigh"  = "#dea562",
  "primlow_LNlow"   = "#35398e"
)

IFN_labels <- list(
  "primhigh_LNhigh" = "Primary High / LN High",
  "primhigh_LNlow"  = "Primary High / LN Low",
  "primlow_LNhigh"  = "Primary Low / LN High",
  "primlow_LNlow"   = "Primary Low / LN Low"
)
condition_labels_x <- c(
  "prim_HLA_combined" = "Primary",
  "ln_met_HLA_combined" = "Lymph node"
)




prim_expression_data <- nCounter_expression_primary %>% 
  mutate(prim_HLA_combined = rowMeans(dplyr::select(., all_of(c("HLA.A","HLA.B","HLA.C"))), na.rm = TRUE)) %>%
  dplyr::select(prim_HLAA = HLA.A, prim_HLAB = HLA.B, prim_HLAC = HLA.C, prim_HLA_combined, Trial, StudyID)

ln_met_expression_data <- nCounter_expression_ln_met %>% 
  mutate(ln_met_HLA_combined = rowMeans(dplyr::select(., all_of(c("HLA.A","HLA.B","HLA.C"))), na.rm = TRUE)) %>%
  dplyr::select(ln_met_HLAA = HLA.A,ln_met_HLAB = HLA.B,ln_met_HLAC = HLA.C, ln_met_HLA_combined, Trial, StudyID)


full_nanostring_expression_mat <- inner_join(prim_expression_data, ln_met_expression_data)
full_nanostring_df <- inner_join(full_nanostring_df, full_nanostring_expression_mat)


full_df_long <- full_nanostring_df %>%
  pivot_longer(
    cols = !c(Sample, annonID, Trial, StudyID, IFN, IFNln),
    names_to = "Condition",
    values_to = "Expression"
  ) %>% mutate(Condition = factor(Condition, levels=c(
    "prim_HLAA", "ln_met_HLAA",
    "prim_HLAB", "ln_met_HLAB",
    "prim_HLAC", "ln_met_HLAC",
    "prim_HLA_combined", "ln_met_HLA_combined"
    )))

my_jitter <- position_jitterdodge(
  jitter.width = 0,    # Use a dot
  dodge.width = 0.75,   # Use a dot
  seed = 42
)


full_df_long %>% 
  filter(Condition %in% c("prim_HLA_combined", "ln_met_HLA_combined")) %>%
  # Put Condition on X to spread them out, and fill by Condition
  ggplot(aes(x = Condition, y = Expression, fill = Condition)) + 
  
  # Boxplot with a controlled width
  geom_boxplot(
    outlier.shape = NA,
    color = "black",
    alpha = 0.75,
    width = 0.6  # Shrinks/widens the box itself to look balanced
  ) +
  
  # Connected lines for paired samples (if applicable)
  geom_line(aes(group = Sample), color = "grey40", alpha = 0.5) +
  
  # Single, clean jitter layer
  geom_jitter(
    color = "black",
    shape = 21,
    size = 1.8,
    alpha = 0.55,
    position = my_jitter,
    show.legend = FALSE
  ) +
  
  # Repel text that actually aligns with the data points
  geom_text_repel(
    aes(label = StudyID),
    size = 3,
    segment.color = 'grey50'
  ) +
  
  # Facet by IFN ONLY if you have multiple different genes/groups to compare side-by-side
  facet_wrap(~ IFN, nrow = 1, scales = "free_x") + 
  
  scale_fill_discrete(name = "Site") +
  theme_minimal() +
  theme(
    # Panel spacing to eliminate massive white gaps
    panel.spacing   = unit(1, "lines"), 
    
    # Clean up facet strip appearance
    strip.background = element_blank(),
    strip.text       = element_text(size = 14, face = "bold"), 
    
    # Text scaling
    axis.title = element_text(size = 18, face = "bold"), 
    axis.text  = element_text(size = 14),                
    legend.title = element_text(size = 16),              
    legend.text = element_text(size = 14)               
  ) + 
  labs(
    x = "Condition",
    y = "Expression Level"
  )





pval_df <- lapply(names(IFN_labels), function(n) {
  sub_df <- full_nanostring_df %>%
    filter(IFN == n) %>%
    select(Sample, prim_HLA_combined, ln_met_HLA_combined) %>%
    arrange(Sample)  # ensure pairs are aligned

  test <- t.test(
    sub_df$prim_HLA_combined,
    sub_df$ln_met_HLA_combined,
    paired = TRUE
  )

  data.frame(
    IFN    = n,
    pval   = test$p.value,
    label  = ifelse(test$p.value < 0.001, "p<0.001",
             ifelse(test$p.value < 0.01,  "p<0.01",
             ifelse(test$p.value <= 0.05,  "p<0.05",
             ifelse(test$p.value < 0.049,  paste0("p=", round(test$p.value, 3)),
                                          paste0("p=", round(test$p.value, 2)))))
  ))
}) %>%
  bind_rows()

pval_df <- pval_df %>%
  mutate(
    x     = 1.5,   # midpoint between the two x conditions
    y     = 6.3    # just below your y limit of 2
  )

mhc_plot <- ggplot(full_df_long %>% filter(Condition %in% c("prim_HLA_combined", "ln_met_HLA_combined")), 
  aes(x = Condition, y = Expression, fill = IFN)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    linewidth = 0.5,
    color = "black",
    alpha = 0.75
  ) +
  geom_jitter(
    aes(fill = Condition),
    position = my_jitter,
    color = "black",
    shape = 21,
    size = 1.8,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  # geom_text_repel(
  #   aes(label = annonID), # Add fill here too!
  #   position = my_jitter,
  #   size = 1.5
  # ) +

  geom_line(aes(group = Sample)) +
  geom_text(
    data = pval_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.5
  ) +
  facet_wrap(~ IFN, nrow = 1, scales = "fixed",
             labeller = labeller(IFN = IFN_labels)) +
  scale_y_continuous(limits = c(0.5, 7)) +
  scale_fill_manual(
    values = condition_colors,
    labels = IFN_labels,
    name = "Condition"
  ) +
  scale_x_discrete(labels = condition_labels_x) +
  labs(
    title = "",
    x     = NULL,
    y     = "MHC class I Expression"
  ) +
  theme_classic(base_size = 4) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_blank(),
    axis.title.y     = element_text(size = 7, face = "bold", margin = margin(r = 8)),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(size = 4),
    axis.text.x      = element_text(size = 4),
    legend.title     = element_text(size = 7, face = "bold"),
    legend.text      = element_text(size = 4),
    legend.position  = "bottom",
    legend.key.size  = unit(0.5, "cm"),
    panel.spacing    = unit(0, "lines"),
    plot.margin      = margin(10, 15, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

mhc_plot


ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/mhc.pdf", 
  plot = mhc_plot,
  width = 13.5,    # Adjusted for a quarter page layout with margins
  height = 9.5,    # Keeps a clean, slightly wider landscape aspect ratio
  units = "cm"
)

