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
  "prim_jak_stat" = "Primary",
  "ln_met_jak_stat" = "Lymph node"
)




# jak_stat_ifng_mhc1 <- c(
#   # IFN-gamma receptor complex
#   "IFNGR1", "JAK1", "JAK2", "STAT1", "IRF2", "STAT5B", "TAP1", "TAP2"
# )
expended_jak_stat <- c("IFNGR1", "JAK1", "JAK2", "JAK3", "PTPRC", "SOCS1", "STAT1", "STAT3", "STAT4", "STAT5B", "STAT6")

prim_expression_data <- nCounter_expression_primary %>% 
  mutate(prim_jak_stat = rowMeans(dplyr::select(., all_of(expended_jak_stat)), na.rm = TRUE)) %>%
  dplyr::select(prim_jak_stat, Trial, StudyID)

ln_met_expression_data <- nCounter_expression_ln_met %>% 
  mutate(ln_met_jak_stat = rowMeans(dplyr::select(., all_of(expended_jak_stat)), na.rm = TRUE)) %>%
  dplyr::select(ln_met_jak_stat,  all_of(expended_jak_stat) , Trial, StudyID)


full_nanostring_expression_mat <- inner_join(prim_expression_data, ln_met_expression_data)
full_nanostring_df <- inner_join(full_nanostring_df, full_nanostring_expression_mat)


full_nanostring_df
full_df_long <- full_nanostring_df %>%
  pivot_longer(
    cols = c(prim_jak_stat, ln_met_jak_stat),
    names_to = "Condition",
    values_to = "Expression"
  ) %>% mutate(Condition = factor(Condition, levels=c("prim_jak_stat", "ln_met_jak_stat")))

my_jitter <- position_jitterdodge(
  jitter.width = 0,    # Use a dot
  dodge.width = 0.75,   # Use a dot
  seed = 42
)

#SIMPLE BOXPLOT 
ggplot(full_df_long, aes(x = IFN, y = Expression, fill = Condition)) + 
  geom_boxplot(outlier.shape = NA) + 
  
  # Crucial: Add 'fill = HLA_Type' here so it knows which box to align with
  geom_jitter(
    aes(fill = Condition), 
    position = my_jitter, 
    color = "black",  
    shape = 21,       
    alpha = 0.5,
    show.legend = FALSE  # Keeps the legend from getting 'dotted'
  ) +  scale_fill_discrete(
    name = "Site",
    labels = c("prim_jak_stat" = "Primary", "ln_met_jak_stat" = "Lymph Node")
  )+  geom_text_repel(
    aes(label = StudyID), # Add fill here too!
    position = my_jitter,
    size = 3
  ) +
  
  theme(
    axis.title = element_text(size = 18, face = "bold"), # X and Y axis labels
    axis.text = element_text(size = 14),                # Axis tick marks
    legend.title = element_text(size = 16),             # Legend title
    legend.text = element_text(size = 14)               # Legend items
  ) + 
  labs(
    title = "",
    x = "IFN group",
    y = "Jak/stat pathway expression"
  )

ggplot(full_nanostring_df, aes(x = IFNln, y = JAK1)) + 
  geom_boxplot(outlier.shape = NA) +  geom_text_repel(aes(label=StudyID))
  


pval_df <- lapply(names(IFN_labels), function(n) {
  sub_df <- full_nanostring_df %>%
    filter(IFN == n) %>%
    select(Sample, prim_jak_stat, ln_met_jak_stat) %>%
    arrange(Sample)  # ensure pairs are aligned

  test <- t.test(
    sub_df$prim_jak_stat,
    sub_df$ln_met_jak_stat,
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
    y     = 1.9    # just below your y limit of 2
  )

jak_stat_plot <- ggplot(full_df_long, aes(x = Condition, y = Expression, fill = IFN)) +
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
  scale_y_continuous(limits = c(-0.8, 2)) +
  scale_fill_manual(
    values = condition_colors,
    labels = IFN_labels,
    name = "Condition"
  ) +
  scale_x_discrete(labels = condition_labels_x) +
  labs(
    title = "",
    x     = NULL,
    y     = "JAK/STAT Pathway Expression"
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

jak_stat_plot





ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/jakstat.pdf", 
  plot = jak_stat_plot,
  width = 13.5,    # Adjusted for a quarter page layout with margins
  height = 9.5,    # Keeps a clean, slightly wider landscape aspect ratio
  units = "cm"
)

results <- lapply(expended_jak_stat, function(gene) {
  ff <- as.formula(paste0(gene, " ~ IFNln"))
  
  test <- t.test(ff, data = full_nanostring_df)
  
  data.frame(
    gene = gene,
    pval = test$p.value,
    W    = test$statistic
  )
}) %>%
  bind_rows() %>%
  mutate(padj = p.adjust(pval, method = "BH"))

print(results)

