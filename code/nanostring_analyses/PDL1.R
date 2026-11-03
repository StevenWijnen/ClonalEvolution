set.seed(1616)
setwd("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/code/nanostring_analyses")
library(dplyr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(ggpubr)
library(cutpointr)
library(rstatix)
library(psych)

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

df <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/old_DATA/sel.bas_2907.csv") %>% dplyr::select(StudyID, IFNLN, PrimIFN)


IFN_IHC_labels <- read_xlsx("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/PDL1_steven_2606.xlsx") %>% 
  select(StudyID, PDL1_IHC = PDL1_binary) %>%
  mutate(PDL1_IHC = factor(PDL1_IHC, levels = c("Low", "High", NA)))

df <- left_join(df, IFN_IHC_labels)



prim_expression_data <- nCounter_expression_primary %>% 
  dplyr::select(PDL1_prim = PDCD1, Trial, StudyID)

ln_met_expression_data <- nCounter_expression_ln_met %>% 
  dplyr::select(PDL1_ln = PDCD1 , Trial, StudyID)

ln_met_expression_data <- left_join(ln_met_expression_data, df)

full_nanostring_expression_mat <- inner_join(prim_expression_data, ln_met_expression_data)
full_nanostring_df <- inner_join(full_nanostring_df, full_nanostring_expression_mat)



lw <- function(pt) pt / (.pt * 0.75)


#Correlations
{
primary <- ggplot(full_nanostring_df, aes(x = PrimIFN, y = PDL1_prim)) +
      theme_classic() +
      geom_point(size = 0.8) +
      geom_smooth(method = "lm", linewidth = lw(0.75)) +
      stat_cor(
        method = "spearman", 
        label.x = -35, 
        label.y = -1,
        r.digits = 2,
        p.accuracy = 0.000001, 
        size = 6 / .pt
        ) +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    text        = element_text(size = 6),
    axis.text   = element_text(size = 6),
    axis.title  = element_text(size = 6),
    axis.line   = element_line(linewidth = lw(0.406)),
    axis.ticks  = element_line(linewidth = lw(0.406))
  ) +     labs(
        x = "IFNy gene score - primary melanoma",
        y = "PD-L1 Expression - primary melanoma"
      ) 



primary



cor.test(
  full_nanostring_df$PrimIFN,
  full_nanostring_df$PDL1_prim,
  method = "spearman",
  conf.level = 0.95
)
ggsave("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/primary_PDL1_IFN_corr.pdf", primary, width = 93.252, height = 65.325, units = "mm")  


lymphnode <- ggplot(full_nanostring_df, aes(x = IFNLN, y = PDL1_ln)) +
      theme_classic() +
      geom_point(size = 0.8) +
      geom_smooth(method = "lm", linewidth = lw(0.75)) +
      stat_cor(
        method = "spearman", 
        label.x = -35, 
        label.y = -1,
        r.digits = 2,
        p.accuracy = 0.000001, 
        size = 6 / .pt
        ) +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    text        = element_text(size = 6),
    axis.text   = element_text(size = 6),
    axis.title  = element_text(size = 6),
    axis.line   = element_line(linewidth = lw(0.406)),
    axis.ticks  = element_line(linewidth = lw(0.406))
  ) + labs(
        x = "IFNy gene score − lymph node metastasis",
        y = "PDL1 Expression − lymph node metastasis"
      )

lymphnode



res_psych <- corr.test(
  full_nanostring_df$PDL1_prim, 
  full_nanostring_df$PDL1_ln, 
  method = "spearman"
)

# View the confidence interval bounds
res_psych$ci$lower
res_psych$ci$upper

ggsave("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/ln_PDL1_IFN_corr.pdf", lymphnode, width = 93.252, height = 65.325, units = "mm")  

# ggsave(
#   filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/ln_PDL1_IFN_corr.pdf", 
#   plot = lymphnode,
# #  width = 97.077,
#   # height = 66.921,
#   # units = "mm"
#   )


pdl1_vs_pdl1 <- ggplot(full_nanostring_df, aes(x = PDL1_prim, y = PDL1_ln)) +
      theme_classic() +
      geom_point(size = 0.8) +
      geom_smooth(method = "lm", linewidth = lw(0.75)) +
      stat_cor(
        method = "spearman", 
        r.digits = 2,
        p.accuracy= 0.01,
        size = 6 / .pt
        ) +
  theme_classic(base_size = 6, base_family = "Helvetica") +
  theme(
    text        = element_text(size = 6),
    axis.text   = element_text(size = 6),
    axis.title  = element_text(size = 6),
    axis.line   = element_line(linewidth = lw(0.406)),
    axis.ticks  = element_line(linewidth = lw(0.406))
  ) + 
      labs(
        x = "PD-L1 - primary melanoma",
        y = "PD-L1 lymph node"
      )
pdl1_vs_pdl1

res <- cor.test(
  full_nanostring_df$PDL1_prim,
  full_nanostring_df$PDL1_ln,
  method = "spearman",
  conf.level = 0.95, 
)
res$p.value

ggsave("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/primary_PDL1_PDL1.pdf", pdl1_vs_pdl1, width = 93.252, height = 65.325, units = "mm")  


ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/primary_PDL1_PDL1.pdf", 
  plot = pdl1_vs_pdl1,
#  width = 97.077,
  # height = 66.921,
  # units = "mm"
)





model <- glm(PDL1_IHC ~ PDL1_ln , data=full_nanostring_df %>% filter(!is.na(PDL1_IHC)), family="binomial")
summary(model)

my_jitter <- position_jitterdodge(
  jitter.width = 0,    # Use a dot
  dodge.width = 0.75,   # Use a dot
  seed = 42
)
filtered_df <- full_nanostring_df %>% 
  filter(!is.na(PDL1_IHC), PDL1_IHC != "NA") %>% 
  mutate(PDL1_IHC = factor(PDL1_IHC, levels=c("Low", "High")))

wilcox.test(PDL1_ln ~ PDL1_IHC, data=filtered_df)

filtered_df %>%
  wilcox_effsize(PDL1_ln ~ PDL1_IHC, ci = TRUE, nboot = 1000, conf.level = 0.95)

IHC_vs_Expression <- filtered_df %>%
 ggplot(aes(x = PDL1_IHC, y = PDL1_ln)) + 
  geom_boxplot(
    outlier.shape = NA,
    color = "black",
    fill = "grey88",
    alpha = 0.75,
    width = 0.5,
    linewidth = 0.6
  ) +
  geom_jitter(
    color = "black",
    fill = "grey40",
    shape = 21,
    size = 1.8,
    alpha = 0.55,
    position = my_jitter,
    show.legend = FALSE
  ) +
  
  # 1. Wilcoxon P-value Only
  stat_compare_means(
    method = "wilcox.test",
    label.x = 1.5,
    label.y = max(filtered_df$PDL1_ln, na.rm = TRUE),
    vjust = -1.8,                  
    size = 5,
    hjust = 0.5
  ) +
  
  stat_cor(
    method = "spearman",
    aes(group = 1, label = ..r.label..), # ..r.label.. explicitly strips out the correlation p-value
    label.x = 1.5,
    label.y = max(filtered_df$PDL1_ln, na.rm = TRUE),
    vjust = -0.3,                  
    size = 5,
    hjust = 0.5
  ) +

  theme_classic() +                   
  theme(
    axis.title   = element_text(size = 18, face = "bold"),
    axis.text    = element_text(size = 14, color = "black"),
    axis.line    = element_line(color = "black", linewidth = 0.6),
    axis.ticks   = element_line(color = "black", linewidth = 0.5),
    panel.grid   = element_blank(),
    plot.margin  = margin(35, 20, 10, 10) 
  ) +
  labs(
    x = "PDL1 IHC",
    y = "PDL1 expression"
  ) +
  
  coord_cartesian(clip = "off")

IHC_vs_Expression

ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/PDL1_IHC_vs_expression.pdf", 
  plot = IHC_vs_Expression,
  width = 9.13,    
  height = 6.22,    
  units = "in"
)




}





#Grouped boxplot
full_df_long <- full_nanostring_df %>%
  pivot_longer(
    cols = c(PDL1_prim, PDL1_ln),
    names_to = "Condition",
    values_to = "Expression"
  ) %>% mutate(Condition = factor(Condition, levels=c("PDL1_prim", "PDL1_ln")))

my_jitter <- position_jitterdodge(
  jitter.width = 0,    # Use a dot
  dodge.width = 0.75,   # Use a dot
  seed = 42
)



pval_df <- lapply(names(IFN_labels), function(n) {
  sub_df <- full_nanostring_df %>%
    filter(IFN == n) %>%
    select(Sample, PDL1_prim, PDL1_ln) %>%
    arrange(Sample)  # ensure pairs are aligned

  test <- t.test(
    sub_df$PDL1_prim,
    sub_df$PDL1_ln,
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
    y     = 0    # just below your y limit of 2
  )

PDL1_grouped_boxplot <- ggplot(full_df_long, aes(x = Condition, y = Expression, fill = IFN)) +
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
  scale_y_continuous(limits = c(-6, 0)) +
  scale_fill_manual(
    values = condition_colors,
    labels = IFN_labels,
    name = "Condition"
  ) +
  scale_x_discrete(labels = condition_labels_x) +
  labs(
    title = "",
    x     = NULL,
    y     = "PDL1 expression"
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

PDL1_grouped_boxplot


ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/PDL1_grouped_boxplot.pdf", 
  plot = PDL1_grouped_boxplot,
  width = 13.5,    # Adjusted for a quarter page layout with margins
  height = 9.5,    # Keeps a clean, slightly wider landscape aspect ratio
  units = "cm"
)
