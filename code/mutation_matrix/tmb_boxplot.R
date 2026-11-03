library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(tidyr)

anonymous_ids <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/Anonymous_ptid_mapper.csv")

pave_tmb_own_primary <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/primary/tmb/pave_tmb.csv")
pave_tmb_own_ln_met <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/ln_met/tmb/pave_tmb.csv")



anonymous_ids <- anonymous_ids %>% 
  inner_join(pave_tmb_own_primary, by=c("sampleID"= "case_id")) %>%
  left_join(pave_tmb_own_ln_met, by=c("sampleID"= "case_id"), suffix=c("_prim", "_ln_met"))


condition_colors <- c(
  "HH" = "#8b1812",
  "HL"  = "#79a8da",
  "LH"  = "#dea562",
  "LL"   = "#35398e"
)

IFN_labels <- list(
  "HH" = "Primary High / LN High",
  "HL"  = "Primary High / LN Low",
  "LH"  = "Primary Low / LN High",
  "LL"   = "Primary Low / LN Low"
)

long_anonymous_ids <- anonymous_ids %>%
  pivot_longer(
    cols = c(tmb_per_mb_prim, tmb_per_mb_ln_met),
    names_to = "Condition",
    values_to = "TMB"
  ) %>% mutate(Condition = factor(Condition, levels=c("tmb_per_mb_prim", "tmb_per_mb_ln_met")))

condition_labels_x <- c(
  "tmb_per_mb_prim" = "Primary",
  "tmb_per_mb_ln_met" = "Lymph node"
)

my_jitter <- position_jitterdodge(
  jitter.width = 0,    # Use a dot
  dodge.width = 0.75,   # Use a dot
  seed = 42
)


plot_tmb <- long_anonymous_ids %>% 
  ggplot(aes(x = Condition, y = TMB, fill = IFN_short)) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.6,
      linewidth = 0.5,
      color = "black",
      alpha = 0.75
    ) +
    scale_fill_manual(
      values = condition_colors,
      labels = IFN_labels,                                      # ← apply labels here
      name = ""
    ) +
#   geom_text_repel(
#     aes(label = annonID, group = IFN_short),    # ← add group
#     size = 1.5,
# ) +   
  geom_jitter(
    aes(fill = Condition), 
    position = my_jitter, 
    color = "black",  
    shape = 21,       
    alpha = 0.5,
    size = 1.5,
    show.legend = FALSE  # Keeps the legend from getting 'dotted'
  ) + 
geom_line(aes(group = sampleID)) +

ggnewscale::new_scale_fill() +  
facet_wrap(~ IFN_short, nrow = 1, scales = "fixed",
             labeller = labeller(IFN = IFN_labels)) +
               scale_fill_manual(
    values = condition_colors,
    labels = IFN_labels,
    name = "Condition"
  ) +
  scale_x_discrete(labels = condition_labels_x) +
  labs(
    title = "",
    x     = NULL,
    y     = "TMB (mut/mb)"
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



ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/tmb_boxplot.pdf", 
  plot = plot_tmb,
  width = 13.5,    # Adjusted for a quarter page layout with margins
  height = 9.5,    # Keeps a clean, slightly wider landscape aspect ratio
  units = "cm"
)
    

  