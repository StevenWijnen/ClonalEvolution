library(ggpubr)
library(tidyr)
library(ggplot2)

source("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/code/tree_analyses/load_tree_struct.R")

df <- get_topologies_oncoanalyser()
table(df$IFN_group, df$seeding_node)
fisher.test(table(df$IFN_group, df$seeding_node))

table(df$IFN_ln, df$seeding_node)
fisher.test(table(df$IFN_ln, df$seeding_node))


table(df$IFN_group, df$shared_prim_ln_subclone)
fisher.test(table(df$IFN_group, df$shared_prim_ln_subclone))

table(df$IFN_group, df$TMB_ln)
chisq.test(table(df$IFN_group, df$shared_prim_ln_subclone))

table(df$IFN_group, df$unused_cluster )
fisher.test(table(df$IFN_group, df$unused_cluster ))

table(df$IFN_ln, df$ln_clonality )
fisher.test(table(df$IFN_ln, df$ln_clonality ))

kruskal.test(IFN_group ~ number_of_nodes, df)


# Data from your table
data <- data.frame(
  IFN_group = c("HH", "HL", "LH", "LL"),
  FALSE_count = c(4, 4, 3, 5),
  TRUE_count  = c(5, 0, 5, 0)
)

# Calculate percentages
data$total <- data$FALSE_count + data$TRUE_count
data$pct_FALSE <- data$FALSE_count / data$total * 100
data$pct_TRUE  <- data$TRUE_count  / data$total * 100

data_long <- pivot_longer(data,
                          cols = c(pct_FALSE, pct_TRUE),
                          names_to = "seeding_node",
                          values_to = "percentage"
)
data_long$seeding_node <- gsub("pct_", "", data_long$seeding_node)

n_labels <- data %>%
  dplyr::select(IFN_group, total) %>%
  distinct() %>%
  mutate(label = paste0(IFN_group, "\n(n=", total, ")")) %>%
  pull(label)

seeding_node_plot<- ggplot(data_long, aes(x = IFN_group, y = percentage, fill = seeding_node)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 2, color = "white", fontface = "bold") +
  scale_fill_manual(
    values = c("FALSE" = "#4A90D9", "TRUE" = "#E05C5C"),
    labels = c("FALSE" = "No seeding node", "TRUE" = "Seeding node")
  ) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_discrete(labels = n_labels) +   # <-- adds n= to x labels
  labs(
    title = "",
    x     = "IFN Group",
    y     = "Percentage",
    fill  = "Seeding node"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title         = element_text(face = "bold"),
    legend.position    = "right",
    panel.grid.major.x = element_blank()
  )
seeding_node_plot


ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/seeding_node_percentages.pdf", 
  plot = seeding_node_plot,
  width = 13.5,    # Adjusted for a quarter page layout with margins
  height = 9.5,    # Keeps a clean, slightly wider landscape aspect ratio
  units = "cm"
)






ggplot(df, aes(x = IFN_ln, y = number_of_private_primary_nodes, fill = IFN_ln, color = IFN_ln)) +
  
  geom_violin(alpha = 0.25, linewidth = 0.7, trim = FALSE) +
  
  geom_boxplot(width = 0.12, alpha = 0.9, outlier.shape = NA,
               color = "black", linewidth = 0.5) +
  
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.4) +
  scale_fill_manual(values  = c("#7F77DD", "#1D9E75", "#D85A30", "#D4537E")) +
  scale_color_manual(values = c("#7F77DD", "#1D9E75", "#D85A30", "#D4537E")) +
  labs(
    x = NULL,
    y = "Number of nodes"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

