#dNdS analyses 
library(dndscv)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(scales)
source("~/packages/R_Helpers/Annotate_gene_names.R")
source("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/code/tree_analyses/load_tree_struct.R")



df_all <- get_tree_nodes_oncoanalyser()
df_all <- build_df(df_all)
genes_in_more_than_two_samples <- df_all %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(gene_name) %>%
  summarise(n_samples = n_distinct(sampleID)) %>%
  filter(n_samples > 2) %>%
  pull(gene_name)


genes_in_more_than_two_samples <- c(genes_in_more_than_two_samples, "CDKN2A.p16INK4a", "CDKN2A.p14arf")

all_dnds <- dndscv(df_all, refdb="hg38", max_muts_per_gene_per_sample = 5, max_coding_muts_per_sample=Inf, outmats=T)
truncal_dnds_genes = dndscv(df_all %>% filter(clonality =="clonal"), refdb="hg38", max_muts_per_gene_per_sample = 5, max_coding_muts_per_sample=Inf, outmats=T)
subclonal_dnds_genes = dndscv(df_all %>% filter(clonality == "subclonal_ln"), refdb="hg38", max_muts_per_gene_per_sample = 5, max_coding_muts_per_sample=Inf, outmats=T)
subclonal_primary_dnds_genes = dndscv(df_all %>% filter(clonality == "subclonal_primary"), refdb="hg38", max_muts_per_gene_per_sample = 5, max_coding_muts_per_sample=Inf, outmats=T)


truncal_genes <- truncal_dnds_genes$sel_cv %>% 
  filter(pallsubs_cv < 0.01, gene_name %in% genes_in_more_than_two_samples) %>% 
  arrange(pallsubs_cv) %>%       # Removed desc() to sort lowest-to-highest
  slice_min(pallsubs_cv, n = 20) %>%  # Safely extracts the 10 lowest values
  pull(gene_name)

subclonal_ln_genes <- subclonal_dnds_genes$sel_cv %>% 
  filter(pallsubs_cv < 0.01, gene_name %in% genes_in_more_than_two_samples) %>% 
  arrange(pallsubs_cv) %>%       # Removed desc() to sort lowest-to-highest
  slice_min(pallsubs_cv, n = 20) %>%  # Safely extracts the 10 lowest values
  pull(gene_name)

subclonal_primary_genes <- subclonal_primary_dnds_genes$sel_cv %>% 
  filter(pallsubs_cv < 0.01, gene_name %in% c(genes_in_more_than_two_samples, "PLD4")) %>% 
  arrange(pallsubs_cv) %>%       # Removed desc() to sort lowest-to-highest
  slice_min(pallsubs_cv, n = 10) %>%  # Safely extracts the 10 lowest values
  pull(gene_name)

gene_list_ln <- c(truncal_genes, subclonal_ln_genes)
gene_list_ln <- gene_list_ln[!duplicated(gene_list_ln)]

gene_list_prim <- c(truncal_genes, subclonal_primary_genes)
gene_list_prim <- gene_list_prim[!duplicated(gene_list_prim)]


ci_truncal <- geneci(truncal_dnds_genes, gene_list = c(gene_list_ln, gene_list_prim ))
ci_sub <- geneci(subclonal_dnds_genes, gene_list = gene_list_ln)
ci_sub_primary <- geneci(subclonal_primary_dnds_genes, gene_list = gene_list_prim)



#Plotting stuff
{
  #Function to combined the w scores of all types of mutations
compute_weighted_dnds <- function(dndsout, ci_obj, label, gene_list) {
    weights <- dndsout$genemuts %>%
      mutate(n_total_point = n_mis + n_non + n_spl) %>%
      dplyr::select(gene_name, exp_mis, exp_non, exp_spl, n_total_point) %>% 
      filter(gene_name %in% gene_list)
    

    
    pts <- dndsout$sel_cv %>%
      dplyr::select(gene_name, wmis_cv, wnon_cv, wspl_cv) #fitted parameters
    
    combined <- weights %>%
      left_join(pts, by = "gene_name") %>%
      left_join(ci_obj, by = c("gene_name" ="gene")) %>%
      mutate(
        exp_tru = exp_non + exp_spl,
        w_combined = (wmis_cv * exp_mis + wnon_cv * exp_non + wspl_cv * exp_spl) /
          (exp_mis + exp_non + exp_spl),
        ci_low_combined = (mis_low * exp_mis + tru_low * exp_tru) /
          (exp_mis + exp_tru),
        ci_high_combined = (mis_high * exp_mis + tru_high * exp_tru) /
          (exp_mis + exp_tru)
      ) %>%
      dplyr::select(gene_name, w_combined, n_total_point, ci_low_combined, ci_high_combined)
    
      
    out <- combined %>%
      dplyr::select(gene = gene_name, 
                    w_comb = w_combined,
                    ci_low = ci_low_combined,
                    ci_high = ci_high_combined, 
                    total = n_total_point
                    )
    
    colnames(out)[-1] <- paste0(colnames(out)[-1], "_", label)
    
    return(out)
  }
  
  
  trunc_df <- compute_weighted_dnds(truncal_dnds_genes, ci_truncal, "t", gene_list_prim )
  sub_df   <- compute_weighted_dnds(subclonal_primary_dnds_genes, ci_sub_primary, "s", gene_list_prim)

  df_plot <- tibble(gene = gene_list_prim) %>%
    left_join(trunc_df, by = "gene") %>%
    left_join(sub_df,   by = "gene") 
   #%>%
    # Clamp lower CIs to 0
    # mutate(across(c(ci_low_t, ci_low_s), ~ pmax(.x, 0)))
  
  # ---- 4. Colour by selection quadrant ----
  df_plot <- df_plot %>%
    mutate(
      selection = case_when(
        w_comb_t > 1 & w_comb_s > 1 ~ "Both",
        w_comb_t > 1               ~ "Truncal only",
        w_comb_s > 1               ~ "Subclonal only",
        TRUE                      ~ "Neither"
      )
    )
  
  sel_colours <- c(
    "Both"           = "#E63946",
    "Truncal only"   = "#2176AE",
    "Subclonal only" = "#F4A261",
    "Neither"        = "grey60"
  )
  df_plot$gene <- gsub("CDKN2A.p16INK4a", "CDKN2A", df_plot$gene)
  # df_plot<-df_plot %>% filter(gene == "PLD4")
  # ---- 5. Plot ----
  p <- ggplot(df_plot, aes(x = w_comb_t, y = w_comb_s , colour = selection)) +
    
    # Reference lines
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
    geom_point(aes(size=total_t + total_s), alpha = 0.92) +

    # CI bars
    geom_errorbar(aes(ymin = ci_low_s, ymax = ci_high_s),
                  width = 0, linewidth = 0.25, alpha = 0.35) +
    geom_errorbarh(aes(xmin = ci_low_t, xmax = ci_high_t),
                   height = 0, linewidth = 0.25, alpha = 0.35) +

    # Points
    ggrepel::geom_text_repel(
          data = ~ subset(., selection == "Truncal only"),
          aes(label = gene),
          show.legend        = FALSE,
          size               = 2,
          fontface           = "italic",
          angle              = 45,           # steeper angle gives more room
          direction          = "x",          
          nudge_y            = -0.4,         # push further down
          vjust              = 1,            
          min.segment.length = 0.2,
          segment.size       = 0.3,
          segment.color      = "grey50",
          segment.linetype   = "dashed",
          box.padding        = 0.8,          # more padding
          point.padding      = 0.4,
          max.overlaps       = Inf,
          force              = 5,            # stronger repulsion
          seed               = 42
        ) +
    # # Gene labels
    ggrepel::geom_text_repel(
      data = ~ subset(., selection != "Truncal only"),
      aes(label = gene),
      show.legend        = FALSE,

      size               = 2,
      fontface           = "italic",
      nudge_x            = -0.085,
      # direction          = "y",
      # hjust              = 1,
      min.segment.length = 0.2,
      segment.size       = 0.3,
      segment.color      = "grey50",
      box.padding        = 0.4,
      point.padding      = 0.3,
      max.overlaps       = Inf,
      seed               = 42
    ) +  scale_colour_manual(values = sel_colours, name = "Selection") +
    scale_x_continuous(
      trans = pseudo_log_trans(base = 10),
      breaks = c(0, 1, 10, 100),
      name   = expression("Truncal dN/dS")
    ) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10),
                       breaks = c(0, 1, 10, 100),
                       name   = expression("Lymph node subclonal dN/dS")
                       ) +
    coord_cartesian(ylim = c(-2, 700), clip = "off") +
    scale_size_continuous(name = "Total mutations") +
    theme_classic(base_size = 4) +
    theme(
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      axis.line       = element_blank(),
      axis.title      = element_text(size = 7),
      axis.text       = element_text(size = 4, colour = "black"),
      axis.ticks      = element_line(colour = "black", linewidth = 0.4),
      legend.position = "right",
      legend.title    = element_text(size = 7, face = "bold"),
      legend.text     = element_text(size = 4),
      panel.grid      = element_blank()
  ) 

  
  print(p)
}

# write.csv(sub_df %>% filter(w_comb_s > 1), "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/genes_of_interest.csv")
ggsave(
  filename = "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/figures/dNdS_biggeronco_PRIM.pdf", 
  plot = p,
  width = 13.5,    # Adjusted for a quarter page layout with margins
  height = 9.5,    # Keeps a clean, slightly wider landscape aspect ratio
  units = "cm"
)




df_truncal <- df_all %>% filter(gene_name %in% gene_list_ln)
all_samples <- df_truncal %>% distinct(sampleID, IFN_group)
genes <- unique(df_truncal$gene_name)

results <- lapply(genes, function(gene) {
  mutated_samples <- df_truncal %>%
    filter(gene_name == gene) %>%
    pull(sampleID) %>%
    unique()
  
  tbl <- all_samples %>%
    mutate(mutated = ifelse(sampleID %in% mutated_samples, "mutated", "not_mutated")) %>%
    count(IFN_group, mutated) %>%
    pivot_wider(names_from = mutated, values_from = n, values_fill = 0)
  
  ct <- tbl %>% dplyr::select(-IFN_group) %>% as.matrix()
  rownames(ct) <- tbl$IFN_group
  print(ct)
  test <- tryCatch(fisher.test(ct), error = function(e) NULL)
  if (is.null(test)) return(NULL)
  
  data.frame(gene = gene, pvalue = test$p.value)
})

results_df <- bind_rows(results)

results_df <- results_df %>%
  mutate(pvalue_adj = p.adjust(pvalue, method = "BH")) %>%  # Benjamini-Hochberg
  arrange(pvalue_adj)

results_df %>% filter(pvalue_adj < 0.2)

table(df_truncal$IFN_group, df_truncal$gene_name == "PON3")

mutated_samples <- df_truncal %>%
    filter(gene_name %in% c("PAX2","NELL1", "NELL1","SHISA3", "FLII" ) ) %>%
    pull(sampleID) %>%
    unique()

tbl <- all_samples %>%
    mutate(mutated = ifelse(sampleID %in% mutated_samples, "mutated", "not_mutated")) %>%
    count(IFN_group, mutated) %>%
    pivot_wider(names_from = mutated, values_from = n, values_fill = 0)
tbl





