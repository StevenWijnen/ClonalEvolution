library(dplyr)
library(tidyr)
library(purrr)
library(tibble)



source("~/packages/R_Helpers/enrichment.R")
source("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/code/tree_analyses/load_tree_struct.R")
source("/home/s.wijnen/packages/R_Helpers/Annotate_gene_names.R")

patient_data <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/filtered_cleaned_df.csv")
maf_path <- "/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/output/mafs/purple_mafs/ln_met/mafs/"

df_all <- get_tree_nodes_oncoanalyser()
LL_genes <- df_all %>% filter(IFN_group == "LL", sampleID != "OpACINneo_41", clonality == "subclonal_ln")

filtered_LL_ALL <- list()
for(ptid in unique(LL_genes$sampleID)) {
    print(ptid)
    tid <- patient_data %>% filter(Sample == ptid) %>% pull(DNA_ln_cegat)
    fp <- paste0(maf_path, tid, ".purple.maf")


    maf_data <- read.delim(fp, sep = "\t", comment.char = "#")

    genes <- LL_genes %>% filter(sampleID == ptid) %>% pull(gene_name)
    
    genes <- maf_data %>% 
    filter(
        Hugo_Symbol %in% genes & 
        (IMPACT %in% c("HIGH", "MODERATE", "MODIFIER") | Consequence == "splice_region_variant")
    ) %>% 
    pull(Hugo_Symbol) %>% unique()
    filtered_LL_ALL[[ptid]] <- genes
}


filtered_LL_genes <- unique(unlist(filtered_LL_ALL))
go_res_LL  <- enrichment_go(filtered_LL_genes, pvalueCutoff = 0.8, qvalueCutoff=0.8, ontology="MF")@result
reactome_res_LL <- enrichment_reactome(filtered_LL_ALL)@result

ECM_genes_LL <- unlist(strsplit(reactome_res_LL["R-HSA-1474244", "geneID"], "/"))


ln_seeding_nodes <- seeding_nodes_onco()
dd_ALL <- get_genes_in_nodes(ln_seeding_nodes)

filtered_dd_ALL <- list()
for(ptid in names(dd_ALL)) {
    tid <- patient_data %>% filter(Sample == ptid) %>% pull(DNA_ln_cegat)
    fp <- paste0(maf_path, tid, ".purple.maf")


    maf_data <- read.delim(fp, sep = "\t", comment.char = "#")

    genes <- dd_ALL[[ptid]]
    genes <- maf_data %>% 
    filter(
        Hugo_Symbol %in% genes & 
        (IMPACT %in% c("HIGH", "MODERATE", "MODIFIER") | Consequence == "splice_region_variant")
    ) %>% 
    pull(Hugo_Symbol) %>% unique()
    filtered_dd_ALL[[ptid]] <- genes
}

filtered_seeding_node_genes <- setdiff(unique(unlist(filtered_dd_ALL)), filtered_LL_genes)
go_res_ALL  <- enrichment_go(unique(unlist(filtered_dd_ALL)), pvalueCutoff = 0.8, qvalueCutoff=0.8, ontology="MF")@result
reactome_res_ALL <- enrichment_reactome(unique(unlist(filtered_dd_ALL)))@result
msig_c2 <- enrichment_msigdbr(unique(unlist(filtered_dd_ALL)), category="C2")@result

per_sample_res <- imap(filtered_dd_ALL, function(genes, sample_name) {
  
  res <- tryCatch({
    enrichment_reactome(genes)@result
  }, error = function(e) {
    message("Failed for sample: ", sample_name, " — ", conditionMessage(e))
    return(NULL)
  })
  
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  res %>%
    filter(pvalue < 0.2) %>%
    dplyr::select(ID, Description, pvalue, geneID, Count) %>%
    mutate(sample = sample_name)
})

## drop samples that returned nothing
per_sample_res <- compact(per_sample_res)

## 2. Combine into one long table
long_df <- bind_rows(per_sample_res)

pval_matrix <- long_df %>%
  dplyr::select(ID, Description, sample, pvalue) %>%
  pivot_wider(names_from = sample, values_from = pvalue)

## 3b. binary hit matrix (1 = pathway significant in that sample, 0/NA otherwise)
hit_matrix <- long_df %>%
  distinct(ID, Description, sample) %>%
  mutate(hit = 1) %>%
  pivot_wider(names_from = sample, values_from = hit, values_fill = 0)

## 4. Recurrence count per pathway
recurrence <- hit_matrix %>%
  rowwise() %>%
  mutate(n_samples = sum(c_across(-c(ID, Description)))) %>%
  ungroup() %>%
  arrange(desc(n_samples))

## 5. Filter: pathway recurrent in more than 5 samples
recurrent_pathways <- recurrence %>%
  filter(n_samples >= 5)

recurrent_pathways %>% dplyr::select(ID, Description, n_samples)


View(go_res_ALL)
View(reactome_res_ALL)
View(msig_c2)

ECM <- unique(unlist(strsplit(reactome_res_ALL[c("R-HSA-1474244") , "geneID"],  "/")))
NCAM1 <- unique(unlist(strsplit(reactome_res_ALL[c("R-HSA-419037") , "geneID"],  "/")))
ECM_pro <- unique(unlist(strsplit(reactome_res_ALL[c("R-HSA-3000178") , "geneID"],  "/")))
Integrin <- unique(unlist(strsplit(reactome_res_ALL[c("R-HSA-216083") , "geneID"],  "/")))
PDGF <- unique(unlist(strsplit(reactome_res_ALL[c("R-HSA-186797") , "geneID"],  "/")))

genes_of_interest <- unique(c(ECM, NCAM1, ECM_pro, Integrin, PDGF))

genes_motoractivity_lh <- unique(unlist(strsplit(go_res_LH[c(1:7),"geneID"],"/")))
go_ids_LH <- go_res_LH[c(1:7),"ID"]

sialyiation <-  unique(unlist(strsplit(reactome_res_ALL[reactome_res_ALL$ID %in% c("R-HSA-3000178", "R-HSA-419037", "R-HSA-375165", "R-HSA-3906995", "O-linked glycosylation") ,"geneID"], "/")))
 


cytoskeletal_motor <- strsplit(go_res[go_res$ID == "GO:0003774","geneID"], "/")[[1]]
ecm <- strsplit(go_res[go_res$ID == "GO:0005201","geneID"], "/")[[1]]
total = 0 
for(pt_id in names(filtered_dd_ALL)) {
    if(any(filtered_dd_ALL[[pt_id]] %in% ECM   )) {
        print(pt_id)
    } 
#  if (any(dd_ALL[[pt_id]] %in% ecm)) {
#     print("ecm")
#        print(pt_id)
#     }
}








ln_seeding_nodes_LH <- seeding_nodes_onco(LH_only=TRUE)
dd_LH <- get_genes_in_nodes(ln_seeding_nodes_LH)
go_res_LH  <- enrichment_go(unique(unlist(dd_LH)), ontology="MF")@result
reactome_res_LH <- enrichment_reactome(unique(unlist(dd_LH)))@result
genes_motoractivity_lh <- unique(unlist(strsplit(go_res_LH[c(1:7),"geneID"],"/")))
go_ids_LH <- go_res_LH[c(1:7),"ID"]



for(pt_id in names(dd_LH)) {
    if(any(dd_LH[[pt_id]] %in% genes_motoractivity_lh   )) {
        print(dd_LH[[pt_id]][which(dd_LH[[pt_id]] %in% genes_motoractivity_lh   )])
        print(pt_id)
    }
}


View(go_res_hl)


View(go_res_LH)
View(reactome_res_LH)


dnds_genes$gene[which(dnds_genes$gene %in%  unlist(strsplit(go_res[9,"geneID"], "/")))]


hl_ln_nodes <- HL_pts_all_ln_nodes()
hl_ALL <- get_genes_in_nodes(hl_ln_nodes)

go_res_hl  <- enrichment_go(unique(unlist(hl_ALL)), ontology="MF")@result
reactome_res_hl <- enrichment_reactome(unique(unlist(hl_ALL)))@result
# msig_c2_hl <- enrichment_msigdbr(unique(unlist(hl_ALL)), category="C2")@result


ll_ln_nodes <- LL_pts_all_ln_nodes()
ll_ALL <- get_genes_in_nodes(ll_ln_nodes)

go_res_ll  <- enrichment_go(unique(unlist(ll_ALL)), ontology="MF")@result
reactome_res_ll <- enrichment_reactome(unique(unlist(ll_ALL)))@result

View(go_res_ll)
genes_minus_end_ll <- unique(unlist(strsplit(go_res_ll[go_ids_LH,"geneID"],"/")))

IMPORTANT_GENES <- setdiff(genes_motoractivity_lh , genes_minus_end_ll)


non_seeding_node_sample_genes <- get_all_nodes() %>% 
    filter(!sampleID %in% names(dd_ALL), clonality == "subclonal_ln") %>% 
    filter(IFN_group != "HH") %>%
    pull(gene_name) %>%
    unique()
go_res_non_seeding  <- enrichment_go(non_seeding_node_sample_genes, ontology="MF")@result
non_seeding_cytoskeletal <- unlist(strsplit(go_res_non_seeding["GO:0003774", "geneID"], "/"))


sum(cytoskeletal_motor %in% non_seeding_cytoskeletal)
