#Using Oncoanalyser pave and purple 
library(dplyr)
library(maftools)
library(ggplot2)

setwd("/DATA/s.wijnen/Projects/ClonalEvolution/clean_code/ClonalEvolution/code")
source("mutation_matrix/build_primary_maf_df.R")
source("mutation_matrix/build_ln_met_maf_df.R")
source("mutation_matrix/drawownoncoplot.R")

clinical_data <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/filtered_cleaned_df.csv")
anonymous_ids <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/Anonymous_ptid_mapper.csv")

clinical_data <- clinical_data %>% left_join(anonymous_ids %>% dplyr::select(sampleID,annonID), by=c("Sample" = "sampleID"))
clinical_data <- clinical_data %>% mutate(IFN = case_when(
    IFN == "primhigh_LNhigh" ~ "HH",
    IFN == "primhigh_LNlow" ~ "HL",
    IFN == "primlow_LNhigh" ~ "LH",
    IFN == "primlow_LNlow" ~ "LL",
))

clinical_features <- c("Trial", "MPR", "EFSstat", "IFN", "Gender", "type", "Tstage", "primary_site")
prim_clinical_data <- clinical_data
prim_clinical_data$Tumor_Sample_Barcode <- clinical_data$DNA_primary_cegat
prim_clinical_data$annonID <- paste0("P", prim_clinical_data$annonID)

ls_driver_genes <- c(
  "BRAF" , "NRAS", "KIT", "NF1",
  "ARID2",  "CDKN2A", "CTNNB1", "DCC", 'DPYD',
  'EZH2','FGFR2', 'GRM3', 'KCNQ5', 'KDR' , 'MAP2K1', 'MET', 
  'PDE1A', 'PTEN', 'RB1', 'TP53', 'TYRP1', "COL2A1",   "COL5A2", "PKHD1L1", "UNC13C", "WNK1", #, "BRAF-non-V600" , "BRAFV600E"
   "SNX31",   "ARID2"
)

#These genes can be removed: STK19 TACC1,RAC1,PPP6C,CD274,SF3B1, LCK,GHR,CHD8,CDK4,TYRP1§

primary_maf_df <- build_primary_maf_df_hmf()
maf_datastructure_primary <- maftools::read.maf(maf = primary_maf_df, clinicalData=prim_clinical_data)

ln_met_clinical_data <- clinical_data
ln_met_clinical_data$Tumor_Sample_Barcode <- clinical_data$DNA_ln_cegat
ln_met_clinical_data$annonID <- paste0("LN", ln_met_clinical_data$annonID)

ln_met_maf_df <- build_ln_met_maf_df_hmf()
maf_datastructure_ln_met <- maftools::read.maf(maf = ln_met_maf_df, clinicalData=ln_met_clinical_data)
df <-somaticInteractions(maf = maf_datastructure_ln_met, genes=ls_driver_genes,  pvalue = c(0.05, 0.01))

combined_data_structure <- rbind(primary_maf_df, ln_met_maf_df)
long_clinical_data <- clinical_data %>%
  pivot_longer(
    cols = c(DNA_primary_cegat, DNA_ln_cegat),
    names_to = "source",
    values_to = "Tumor_Sample_Barcode"
  ) %>%
  mutate(
    source = if_else(source == "DNA_primary_cegat", "primary", "ln_met"),    
    name = if_else(source == "primary", "case_id_prim", "case_id_ln_met")) %>%
  dplyr::select(case_id=Sample, name, Tumor_Sample_Barcode, source, everything())






maf_ds_combined <- maftools::read.maf(maf = combined_data_structure, long_clinical_data)
# df <-somaticInteractions(maf = maf_ds_combined, genes = c("BRAF", "STAB1", "CSMD2", "NRAS", "COL19A1", "MROH2B", "DNAH6"), pvalue = c(0.05, 0.1))
driver_genes <- read.csv("/data_pool_1/ClonalEvolution/oncoanalyser/input/custom_panel_self/DriverGenePanel.38.tsv", sep="\t")

driver_gene_names <- c(driver_genes$gene, ls_driver_genes)
gene_counts <- maf_ds_combined@data[
  Hugo_Symbol %in% driver_gene_names, 
  .(n_samples = data.table::uniqueN(Tumor_Sample_Barcode)), 
  by = Hugo_Symbol
]



genes_to_show <- gene_counts[n_samples > 4, Hugo_Symbol]
#genes_to_show <- c(genes_to_show, "TERT") #Christian wanted TERT added TERT was only present in one sample
somaticInteractions(maf = maf_datastructure_ln_met, genes=genes_to_show,  pvalue = c(0.1, 0.05))


publication_oncoplot(maf_ds_combined, genes=genes_to_show, show_column_names = TRUE, output="../output/figures/mut_matrix.pdf")



primary_ <- combined_data_structure %>% filter(origin == "primary")
ln_met_ <- combined_data_structure %>% filter(origin == "ln_met")

maf_obj_primary_ <- maftools::read.maf(maf=primary_)
maf_obj_ln_met_ <- maftools::read.maf(maf=ln_met_)

res <- mafCompare(maf_obj_primary_, maf_obj_ln_met_, m1Name = "primary", minMut=5, m2Name = "ln_met", useCNV = FALSE)
View(res$results)


combined_data_structure <- combined_data_structure %>% left_join( clinical_data %>% dplyr::select(Sample, IFN),by=c("sample" = "Sample"))


ls_all_results <- list()
for(group in c("LL", "LH", "HL", "HH")) {
  df_maf <- combined_data_structure %>% filter(IFN == group)
  primary_ <- df_maf %>% filter(origin == "primary")
  ln_met_ <- df_maf %>% filter(origin == "ln_met")

  maf_obj_primary_ <- maftools::read.maf(maf=primary_)
  maf_obj_ln_met_ <- maftools::read.maf(maf=ln_met_)


  res <- mafCompare(maf_obj_primary_, maf_obj_ln_met_, m1Name = "primary", minMut=3, m2Name = "ln_met", useCNV = FALSE)

  ls_all_results[[group]] <- res
}


maf1 <- maftools::read.maf(combined_data_structure %>% filter(origin == "ln_met") %>% filter(IFN %in% c("HH", "LH")))
maf2 <- maftools::read.maf(combined_data_structure %>% filter(origin == "ln_met") %>% filter(IFN %in% c("LL", "HL")))



test<-mafCompare(maf1, maf2, m1Name = "HH", minMut=3, m2Name = "LL", useCNV = FALSE)
