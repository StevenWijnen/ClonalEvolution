full_nanostring_df <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/cleaned_full_df.csv") %>% dplyr::select(Sample, StudyID, Trial, IFN, IFNln)
anon_ids <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/output/Anonymous_ptid_mapper.csv")
full_nanostring_df <- left_join(full_nanostring_df, anon_ids %>% dplyr::select(sampleID, annonID), by=c("Sample"= "sampleID"))

nCounter_expression_primary <- read.csv( "/DATA/s.wijnen/Projects/ClonalEvolution/DATA/nCounter/primary_nCounter_processed.csv")
nCounter_expression_ln_met <- read.csv("/DATA/s.wijnen/Projects/ClonalEvolution/DATA/nCounter/ln_met_nCounter_processed.csv")


