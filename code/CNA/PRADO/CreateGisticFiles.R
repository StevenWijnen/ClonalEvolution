library(dplyr)
library(readxl)
library(readr)
library(tidyr)
library(survminer)
library(survival)
library(ggpubr)

ff<-list.files("/data_pool_1/PRADO/processed/dna/SAREK_pipeline/output/variant_calling/facets/results", recursive=TRUE, pattern = "*._hisens.rds", full.names=TRUE)

make_gistic_seg <- function(ssegs, dipLogR, sample_id) {
  ssegs %>%
    transmute(
      Sample = sample_id,
      Chromosome = chrom,
      Start_Position = start,
      End_Position = end,
      Num_Markers = num.mark,
      Seg_CN = cnlr.median - dipLogR
    )
}

res_ls <- list()
for(file in ff) {
    
    tid <- unlist(strsplit(file, "/"))[11]
    facetsObj <- readRDS(file)
    purity <- facetsObj$purity
    if(purity >= 0.3 & !is.na(purity)) {
        dipLogR <- facetsObj$dipLogR
        ssegs <- facetsObj$segs %>% filter(nhet >= 10, !is.na(mafR))
        res_ls[[tid]] <- make_gistic_seg(ssegs, dipLogR, tid)

    } else {
        print(tid)
    }
    
        

}


gistic_prado <- do.call(rbind, res_ls)

write.table(gistic_prado, "/data_pool_1/PRADO/processed/dna/SAREK_pipeline/output/variant_calling/GISTIC2/input/prado.seg", quote=FALSE, row.names=FALSE, sep="\t")
