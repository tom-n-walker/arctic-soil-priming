################################################################################
#### Project: Arctic priming
#### Title:   Function | Wrangle | Load MCC data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    15 July 2020
#### ---------------------------------------------------------------------------

load_mcc <- function(){
  # setup options
  options(stringsAsFactors = F)
  # load directory
  dir <- "./data/raw_data/"
  # load maps
  funmap1 <- fread(paste0(dir, "fungi_mapping.txt"), data.table = F)
  bacmap1 <- fread(paste0(dir, "bac_mapping.txt"), data.table = F)
  bacmap2 <- fread(paste0(dir, "bac_mapping_second.txt"), data.table = F)
  bacmap3 <- fread(paste0(dir, "bac_mapping_third.txt"), data.table = F)
  # load data
  fundat1 <- fread(paste0(dir, "fungi_priority_otu.txt"), data.table = F, header = T)
  fundat2 <- fread(paste0(dir, "fungi_second_otu.txt"), data.table = F, header = T)
  fundat3 <- fread(paste0(dir, "fungi_third_otu.txt"), data.table = F, header = T)
  bacdat1 <- fread(paste0(dir, "bac_priority_otu.txt"), data.table = F, header = T)
  bacdat2 <- fread(paste0(dir, "bac_second_otu.txt"), data.table = F, header = T)
  bacdat3 <- fread(paste0(dir, "bac_third_otu.txt"), data.table = F, header = T)
  # collate data and maps separately
  bacmap <- bind_rows(bacmap1, bacmap2, bacmap3)
  bacdat <- bacdat1 %>% 
    left_join(., bacdat2, by = "#OTU ID") %>%
    left_join(., bacdat3, by = "#OTU ID")
  fundat <- fundat1 %>% 
    left_join(., fundat2, by = "#OTU ID") %>%
    left_join(., fundat3, by = "#OTU ID")
  # compile and export
  bac <- list(data = bacdat,
              map = bacmap)
  fun <- list(data = fundat,
              map = funmap1)
  out <- list(bac = bac, fun = fun)
  return(out)
}