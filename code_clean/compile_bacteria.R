################################################################################
#### Project: Arctic priming
#### Title:   Function | Wrangle | Compile bacteria community data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    15 July 2020
#### ---------------------------------------------------------------------------

compile_bacteria <- function(raw_mcc, distance){
  # isolate and clean mapping
  map <- raw_mcc$bac$map
  map <- map[!duplicated(map$Description), ]
  # transpose and clean abundance data
  dat <- raw_mcc$bac$data %>%
    select(-"#OTU ID") %>%
    t %>%
    as.data.frame %>%
    mutate("#SampleID" = rownames(.))
  # collate and filter abundance data
  ready <- map %>%
    right_join(., dat, by = "#SampleID") %>%
    filter(temp == 15 & treatment == "control") %>%
    select(insitu_code, V1:V1381)
  # change column names
  columns <- raw_mcc$bac$data$`#OTU ID`
  colnames(ready) <- c("sampleID", columns)
  # remove OTUs that are totally absent
  nums <- ready[, -1]
  nums <- nums[, colSums(nums) > 0]
  
  ## Perform ordinations ----
  # do ordinations
  ordinations <- ordinate(
    nums = nums,
    distance = distance
  )
  
  ## Derive diversity indices ----
  # information theory calculations
  infoThry <- apply_info_theory(nums)
  # classical diversity indices
  div <- apply_div(nums) %>%
    # add info theory metrics and ordination scores to div output
    mutate(entropy = infoThry$obsEntropy,
           specificity = infoThry$obsSpecificity,
           nmds1 = ordinations$nmds$points[, 1],
           nmds2 = ordinations$nmds$points[, 2],
           pcoa1 = ordinations$pcoa[, 1],
           pcoa2 = ordinations$pcoa[, 2])
  
  ## Generate taxonomy and add specificity ----
  # create taxonomy
  taxo <- substr(colnames(nums), 2, 1000) %>%
    strsplit(., ";") %>%
    do.call(rbind, .) %>%
    apply(., 2, substr, 4, 100) %>%
    as.data.frame %>%
    # add specificity from info theory and ordination loadings metrics
    mutate(specificity = as.numeric(infoThry$varSpecificity),
           nmds1loadings = ordinations$nmds$species[, 1],
           nmds2loadings = ordinations$nmds$species[, 2])
  # rename columns
  colnames(taxo) <- c("kingdom", "phylum", "clade", "order", 
                      "family", "genus", "specificity", 
                      "nmds1loadings", "nmds2loadings")
  
  ## Collate and return ----
  # create data frame of sampleID
  sampleID <- ready %>% select(sampleID)
  # bind columns together
  numsOut <- bind_cols(
    sampleID,
    nums
  )
  infoOut <- bind_cols(
    sampleID,
    div
  )
  # collate/return
  out <- list(
    relAbun = numsOut,
    nmds = ordinations$nmds,
    diversity = infoOut,
    taxonomy = taxo
  )
  return(out)
}