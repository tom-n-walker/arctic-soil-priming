################################################################################
#### Project: Arctic priming
#### Title:   Functions | Wrangling | Clean metabolite data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    15 July 2020
#### ---------------------------------------------------------------------------

## METS: REPLACE NA, RELATIVE ABUNDANCE, RENAME ##

clean_metabolites <- function(other_data, distance){
  ## Peaks house-keeping ----
  # select and subset mtbs peak areas
  peaks <- other_data$om_composition$peaks %>%
    select(Pyridine:Naphthalene..1.4.6.trimethyl.)
  # change missing to 0
  peaks[is.na(peaks)] <- 0
  # calculate peak relative abundances
  totals <- rep(
    rowSums(peaks), 
    each = ncol(peaks)
  ) %>%
    matrix(
      nrow = nrow(peaks), 
      byrow = T
    )
  relAbun <- peaks / totals
  
  ## Perform ordinations ----
  ordinations <- ordinate(
    nums = relAbun,
    distance = distance
  )
  
  ## Mtb info house-keeping ----
  # select and format
  mtbs <- other_data$om_composition$compounds %>%
    # format to make compound names match peaks column name
    mutate(column_name = make.names(column_name))
  # get smiles and add to data frame
  smiles <- cir_query(inchikey, match = "first")
  smiles[names(smiles) == "VRQGSFPALWCJBE-RQYWWVFRSA-N"] <- "CO[C@H]1[C@@H]([C@H]([C@H](O[C@@H]1C(=O)OC)OC)O)O"
  smiles[names(smiles) == "OCJFQVCIXIGACT-RBXMUDONSA-N"] <- "CC(=O)O[C@@H]1[C@H](CO[C@@H]([C@H]1O)OC(=O)C)O"
  smiles[names(smiles) == "JVBICOCAMPYKMA-UIAUGNHASA-N"] <- "CC(=O)O[C@H]1[C@H]2[C@H](O2)[C@H]3CO[C@@H]1O3"
  smilesDF <- data.frame(
    inchikey = names(smiles),
    smiles = smiles
  )
  mtbs <- left_join(mtbs, smilesDF)
  
  ## Derive sample-wise chemical information and info theory metrics ----
  chemInfo <- get_atom_counts(relAbun, mtbs)
  infoThry <- apply_info_theory(relAbun)
  
  ## Get molecular fingerprints ----
  # generate smiles vectors
  smilesVector <- mtbs$smiles
  names(smilesVector) <- mtbs$name
  # parse to correct atom format
  smilesParsed <- parse.smiles(smilesVector, omit.nulls = T)
  # generate molecular fingerprints
  fingerprints <- lapply(smilesParsed, get.fingerprint, type = "extended")
  # build similarities
  similarities <- fingerprint::fp.sim.matrix(fingerprints, method = "tanimoto")
  # convert to dissimilarities
  dissim <- 1 - similarities
  # add row and column names
  colnames(dissim) <- rownames(dissim) <- names(smilesParsed)
  # remove large redundant objects (occasional hang)
  rm(smilesParsed, fingerprints, similarities)
  # convert to distance object
  dissim <- as.dist(dissim)
  # cluster
  mtbClust <- hclust(dissim)
  
  ## Collate and return ----
  # create sample ID column
  sampleID <- other_data$om_composition$peaks %>%
    select(sampleID = Site_pit)
  # add sampleID to both rel abund and chem info data
  relAbunOut <- bind_cols(
    sampleID,
    relAbun
  )
  chemInfoOut <- bind_cols(
    sampleID,
    chemInfo
  )
  # add info theory metrics to chem info data
  chemInfoOut <- chemInfoOut %>%
    mutate(entropy = infoThry$obsEntropy,
           specificity = infoThry$obsSpecificity,
           nmds1 = ordinations$nmds$points[, 1],
           nmds2 = ordinations$nmds$points[, 2],
           pcoa1 = ordinations$pcoa[, 1],
           pcoa2 = ordinations$pcoa[, 2])
  # add compound specificity to mtbs dataset
  mtbsOut <- mtbs %>%
    select(inchikey, smiles, column_name, name, C:complexity) %>%
    mutate(specificity = infoThry$varSpecificity,
           nmds1loadings = ordinations$nmds$species[, 1],
           nmds2loadings = ordinations$nmds$species[, 2])
  # collate output
  out <- list(
    peakAreas = relAbunOut,
    nmds = ordinations$nmds,
    clustered_compounds = mtbsClust,
    chemistry = chemInfoOut,
    compounds = mtbsOut
  )
  return(out)
}
