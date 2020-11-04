# 
# 
# 
# ### ORDINATE FUNCTION TO USE ---------------------------------------------------
# 
# ordinate <- function(input){
#   # packages
#   require(vegan)
#   require(parallelDist)
#   # split nums and cats
#   nums <- numsPA <- input %>% select(-sampleID)
#   # make presence/absence
#   numsPA[numsPA > 0] <- 1
#   # generate bray curtis and cosine distance matrices
#   bray <- vegdist(nums)
#   brayPA <- vegdist(numsPA)
#   csne <- parDist(as.matrix(nums), method = "cosine")
#   csnePA <- parDist(as.matrix(numsPA), method = "cosine")
#   # perform nmds, pcoas
#   nmds <- metaMDS(nums, k = 2, tries = 200)
#   nmdsPA <- metaMDS(numsPA, k = 2, tries = 200)
#   pcoaBray <- cmdscale(bray, k = 2, eig = T)
#   pcoaCosine <- cmdscale(csne, k = 2, eig = T)
#   pcoaBrayPA <- cmdscale(brayPA, k = 2, eig = T)
#   pcoaCosinePA <- cmdscale(csnePA, k = 2, eig = T)
#   # collate scores
#   scores <- data.frame(sampleID = input$sampleID,
#                        brayPCOA1 = pcoaBray$points[, 1],
#                        brayPCOA2 = pcoaBray$points[, 2],
#                        brayBinPCOA1 = pcoaBrayPA$points[, 1],
#                        brayBinPCOA2 = pcoaBrayPA$points[, 2],
#                        csnePCOA1 = pcoaCosine$points[, 1],
#                        csnePCOA2 = pcoaCosine$points[, 2],
#                        csneBinPCOA1 = pcoaCosinePA$points[, 1],
#                        csneBinPCOA2 = pcoaCosinePA$points[, 2],
#                        NMDS1 = nmds$points[, 1],
#                        NMDS2 = nmds$points[, 2],
#                        NMDSBin1 = nmdsPA$points[, 1],
#                        NMDSBin2 = nmdsPA$points[, 2])
#   # collate ordinations
#   ords <- list(nmds,
#                nmdsBin = nmdsPA,
#                pcoaBray = pcoaBray,
#                pcoaBinBray = pcoaBrayPA,
#                pcoaCosine = pcoaCosine,
#                pcoaBinCosine = pcoaCosinePA)
#   # return
#   out <- list(scores = scores,
#               ordinations = ords)
#   return(out)
# }
# 
# 
# 
# 
# 
# ### COLLAPSE ALL ---------------------------------------------------------------
# 
# collate_all <- function(abiotic, priming, metabolites, fungi, bacteria){
#   # apply ordination and imputations to mtbs
#   mtb_fin <- ordinate(metabolites)
#   fun_fin <- ordinate(fungi)
#   bac_fin <- ordinate(bacteria)
#   # apply function to detrend abiotic and imputs
#   abi_fin <- detrend_abiotic(abiotic)
#   # impute priming
#   res_fin <- impute_all(priming)
#   # build categoricl skeleton
#   skeleton <- abiotic %>%
#     select(sampleID, horizon) %>%
#     mutate(site = substr(sampleID, 1, 2),
#            pit = paste0("P", sprintf("%02s", substr(sampleID, 5, 7)))) %>%
#     select(sampleID, site, pit, horizon)
#   # expand all data to it
#   expanded <- skeleton %>%
#     bind_cols(., abi_fin) %>%
#     left_join(., mtb_fin$nmds_scores, by = "sampleID") %>%
#     left_join(., fun_fin$nmds_scores, by = "sampleID") %>%
#     left_join(., bac_fin$nmds_scores, by = "sampleID") %>%
#     bind_cols(., res_fin) %>%
#     select(-sampleID)
#   # change column names
#   colnames(expanded) <- c(colnames(expanded)[1:17],
#                           "om#1", "om#2", "omRic", "omShn", "omSmp",
#                           "fun#1", "fun#2", "funRic", "funShn", "funSmp", 
#                           "bac#1", "bac#2", "bacRic", "bacShn", "bacSmp",
#                           colnames(expanded)[33:37])
#   # get nmds outputs
#   all_nmds <- list(om_composition = mtb_fin$nmds_metadata,
#                    fungi = fun_fin$nmds_metadata,
#                    bacteria = bac_fin$nmds_metadata)
#   # get original data
#   original <- list(categorical = skeleton,
#                    physical = abi_fin,
#                    om_composition = mtb_fin$values,
#                    fungi = fun_fin$values,
#                    bacteria = bac_fin$values,
#                    response = res_fin)
#   # combine to list
#   out <- list(collated = expanded,
#               separate = original)
#   return(out)
# }
# 
# 
# ### CONVENIENCE WRAPPERS -------------------------------------------------------
# 
# split_collated <- function(finished){
#   out <- finished$collated
#   return(out)
# }
# 
# split_ordinations <- function(finished){
#   out <- finished$ordinations
#   return(out)
# }
# 
# split_original <- function(finished){
#   out <- finished$original
#   return(out)
# }
# 
# 
# 
# 
