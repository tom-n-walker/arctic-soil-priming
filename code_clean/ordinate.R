################################################################################
#### Project: Arctic priming
#### Title:   Function | Wrangle | Ordinate
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    27 July 2020
#### ---------------------------------------------------------------------------

ordinate <- function(nums, distance){
  # require parallel dist package
  require(parallelDist)
  # extract data, do distance matrix
  dataDist <- nums %>%
    as.matrix %>%
    parDist(method = distance)
  # do NMDS
  nmds <- metaMDS(
    comm = nums,
    distance = "bray",
    k = 2,
    try = 200,
    wascores = T
  )
  # do PCoA
  pcoa <- cmdscale(
    d = dataDist,
    k = 2
  )
  # collate and return
  output <- list(
    nmds = nmds,
    pcoa = pcoa
  )
  return(output)
}