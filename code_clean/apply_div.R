################################################################################
#### Project: Arctic priming
#### Title:   Function | Wrangle | Calculate a range of diversity indices
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    16 July 2020
#### ---------------------------------------------------------------------------

apply_div <- function(input){
  out <- data.frame(
    richness = rowSums(input > 0),
    shannon = diversity(input),
    simpson = diversity(input, "simpson")
  )
  return(out)
}
