################################################################################
#### Project: Arctic priming
#### Title:   Function | Wrangling | Calculate Shannon entropy and specificity
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    14 July 2020
#### ---------------------------------------------------------------------------


apply_info_theory <- function(input){
  ## Pre-amble ----
  # packages
  require(entropy)
  # define function for specificity
  specificity <- function(x){
    # same as for entropy, but divide by N and non-negative sum
    out <- sum(ifelse(x > 0, x * log2(x), 0)) / length(x)
    return(out)
  }
  ## Calculate entropy ----
  obsEntropy <- apply(
    input,
    1, 
    entropy,
    unit = "log2"
  )
  
  ## Calculate specificity ----
  # calculate variable mean abundance
  varMean <- apply(
    input,
    2,
    mean
  )
  # calculate difference between abundance and mean abundance
  varDiffs <- apply(
    input,
    1,
    function(x){x/varMean}
  ) %>%
    t
  # calculate variable specificity
  varSpecificity <- apply(
    varDiffs,
    2, 
    specificity
  )
  # calculate sample specificity
  obsSpecificity <- apply(
    input,
    1,
    function(x){
      sum(x * varSpecificity)
    }
  )
  ## Return ----
  out <- list(obsEntropy = obsEntropy,
              varSpecificity = varSpecificity,
              obsSpecificity = obsSpecificity)
  return(out)
}
