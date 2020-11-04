################################################################################
#### Project: Arctic priming
#### Title:   Function | Wrangle | Load physical, OM and response data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    15 July 2020
#### ---------------------------------------------------------------------------

load_other_data <- function(){
  # no factors
  options(stringsAsFactors = F)
  # directory
  dir <- "./data/raw_data/"
  # load metabolite information
  mtbs <- list(
    peaks = read.csv(paste0(dir, "area_peaks_named.csv"), 
                     sep = ";", 
                     dec = ","),
    compounds = read.csv(paste0(dir, "mtbs_info.csv"))
  )
  # put together into list
  out <- list(priming = read.csv(paste0(dir, "respiration.csv")),
              om_composition = mtbs,
              minerals = read.csv(paste0(dir, "minerals.csv")),
              physical = read.csv(paste0(dir, "other_physical.csv")))
  # return list
  return(out)
}
