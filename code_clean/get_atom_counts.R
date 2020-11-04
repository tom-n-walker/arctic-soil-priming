################################################################################
#### Project: Arctic priming
#### Title:   Function | Small function | Get atom counts for mtbs data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    16 July 2020
#### ---------------------------------------------------------------------------

get_atom_counts <- function(relAbun, mtbs){
  ## Per sample mean atomic info based on presence-absence ----
  chemInfo <- apply(relAbun, 1, function(x){
    # index of samples containing each compound
    present <- x > 0
    # filter columns for these samples
    out <- data.frame(
      mtbCarbon = mean(mtbs$C[present], na.rm = T),
      mtbHydrogen = mean(mtbs$H[present], na.rm = T),
      mtbNitrogen = mean(mtbs$N[present], na.rm = T),
      mtbOxygen = mean(mtbs$O[present], na.rm = T),
      mtbMolMass = mean(mtbs$mol_mass[present], na.rm = T),
      mtbComplex = mean(mtbs$complexity[present], na.rm = T)
    )
  }) %>%
    # bind into dataframe
    do.call(rbind, .)
  ## Per sample atom counts based on relative abundance ----
  cwmChemInfo <- apply(relAbun, 1, function(x){
    out <- data.frame(
      cwmCarbon = weighted.mean(mtbs$C, x, na.rm = T),
      cwmHydrogen = weighted.mean(mtbs$H, x, na.rm = T),
      cwmNitrogen = weighted.mean(mtbs$N, x, na.rm = T),
      cwmOxygen = weighted.mean(mtbs$O, x, na.rm = T),
      cwmMolMass = weighted.mean(mtbs$mol_mass, x, na.rm = T),
      cwmComplex = weighted.mean(mtbs$complexity, x, na.rm = T)
    )
  }) %>%
    # bind into dataframe
    do.call(rbind, .)
  ## Bind to output ----
  finished <- bind_cols(chemInfo, cwmChemInfo)
  return(finished)
}