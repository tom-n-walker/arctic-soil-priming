################################################################################
#### Project: Arctic priming
#### Title:   Function | Wrangle | Clean priming response data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    15 July 2020
#### ---------------------------------------------------------------------------

clean_priming <- function(other_data){
  ## Clean data ----
  # rename and select
  out <- other_data$priming %>%
    # divide to get good units
    mutate(sampleID = paste(site, code, sep = "-"),
           control = control/1000,
           cellulose = cellulose/1000,
           protein = protein/1000) %>%
    # select variables of interest and rearrange
    select(sampleID, site, horizon, control, cellulose, protein,
           Rcellulose, Rprotein) %>%
    arrange(sampleID)
  ## Impute missing values ---- 
  # isolate numeric
  nums <- out %>% 
    select(control:Rprotein)
  # impute and bind column names
  imputed <- apply_mice(
    nums, 
    5
  ) %>%
    bind_cols(
      select(out, sampleID),
      .
    ) %>%
    arrange(sampleID)
  # return
  return(imputed)
}
