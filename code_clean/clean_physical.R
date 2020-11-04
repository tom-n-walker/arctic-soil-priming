################################################################################
#### Project: Arctic priming
#### Title:   Function | Wrangle | Clean physical data
#### Author:  Tom Walker (thomas.walker@usys.ethz.ch)
#### Date:    15 July 2020
#### ---------------------------------------------------------------------------

clean_physical <- function(other_data){
  ## Basic house-keeping ----
  # subset physicaldata
  phys <- other_data$physical %>%
    filter(Temp == 15) %>%
    select(-Site, -Code, -Temp)
  # summarise where more than one value available
  phys <- phys %>%
    group_by(site.code) %>%
    summarise(Cmic = mean(Cmic, na.rm = T),
              pH = mean(pH, na.rm = T), 
              CN_ratio = mean(CN_ratio, na.rm = T),
              Fep = mean(Fep, na.rm = T)) %>%
    ungroup %>%
    as.data.frame(stringsAsFactors = F)
  # subset mineral data
  mins <- other_data$minerals %>%
    mutate(site.code = paste(site, code, sep = "-")) %>%
    select(-site, -code)
  # calculate percentages, select columns and rearrange
  abiotic <- left_join(mins, phys, by = "site.code") %>%
    mutate(soc_perc = SOM_OC/10,
           poc_perc = POM_OC/10,
           moc_perc = MOM_OC/10,
           cmic_perc = Cmic/1000/10,
           fe_perc = Fep/10) %>%
    select(sampleID = site.code, horizon, 
           pH, clay_perc, silt_perc, sand_perc, 
           soil_cn = CN_ratio,
           soc_perc, poc_perc, moc_perc, cmic_perc, fe_perc) %>%
    arrange(sampleID)
  ## Impute missing values ----
  # subset for numbers
  nums <- abiotic %>% select(-sampleID, -horizon)
  imputed <- apply_mice(nums, 5)
  # remove values from O horizon (these are never there, so wrong to have them)
  o.fails <- c("clay_perc", "silt_perc", "sand_perc", 
               "poc_perc", "moc_perc", "fe_perc")
  imputed[which(abiotic$horizon == "O"), o.fails] <- NA
  ## Detrend carbon pools for soil C effect ----
  # detrend with linear models
  m_cmic <- lm(cmic_perc ~ poly(soc_perc, 2), imputed, na.action = "na.exclude")
  m_poc <- lm(poc_perc ~ poly(soc_perc, 2), imputed, na.action = "na.exclude")
  m_moc <- lm(moc_perc ~ poly(soc_perc, 2), imputed, na.action = "na.exclude")
  m_fe <- lm(fe_perc ~ poly(soc_perc, 2), imputed, na.action = "na.exclude")
  # extract outputs
  imputed$cmic_detrend <- residuals(m_cmic)
  imputed$poc_detrend <- residuals(m_poc)
  imputed$moc_detrend <- residuals(m_moc)
  imputed$fe_detrend <- residuals(m_fe)
  ## Return output ----
  out <- bind_cols(
    select(abiotic, sampleID, horizon),
    imputed
  )
  # return
  return(out)
}