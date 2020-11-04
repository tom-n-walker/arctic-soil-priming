
source("packages.R")
loadd()
abio <- abiotic
mtbs <- metabolites$chemistry
prim <- priming
fun <- fungi$diversity
bac <- bacteria$diversity

metsPeaks <- select(metabolites$peakAreas, -sampleID)
funPeaks <- select(fungi$relAbun, -sampleID)
bacPeaks <- select(bacteria$relAbun, -sampleID)

mtbs2fun <- match(metabolites$peakAreas$sampleID,
                  fungi$relAbun$sampleID)
metsPeaksFun <- metsPeaks[!is.na(mtbs2fun), ]
funPeaks <- funPeaks[mtbs2fun[!is.na(mtbs2fun)], ]

mtbs2bac <- match(metabolites$peakAreas$sampleID,
                  bacteria$relAbun$sampleID)
metsPeaksBac <- metsPeaks[!is.na(mtbs2bac), ]
bacPeaks <- bacPeaks[mtbs2bac[!is.na(mtbs2bac)], ]

mtbs2bac <- match(metabolites$peakAreas$sampleID,
                  bacteria$relAbun$sampleID)
metsPeaksBac <- metsPeaks[!is.na(mtbs2bac), ]
bacPeaks <- bacPeaks[mtbs2bac[!is.na(mtbs2bac)], ]

protest(X = vegdist(metsPeaksFun, "bray"),
        Y = vegdist(funPeaks, "bray"))
protest(X = vegdist(metsPeaksBac, "bray"),
        Y = vegdist(bacPeaks, "bray"))




adonis2(drivers ~ control, "euclidean")



all <- mtbs %>%
  left_join(., abio, by = "sampleID") %>%
  left_join(., prim, by = "sampleID") %>%
  left_join(., fun, by = "sampleID") %>%
  left_join(., bac, by = "sampleID") %>%
  filter(horizon != "O") %>%
  filter(!is.na(clay_perc)) %>%
  filter(!is.na(mtbHydrogen)) %>%
  filter(!is.na(richness.x)) %>%
  filter(!is.na(shannon.y)) %>%
  mutate(binary_horizon = as.numeric(as.factor(horizon))) %>%
  select(-sampleID, -horizon)

drivers <- all %>%
  select(-control, -cellulose, -protein, -Rcellulose, -Rprotein)
control <- all$Rprotein

index <- createDataPartition(control, p = 0.75, list = F)
trainDrivers <- drivers[index, ]
trainResponse <- control[index]

rf <- randomForest(x = trainDrivers,
                   y = trainResponse)
preds <- predict(rf, drivers[-index, ])
plot(preds ~ control[-index])

varImpPlot(rf)

ncol(all)
