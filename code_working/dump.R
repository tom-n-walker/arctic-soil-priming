rm(list = ls())

#### ATTEMPTING ENSEMBLE LEARNING ----------------------------------------------


#### PACKAGES AND OPTIONS ----

library(tidyverse)
library(caret)
set.seed(1)


#### DATA ----

# load data
collated <- drake::readd(collated)

# preprocess by imputing, centring, scaling
prepValues <- preProcess(collated, 
                         method = c("medianImpute", "center", "scale"))
noMissing <- predict(prepValues, collated)

# split data into training/test sets
index <- createDataPartition(noMissing$control, p = 0.75, list = F)
trainSet <- noMissing[index, ]
testSet <- noMissing[-index, ]


#### PERFORM ENSEMBLE MODEL ----------------------------------------------------

## Set up control, predictors, response
fitControl <- trainControl(method = "cv",
                           number = 5)
predictors <- c("pH", "soc_perc", "moc_detrend", "poc_detrend", "cmic_detrend",
                "fe_detrend", "om#1", "om#2", "fun#1", "fun#2", "bac#1", "bac#2")
response <- "control"


## Random forest ----

# run model
modelRF <- train(trainSet[, predictors], 
                 trainSet[, response], 
                 method = "rf", 
                 trControl = fitControl, 
                 tuneLength = 3)

get_best_result(modelRF)

# predict
testSet$predRF <- predict(modelRF, testSet[, predictors])

# visualise
par(mfrow = c(1, 2))
plot(testSet$control, testSet$predRF)
plot(modelRF)
