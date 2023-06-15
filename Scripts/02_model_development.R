## Model generation: Random Forest and empirical (NDCI) algorithm for (Chl-a) 

# loading require packages

require(data.table)
require(dplyr)
require(Metrics)
require(randomForest)

source("Scripts/Functions.R")
## Loading data

data = fread("Outputs/sentinel2_simulated_filtered.csv")



## Two models: Chl-a Random Forest and NDCI empirical


# Steps:
# 1) Separate data into training / test 
# 2) Calibrate the model with training data
# 3) Validate the model with test data
# 4) Create the full model based on all data
# 5) Apply to satellite image
# 6) If available, validate the image predictions


## Creating the validation and training datasets (70% train / 30% validation)

set.seed(13) # To allow replicability

samples = sample(x = 1:nrow(data),
                 size = 0.7*nrow(data), 
                 replace = F)


train = data[samples,]
valid = data[-samples,]


dim(train)
dim(valid)


## Creating an empirical (simple model) based on Normalized Difference Chlorophyll Index


emp.chla = lm(log(Chla)~log(NDCI), data = train)

summary(emp.chla)


valid$NDCI_CHLA = exp(predict(emp.chla, valid))

plot(valid$Chla, valid$NDCI_CHLA , xlim = c(0,400), ylim = c(0,400))
abline(0,1)


###### Create a random forest algorithm for Chl-a


chla.rf = randomForest(Chla~ x490+x560+x660+x705+NDCI, data = train, 
                       ntree = 200, mtry = 4, importance = T)

varImpPlot(chla.rf)

valid$RF_CHLA = (predict(chla.rf, valid))

plot(valid$Chla, valid$RF_CHLA ,pch = 20, xlab = "Measured Chla",
     ylab = "Predicted Chla", xlim = c(0,400), ylim = c(0,400))

abline(0,1)

# Lets compare both models

stat_calc(real = valid$Chla, predicted = valid$NDCI_CHLA)

stat_calc(real = valid$Chla, predicted = valid$RF_CHLA)

## Random Forest presented better results. 


## Calculate empirical and RF final models

RF_FINAL = randomForest(Chla~ x490+x560+x660+x705+NDCI, data = data, 
             ntree = 200, mtry = 4, importance = T)

EMP_FINAL = lm(log(Chla)~log(NDCI), data = data)

saveRDS(RF_FINAL, file = 'Outputs/rf_chla.R')
saveRDS(EMP_FINAL, file = 'Outputs/emp_chla.R')






