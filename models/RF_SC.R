library(e1071)
library(devtools)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cluster)
library(factoextra)

if (args$trainingData == "F2")
    M = as.data.frame(pullSegSiteGeno(F2)
    y = as.data.frame(pheno(F2)))

if (args$trainingData == "F5")
    M = as.data.frame(pullSegSiteGenoF5)
    y = as.data.frame(pheno(F5))

if (args$trainingData == "F2_and_F5")
    F2M = as.data.frame(pullSegSiteGeno(F2))
    F2y = as.data.frame(pheno(F2))
    F5M = as.data.frame(pullSegSiteGeno(F5))
    F5y = as.data.frame(pheno(F5))

    M = rbind(F2M,F5M)
    y = rbind(F2y,F5y)

StratClusTRN(y,M) #calls function for stratified clustering algorithm

BV <- OptimPheno

trainingset= as.data.frame(cbind(BV,OptimGeno))
colnames(trainingset) <- paste("ID",1:(ncol(y) + ncol(M)), sep="")




## create cross validation strategy ##
control <- trainControl(method='repeatedcv', 
                        number=10, ##will test 10 different values for mtry (number of variables for splitting) ##
                        repeats=3,
                        search = "random")  

##build model##

rf_fit = train(ID1 ~ ., 
               data = trainingset, 
               method = "ranger",
               tuneLength = 10,
               trControl=control) ## search a random tuning grid ##

### This command takes about 90 minutes in an compute canada interactive session ###

## look at the parameters of the model ##
print(rf_fit) 
