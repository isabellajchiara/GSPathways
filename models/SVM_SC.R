library(devtools)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cluster)
library(factoextra)
##Build GS model using F2 as TRN to get EBVs##

set.seed(23489)

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

Training = as.data.frame(cbind(y,M))

colnames(Training) <- paste("ID",1:ncol(Training), sep="")

StratClusTRN(y,M) #calls function for stratified clustering algorithm

Optim_train <- cbind(OptimPheno, OptimGeno)
dim(Optim_train)
colnames(Optim_train) <- paste("ID",1:ncol(Optim_train), sep="")


##fit model, predict pheno on all markers
SVMfit = svm(ID1 ~ ., data = Optim_train, kernel = "radial", cost = 10, scale = FALSE)
