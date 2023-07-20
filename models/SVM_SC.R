library(devtools)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cluster)
library(factoextra)
##Build GS model using F2 as TRN to get EBVs##

set.seed(23489)

Training = as.data.frame(cbind(y,M))

colnames(Training) <- paste("ID",1:ncol(Training), sep="")

StratClusTRN(y,M) #calls function for stratified clustering algorithm

Optim_train <- cbind(OptimPheno, OptimGeno)
dim(Optim_train)
colnames(Optim_train) <- paste("ID",1:ncol(Optim_train), sep="")


##fit model, predict pheno on all markers
SVMfit = svm(ID1 ~ ., data = Optim_train, kernel = "radial", cost = 10, scale = FALSE)
