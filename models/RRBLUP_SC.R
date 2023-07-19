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

EBVans <-mixed.solve(BV, Z=OptimGeno, K=NULL, X=NULL, SE=FALSE, return.Hinv=FALSE)

markerEffects <- as.matrix(EBVans$u)
markerEffects <- as.vector(markerEffects)
