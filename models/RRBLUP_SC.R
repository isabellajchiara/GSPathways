library(devtools)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cluster)
library(factoextra)


StratClusTRN(y,M) #calls function for stratified clustering algorithm

BV <- OptimPheno

EBVans <-mixed.solve(BV, Z=OptimGeno, K=NULL, X=NULL, SE=FALSE, return.Hinv=FALSE)

markerEffects <- as.matrix(EBVans$u)
markerEffects <- as.vector(markerEffects)
