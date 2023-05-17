library(devtools)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cluster)
library(factoextra)
##Build GS model using F2 as TRN to get EBVs##

set.seed(23489)
y <- as.data.frame(TrainingPheno)
x <- as.data.frame(TrainingGeno)
Training = as.data.frame(cbind(y,x))
colnames(Training) <- paste("ID",1:ncol(Training), sep="")


newgeno <- x %>%  select(where(~ n_distinct(.) > 1))

colnames(newgeno) =NULL

PCAgeno <- prcomp(newgeno, center=TRUE, scale=TRUE) ##take out categorical columns##

PCAselected = as.data.frame(-PCAgeno$x[,1:3])

silhouette <- fviz_nbclust(PCAselected, kmeans, method = 'silhouette')
silhouette$data ##largest value tells how many clusters are optimal ##

k=7
kmeans_geno = kmeans(PCAselected, centers = k, nstart = 50)
clusters <- fviz_cluster(kmeans_geno, data = PCAselected)

clusterData <- clusters$data

clusterData <- clusterData[order(clusterData$cluster),]

nclusters <- as.numeric(clusterData[as.numeric(nrow(clusterData)),as.numeric(ncol(clusterData))])

i = 1
datalist = vector("list", length = nclusters)
for (i in 1:nclusters) {
  clustername <- paste0("cluster",i)
  clustername <- clusterData[clusterData$cluster==i,] 
  
  assign(paste0("cluster",i), clustername)
  
  trnname <- paste0("trn",i)
  trnname <- clustername[sample(0.3*nrow(clustername)),]
  datalist[[i]] <- trnname

  i = i + 1
  if (i > nclusters){ 
    break
  }
  }
  
TRN <- do.call(rbind, datalist)


TRN <- TRN[,1]

M <- as.data.frame(x)   
rownames(M) <- c(1:nrow(M))
OptimGeno <- M[TRN,]
y <- as.data.frame(y)
OptimPheno <- as.data.frame(y[TRN,])

Optim_train <- cbind(OptimPheno, OptimGeno)
dim(Optim_train)
colnames(Optim_train) <- paste("ID",1:ncol(Optim_train), sep="")


##fit model, predict pheno on all markers
SVMfit = svm(ID1 ~ ., data = Optim_train, kernel = "radial", cost = 10, scale = FALSE)
