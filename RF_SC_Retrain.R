library(e1071)
library(devtools)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cluster)
library(factoextra)

library(caret)
library(ranger)
library(randomForest)

TrainingPheno2 <- pheno(F5)
TrainingGeno2 <- pullSegSiteGeno(F5)

y <- as.data.frame(TrainingPheno2)
M <- as.data.frame(TrainingGeno2)

newgeno <- M %>%  select(where(~ n_distinct(.) > 1))

colnames(newgeno) =NULL

PCAgeno <- prcomp(newgeno, center=TRUE, scale=TRUE) ##take out categorical columns##

PCAselected = as.data.frame(-PCAgeno$x[,1:3])

silhouette <- fviz_nbclust(PCAselected, kmeans, method = 'silhouette')
kvalues <- silhouette$data ##largest value tells how many clusters are optimal ##
kvalues <- kvalues[order(-kvalues$y),]

k=as.numeric(kvalues[1,1])

kmeans_geno = kmeans(PCAselected, centers = k, nstart = 50)
clusters <- fviz_cluster(kmeans_geno, data = PCAselected)

clusterData <- clusters$data

clusterData <- clusterData[order(clusterData$cluster),]

nclusters = nrow(clusterData)

for (i in nclusters) {
  clustername <- paste0("cluster",i)
  clustername <- clusterData[clusterData$cluster=i] 
  
  trnname <- paste0("trn",i)
  trnname <- clustername[sample(0.3*nrow(clustername)),]
  
  i = i + 1
  if (i > ncluster){ 
    break
  }
  
TRN <- rbind(trn1, trn2,trn3,trn4,trn5, trn6, trn7,trn8,trn9,trn10)
TRN <- na.omit(TRN)

TRN <- TRN[,1]

OptimGeno <- M[TRN,]
y <- as.data.frame(y)
OptimPheno <- y[TRN,]

BV <- OptimPheno

trainingset= as.data.frame(cbind(BV,OptimGeno))
colnames(trainingset) <- paste("ID",1:(ncol(y) + ncol(M)), sep="")

## create cross validation strategy ##
control <- trainControl(method='repeatedcv', 
                        number=10, ##will test 10 different values for mtry (number of variables for splitting) ##
                        repeats=3,
                        search = "random")  

##build model##

rf_fit2 = train(ID1 ~ ., 
                data = trainingset, 
                method = "ranger",
                tuneLength = 10,
                trControl=control) ## search a random tuning grid ##

### This command takes about 90 minutes in an compute canada interactive session ###

## look at the parameters of the model ##
print(rf_fit2) 
