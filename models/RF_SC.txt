library(e1071)
library(devtools)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(cluster)
library(factoextra)

y <- as.data.frame(TrainingPheno)
M <- as.data.frame(TrainingGeno)

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

rf_fit = train(ID1 ~ ., 
               data = trainingset, 
               method = "ranger",
               tuneLength = 10,
               trControl=control) ## search a random tuning grid ##

### This command takes about 90 minutes in an compute canada interactive session ###

## look at the parameters of the model ##
print(rf_fit) 
