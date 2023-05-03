library(AlphaSimR)
library(readxl)
library(writexl)
library(rrBLUP)
library(keras)
library(tensorflow)
library(readr)
library(BMTME)


## create GRM ##
geno <- as.matrix(TrainingGeno)
GM <- tcrossprod(geno)/dim(geno)
LG <- cholesky(GM)

Y <- as.matrix(TrainingPheno)
X = LG

## create data frame with geno and pheno data ##

phenotypes <- as.data.frame(Y)
genotypes <- as.data.frame(geno)
data <-cbind(phenotypes,genotypes)

## stratified clustering training testing partition##

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

cluster1 <- clusterData[clusterData$cluster==1,]
cluster2 <- clusterData[clusterData$cluster==2,]
cluster3 <- clusterData[clusterData$cluster==3,]
cluster4 <- clusterData[clusterData$cluster==4,]
cluster5 <- clusterData[clusterData$cluster==5,]
cluster6 <- clusterData[clusterData$cluster==6,]
cluster7 <- clusterData[clusterData$cluster==7,]
cluster8 <- clusterData[clusterData$cluster==8,]
cluster9 <- clusterData[clusterData$cluster==9,]
cluster10 <- clusterData[clusterData$cluster==10,]

trn1 <- cluster1[sample(0.75*nrow(cluster1)),]
trn2 <- cluster2[sample(0.75*nrow(cluster2)),]
trn3 <- cluster3[sample(0.75*nrow(cluster3)),]
trn4 <- cluster4[sample(0.75*nrow(cluster4)),]
trn5 <- cluster5[sample(0.75*nrow(cluster5)),]
trn6 <- cluster6[sample(0.75*nrow(cluster6)),]
trn7 <- cluster7[sample(0.75*nrow(cluster7)),]
trn8 <- cluster8[sample(0.75*nrow(cluster8)),]
trn9 <- cluster9[sample(0.75*nrow(cluster9)),]
trn10 <- cluster9[sample(0.75*nrow(cluster10)),]
TRN <- rbind(trn1, trn2,trn3,trn4,trn5, trn6, trn7,trn8,trn9,trn10)
TRN <- na.omit(TRN)

TRN <- TRN[,1]

OptimGeno <- M[TRN,]
y <- as.data.frame(y)
OptimPheno <- y[TRN,]

BV <- OptimPheno

trainingset= as.data.frame(cbind(BV,OptimGeno))
colnames(trainingset) <- paste("ID",1:(ncol(y) + ncol(M)), sep="")


## one holdout CV ##

No_Epoch = 1000
N_Units = 33
X_trn = M[TRN,]
X_tst = M[-TRN,]
Y_trn = Y[TRN,]
Y_tst = Y[-TRN,]

## Build model using pipe operator ##

build_model <- function() {
  model <- keras_model_sequential()
  model %>%
    layer_dense(units =N_Units, activation = "relu", input_shape = c(dim
                                                                     (X_trn)[2])) %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 1, activation = "linear")
  model %>% compile(
    loss = "mse",
    optimizer = "rmsprop",
    metrics = c("mse"))
  model}

## build and view model ##

model <- build_model()
model %>% summary()

## fitting and plotting the model ##

print_dot_callback <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 20 == 0) cat("\n")
    cat(".")
  })

## inner training and validation ##

model_fit <- model %>% fit(
  X_trn, Y_trn,
  shuffle=F,
  epochs = No_Epoch, batch_size = 640,
  validation_split = 0.2,
  verbose = 0,
  callbacks = list(print_dot_callback))

## plot ## 

plot(model_fit)

## Refit model using early stopping ##

early_stop <- callback_early_stopping(monitor="val_loss", mode='min', patience =50)

model_Final <- build_model()
model_fit_Final <- model_Final%>%fit(
  X_trn, Y_trn,
  shuffle=F,
  epochs = No_Epoch, batch_size=640,
  validation_split = 0.2,
  verbose=0, callbacks = list(early_stop, print_dot_callback)
)

## plot history of training process ##

length(model_fit_Final$metrics$mean_squared_error)
plot(model_fit_Final)

## predition to see how model performs##

prediction = model_Final %>% predict(X_tst)
Predicted = c(prediction)
Observed = Y_tst
plot(Observed, Predicted)
MSE = mean((Observed-Predicted)^2)
MSE
Obs_Pred = cbind(Observed, Predicted)
colnames(Obs_Pred) = c("Observed", "Predicted")
Obs_Pred

###### Model is built, continue with breeding simulation #####

## Use model to predict F2 EBV ##

geno <- pullSnpGeno(F2)
geno <- as.matrix(geno)
GM <- tcrossprod(geno)/dim(geno)
LGF2 <- cholesky(GM)