library(AlphaSimR)
library(readxl)
library(writexl)
library(rrBLUP)
library(keras)
library(tensorflow)
library(readr)


## create GRM ##
geno <- as.matrix(M)
GM <- tcrossprod(geno)/dim(geno)
LG <- GM

Y <- as.matrix(y)
X = LG

## create data frame with geno and pheno data ##

phenotypes <- as.data.frame(Y)
genotypes <- as.data.frame(geno)
data <-cbind(phenotypes,genotypes)

## stratified clustering training testing partition##

y <- as.data.frame(y)
M <- as.data.frame(M)

StratClusTRN(y,M) #calls function for stratified clustering algorithm

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
