library(AlphaSimR)
library(readxl)
library(writexl)
library(rrBLUP)
library(keras)
library(tensorflow)
library(readr)
library(BMTME)

#BUILD ANN using F1 as TP ##

ReTrainingPheno5 <- pheno(F5)
ReTrainingGeno5 <- pullSegSiteGeno(F5)

ReTrainingPheno4 <- pheno(F4)
ReTrainingGeno4 <- pullSegSiteGeno(F4)


ReTrainingPheno3 <- pheno(F3)
ReTrainingGeno3 <- pullSegSiteGeno(F3)

ReTrainingPheno2 <- pheno(F2)
ReTrainingGeno2 <- pullSegSiteGeno(F2)



##Create Matrices

Y_ReTrain <- as.matrix(rbind(TrainingPheno, ReTrainingPheno5,ReTrainingPheno4,ReTrainingPheno3,ReTrainingPheno2))
X_ReTrain = as.matrix(rbind(TrainingGeno, ReTrainingGeno5,ReTrainingGeno4,ReTrainingGeno3,ReTrainingGeno2))

## Build model using pipe operator ##
N_Units = ncol(X)
No_Epoch = 100


build_model <- function() {
  model <- keras_model_sequential()
  model %>%
    layer_dense(units =N_Units, activation = "relu", input_shape = c(dim
                                                                     (X)[2])) %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(units = 1, activation = "linear")
  model %>% compile(
    loss = "mse",
    optimizer = "rmsprop",
    metrics = c("mse"))
  model}

## build and view model ##

model <- build_model()
model %>% summary()

## inner training and validation ##

model_fit <- model %>% fit(
  X_ReTrain, Y_ReTrain,
  shuffle=F,
  epochs = No_Epoch, batch_size = 10,
  validation_split = 0.2,
  verbose = 0)

## Refit model using early stopping ##

early_stop <- callback_early_stopping(monitor="val_loss", mode='min', patience =50)

model_Final2 <- build_model()
model_fit_Final <- model_Final%>%fit(
  X, Y,
  shuffle=F,
  epochs = No_Epoch, batch_size=10,
  validation_split = 0.2,
  verbose=0, callbacks = list(early_stop)
)

###### Model is built, continue with breeding simulation #####

## Use model to predict ##

# prediction = model_Final2 %>% predict(genoMatrix)



