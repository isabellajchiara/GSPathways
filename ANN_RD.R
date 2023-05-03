library(tidyverse)
library(keras)
library(readr)
library(BMTME)

#BUILD ANN using PYT TP ##

## create GRM ##

X <- as.matrix(TrainingGeno)
Y <- as.matrix(TrainingPheno)

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
  X, Y,
  shuffle=F,
  epochs = No_Epoch, batch_size = 10,
  validation_split = 0.2,
  verbose = 0)

## Refit model using early stopping ##

early_stop <- callback_early_stopping(monitor="val_loss", mode='min', patience =50)

model_Final <- build_model()
model_fit_Final <- model_Final%>%fit(
  X, Y,
  shuffle=F,
  epochs = No_Epoch, batch_size=10,
  validation_split = 0.2,
  verbose=0, callbacks = list(early_stop)
)

###### Model is built, continue with breeding simulation #####

## Use model to predict ##

# prediction = model_Final %>% predict(genoMatrix)


