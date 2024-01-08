library(tidyverse)
library(keras)
library(tensorflow)
library(readr)

## create GRM ##
geno <- as.matrix(M)
GM <- tcrossprod(geno)/dim(geno)
LG <- GM

Y <- as.matrix(y)
X = LG

X_train <- as.matrix(M)
Y_train <- as.matrix(y)


## define hyperparameters##

cutoff =200

nEpoch = 100

if (nrow(geno) < cutoff){
  batchSize = 20
}else{
batchSize=500
}

valSplit = 0.2


# add layers

inputs = layer_input(shape=(ncol(X_train))) 

predictions <- inputs %>% 
  layer_dense(units = ncol(X_train), activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = ncol(X_train), activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 1)

# create and compile model 

model <- keras_model(inputs = inputs, outputs = predictions) 

earlyStopping <- callback_early_stopping(monitor = "mean_squared_error", min_delta = 0.1, patience = 10,
                                         mode = "auto")

model %>% compile( 
  optimizer = 'Adam', 
  loss = 'mean_squared_error',
  metrics = c('mean_squared_error')
) 

fit(
  model,
  x = X_train,
  y = Y_train,
  batch_size = batchSize,
  epochs = nEpoch,
  verbose = 0,
  validation_split = valSplit,
  callbacks = list(earlyStopping)
)

cli_alert_success("Fit Neural Net at {args$trainGen} using {args$trainingData} data")









