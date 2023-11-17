library(tidyverse)
library(keras)
library(readr)

ar <- array(M, c(nInd,3547,2)
            
Y <- as.matrix(y)

## define hyperparameters##

nEpoch = 100
batchSize = 100
valSplit = 0.2

# add layers

inputs = layer_input(shape=c(dim(ar)))

If we create an array of dimension (2, 3, 4) then it creates 4 rectangular matrices 
each with 2 rows and 3 columns. 

predictions <- inputs %>% 
  layer_flatten(ar) %>% 
  layer_dense(units = ncol(ar), activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = ncol(ar), activation = 'relu') %>% 
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



