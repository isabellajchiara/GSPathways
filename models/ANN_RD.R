library(tidyverse)
library(keras)
library(readr)

## create GRM ##
geno <- as.matrix(M)
GM <- tcrossprod(geno)/dim(geno)
LG <- GM

Y <- as.matrix(y)
X = LG

X <- as.matrix(M)
Y <- as.matrix(y)

train_index <- sample(1:nrow(Y), 0.75 * nrow(Y))
X_train <- X[train_index, ]
Y_train <- Y[train_index,]

X_test <- X[-train_index, ]
Y_test <- Y[-train_index, ]

## define hyperparameters##

nEpoch = 100
batchSize = 100
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





