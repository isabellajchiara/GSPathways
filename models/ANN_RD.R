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


## define input and number of epochs ##
inputs = layer_input(shape=(ncol(X_train))) 
nEpoch = 100
batchSize = 100
valSplit = 0.2

# add layers
predictions <- inputs %>% 
  layer_dense(units = ncol(X_train), activation = 'relu') %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1)

# create and compile model 
model <- keras_model(inputs = inputs, outputs = predictions) 
model %>% compile( 
  optimizer = 'rmsprop', 
  loss = 'binary_crossentropy',
  metrics = c('accuracy')
) 

fit(
  model,
  x = X_train,
  y = Y_train,
  batch_size = batchSize,
  epochs = nEpoch,
  verbose = 0,
  validation_split = valSplit 
)



