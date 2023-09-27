library(tidyverse)
library(keras)
library(readr)


# add layers
predictions <- inputs %>% 
  layer_dense(units = ncol(X_train), activation = 'relu') %>% 
  layer_dropout(rate = 0.5) %>%
  layer_dense(units = 1)

# create and compile model 
model <- keras_model(inputs = inputs, outputs = predictions) 
model %>% compile( 
  optimizer = 'rmsprop', 
  loss = 'sparse_categorical_crossentropy',
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



