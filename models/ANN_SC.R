library(AlphaSimR)
library(readxl)
library(writexl)
library(rrBLUP)
library(keras)
library(tensorflow)
library(readr)
library(factoextra)


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

Y <- as.data.frame(OptimPheno)
X <- as.data.frame(OptimGeno)

train_index <- sample(1:nrow(Y), 0.75 * nrow(Y))
X_train <- as.matrix(X[train_index, ])
Y_train <- as.matrix(Y[train_index,])

X_test <- as.matrix(X[-train_index, ])
Y_test <- as.matrix(Y[-train_index, ])


## one holdout CV ##

nEpoch = 100
batchSize = 100
valSplit = 0.2


# add layers

inputs = layer_input(shape=(ncol(X_train))) 

predictions <- inputs %>% 
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



