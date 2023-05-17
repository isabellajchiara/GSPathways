set.seed(123)
y <- TrainingPheno
x <- TrainingGeno
Training = as.data.frame(cbind(y,x))
colnames(Training) <- paste("ID",1:(ncol(y) + ncol(x)), sep="")

train_index <- sample(1:nrow(Training), 0.75 * nrow(Training))
train <- Training[train_index, ]
test <- Training[-train_index, ]

##fit model, predict pheno on all markers
SVMfit = svm(ID1 ~ ., data = train, kernel = "radial", cost=10, scale=FALSE)

