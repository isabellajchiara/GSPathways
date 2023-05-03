set.seed(123)
y <- as.data.frame(pheno(F5))
x <- as.data.frame(pullSegSiteGeno(F5))
Training = as.data.frame(cbind(y,x))
colnames(Training) <- paste("ID",1:(ncol(y) + ncol(x)), sep="")

train_index <- sample(1:nrow(Training), 0.9 * nrow(Training))
SR_train <- Training[train_index, ]
SR_test <- Training[-train_index, ]

##fit model, predict pheno on all markers
SVMfit = svm(ID1 ~ ., data = Training, kernel = "radial", cost=10, scale=FALSE)