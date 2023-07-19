set.seed(123)
if (args$trainingData == "F2")
    M = as.data.frame(pullSegSiteGeno(F2)
    y = as.data.frame(pheno(F2)))

if (args$trainingData == "F5")
    M = as.data.frame(pullSegSiteGenoF5)
    y = as.data.frame(pheno(F5))

if (args$trainingData == "F2_and_F5")
    F2M = as.data.frame(pullSegSiteGeno(F2))
    F2y = as.data.frame(pheno(F2))
    F5M = as.data.frame(pullSegSiteGeno(F5))
    F5y = as.data.frame(pheno(F5))

    M = rbind(F2M,F5M)
    y = rbind(F2y,F5y)

Training = as.data.frame(cbind(y,M))
colnames(Training) <- paste("ID",1:(ncol(y) + ncol(M)), sep="")

train_index <- sample(1:nrow(Training), 0.75 * nrow(Training))
train <- Training[train_index, ]
test <- Training[-train_index, ]

##fit model, predict pheno on all markers
SVMfit = svm(ID1 ~ ., data = train, kernel = "radial", cost=10, scale=FALSE)

