pheno <- as.data.frame(TrainingPheno)
geno <- as.data.frame(TrainingGeno)
trainingdata <- cbind(pheno, geno)
colnames(trainingdata) <- paste("ID",1:ncol(trainingdata), sep="") ##1-605 because the SNP chip has 605 SNPs + phenotypes may have to change if you have a different num. SNPS###
##note ID1 will be the phenotype, IDs 2-606 are genotypes##


## create TRN TST split ###
train_index <- sample(1:nrow(trainingdata), 0.75 * nrow(trainingdata))
trainingset <- trainingdata[train_index, ]
testingset <- trainingdata[-train_index, ]

## create cross validation strategy ##
control <- trainControl(method='repeatedcv', 
                        number=5, ##will test 10 different values for mtry (number of variables for splitting) ##
                        repeats=1,
                        search = "random")  

trainMethod = "rf"
if (nCores > 1){
    print(paste("Creating cluster with", nCores-1, "cores..."))
    cl <- makeForkCluster(nCores - 1)
    registerDoParallel(cl)
    trainMethod = "parRF"
    print("Training in parallel...")
}

##build model##

rf_fit = train(ID1 ~ ., 
               data = trainingset, 
               method = trainMethod,
               tuneLength = 10,
               trControl=control) ## search a random tuning grid ##

if (nCores > 1){
    print("Training finished.")
    stopCluster(cl)
}

### This command takes about 90 minutes in an compute canada interactive session ###
