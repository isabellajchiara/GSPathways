pheno <- as.data.frame(TrainingPheno)
geno <- as.data.frame(TrainingGeno)
trainingdata <- cbind(pheno, geno)
colnames(trainingdata) <- paste("ID",1:ncol(trainingdata), sep="") ##1-605 because the SNP chip has 605 SNPs + phenotypes may have to change if you have a different num. SNPS###
##note ID1 will be the phenotype, IDs 2-606 are genotypes##


## create TRN TST split ###
train_index <- sample(1:nrow(trainingdata), min(1000, 0.75 * nrow(trainingdata)))
trainingset <- trainingdata[train_index, ]
testingset <- trainingdata[-train_index, ]

## create cross validation strategy ##
control <- trainControl(method='repeatedcv', 
                        number=10, ##will test 10 different values for mtry (number of variables for splitting) ##
                        repeats=3,
                        search = "random")  

trainMethod <- "rf"
if (args$nCores > 1){
    cli_alert_info("Creating cluster with {args$nCores} cores...")
    cl <- makePSOCKcluster(args$nCores)
    registerDoParallel(cl)
    trainMethod <- "parRF"
}

##build model##
cat ("Training...\n")

rf_fit = train(ID1 ~ ., 
               data = trainingset, 
               method = trainMethod,
               tuneLength = 10,
               trControl=control) ## search a random tuning grid ##

# deactivate cluster
if (args$nCores > 1)
    stopCluster(cl)

### This command takes about 90 minutes in an compute canada interactive session ###
