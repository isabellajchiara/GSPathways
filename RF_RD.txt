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
                        number=10, ##will test 10 different values for mtry (number of variables for splitting) ##
                        repeats=3,
                        search = "random")  

##build model##

rf_fit = train(ID1 ~ ., 
               data = trainingset, 
               method = "rf",
               tuneLength = 10,
               trControl=control) ## search a random tuning grid ##

### This command takes about 90 minutes in an compute canada interactive session ###
