# Defines some variables and functions based on what model was chosen

# ADD NEW MODELS HERE
if (args$model == "rf with randomly sampled training set"){
    fileTrain <- "RF_RD.R"
    modelLibs <- c("caret","ranger","tidyverse","e1071","randomForest","foreach","import")
    hasParallelVersion <- TRUE

    getEBV <- function(gen){
        M = as.data.frame(pullSegSiteGeno(gen))
        colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
        EBV <- as.numeric(predict(rf_fit, M))
        as.matrix(EBV)
    }
}

if (args$model == "rf with stratified clustering training set"){
    fileTrain <- "RF_SC.R"
    modelLibs <- c("caret","ranger","tidyverse","e1071","randomForest","foreach","import")
    hasParallelVersion <- TRUE

    getEBV <- function(gen){
        M = as.data.frame(pullSegSiteGeno(gen))
        colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
        EBV <- as.numeric(predict(rf_fit, M))
        as.matrix(EBV)
    }
}

if (args$model == "rrblup with randomly sampled set") {
    fileTrain <- "RRBLUP_RD.R"
    modelLibs <- c("rrBLUP","devtools","dplyr","tidyverse","ggplot2","cluster","factoextra")
    hasParallelVersion <- FALSE

    getEBV <- function(gen){
        genMat <- pullSegSiteGeno(gen) 
        genMat <- genMat-1
        genMat %*% markerEffects
    }
}

if (args$model == "rrblup with stratified clustering training set") {
    fileTrain <- "RRBLUP_SC.R"
    modelLibs <- c("rrBLUP","devtools","dplyr","tidyverse","ggplot2","cluster","factoextra")
    hasParallelVersion <- FALSE

    getEBV <- function(gen){
        genMat <- pullSegSiteGeno(gen) 
        genMat <- genMat-1
        genMat %*% markerEffects
    }
}

if (args$model == "svm with randomly sampled training set"){
    fileTrain <- "SVM_RD.R"
    modelLibs <- c("e1071")
    hasParallelVersion <- FALSE

    getEBV <- function(gen){
        M = as.data.frame(pullSegSiteGeno(gen))
        colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
        EBV <- as.numeric(predict(SVMfit, M))
        as.matrix(EBV)
    }
}

if (args$model == "svm with stratified clustering training set"){
    fileTrain <- "SVM_SC.R"
    modelLibs <- c("e1071")
    hasParallelVersion <- FALSE

    getEBV <- function(gen){
        M = as.data.frame(pullSegSiteGeno(gen))
        colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
        EBV <- as.numeric(predict(SVMfit, M))
        as.matrix(EBV)
    }
}
