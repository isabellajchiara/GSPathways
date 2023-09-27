# Defines some variables and functions based on what model was chosen

# ADD NEW MODELS HERE
if (args$model == "rf_random"){
    fileTrain <- "RF_RD.R"
    modelLibs <- c("caret","ranger","tidyverse","e1071","randomForest","foreach","import")
    modelParallelism <- TRUE

    getEBV <- function(gen){
        M = as.data.frame(pullSegSiteGeno(gen))
        colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
        EBV <- as.numeric(predict(rf_fit, M))
        as.matrix(EBV)
    }
}

if (args$model == "rf_stratifiedclusters"){
    fileTrain <- "RF_SC.R"
    modelLibs <- c("caret","ranger","tidyverse","e1071","randomForest","foreach","import","factoextra")
    hasParallelVersion <- TRUE

    getEBV <- function(gen){
        M = as.data.frame(pullSegSiteGeno(gen))
        colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
        EBV <- as.numeric(predict(rf_fit, M))
        as.matrix(EBV)
    }
}

if (args$model == "rrblup_random") {
    fileTrain <- "RRBLUP_RD.R"
    modelLibs <- c("rrBLUP","devtools","dplyr","tidyverse","ggplot2","cluster","factoextra")
    modelParallelism <- FALSE

    getEBV <- function(gen){
        genMat <- pullSegSiteGeno(gen) 
        genMat <- genMat-1
        genMat %*% markerEffects
    }
}

if (args$model == "rrblup_stratifiedclusters") {
    fileTrain <- "RRBLUP_SC.R"
    modelLibs <- c("rrBLUP","devtools","dplyr","tidyverse","ggplot2","cluster","factoextra")
    hasParallelVersion <- FALSE

    getEBV <- function(gen){
        genMat <- pullSegSiteGeno(gen) 
        genMat <- genMat-1
        genMat %*% markerEffects
    }
}

if (args$model == "svm_random"){
    fileTrain <- "SVM_RD.R"
    modelLibs <- c("e1071")
    modelParallelism <- FALSE

    getEBV <- function(gen){
        M = as.data.frame(pullSegSiteGeno(gen))
        colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
        EBV <- as.numeric(predict(SVMfit, M))
        as.matrix(EBV)
    }
}

if (args$model == "svm_stratifiedclusters"){
    fileTrain <- "SVM_SC.R"
    modelLibs <- c("e1071", "factoextra")
    hasParallelVersion <- FALSE

    getEBV <- function(gen){
        M = as.data.frame(pullSegSiteGeno(gen))
        colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
        EBV <- as.numeric(predict(SVMfit, M))
        as.matrix(EBV)
    }
}

if (args$model == "ann_random"){
    fileTrain <- "ANN_RD.R"
    modelLibs <- c("tidyverse","keras","tensorflow","readr","devtools")
    modelParallelism <- FALSE

    getEBV <- function(gen){
        geno <- as.matrix(pullSegSiteGeno(gen))
        predict(model,geno)
    }
}

if (args$model == "ann_stratifiedclusters"){
    fileTrain <- "ANN_SC.R"
    modelLibs <- c("tidyverse","keras","tensorflow","readr","devtools","writexl","factoextra")
    modelParallelism <- FALSE

    getEBV <- function(gen){
        geno <- as.matrix(pullSegSiteGeno(gen))
        predict(model,geno)
    }
}
