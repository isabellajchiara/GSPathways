# Defines some variables and functions based on what model was chosen

# ADD NEW MODELS HERE
if (args$model == "rf"){
    cli_alert_info("args$model chosen: Random Forest")
    fileTrain <- "RF_RD.R"
    fileRetrain <- "RF_RD_retrain.R"
    modelLibs <- c("caret","ranger","tidyverse","e1071","randomForest","foreach","import")
    isParallel <- TRUE

    getEBV <- function(gen){
        M = as.data.frame(pullSegSiteGeno(gen))
        colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
        EBV <- as.numeric(predict(rf_fit, M))
        as.matrix(EBV)
    }
}

if (args$model == "rrblup") {
    cli_alert_info("Model chosen: Ridge Regression")
    fileTrain <- "rrblup_sc.R"
    fileRetrain <- "rrblup_sc_retrain.R"
    modelLibs <- c("rrBLUP","devtools","dplyr","tidyverse","ggplot2","cluster","factoextra")
    isParallel <- FALSE

    getEBV <- function(gen){
        genMat <- pullSegSiteGeno(gen) 
        genMat <- genMat-1
        genMat %*% markerEffects
    }
}
