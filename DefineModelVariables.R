# Defines some variables and functions based on what is set 
# in the model chosen

# - rf: Random Forest
# - rrblup: RRBLUP

# ADD NEW MODELS HERE
if (model == "rf"){
    print("Model chosen: Random Forest")
    fileTrain <- "RF_RD.R"
    fileRetrain <- "RF_RD_retrain.R"

    getEBV <- function(gen){
        M = as.data.frame(pullSegSiteGeno(gen))
        colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
        EBV <- as.numeric(predict(rf_fit, M))
        as.matrix(EBV)
    }
}

if (model == "rrblup") {
    print("Model chosen: RRBLUP")
    fileTrain <- "rrblup_sc.R"
    fileRetrain <- "rrblup_sc_retrain.R"

    getEBV <- function(gen){
        genMat <- pullSegSiteGeno(gen) 
        genMat <- genMat-1
        genMat %*% markerEffects
    }
}
