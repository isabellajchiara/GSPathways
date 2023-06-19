if (model == "rf"){
    print("Model chosen: Random Forest")
    fileTrain <- "RF_RD.R"
    fileRetrain <- "RF_RD_retrain.R"
    getEBV <- GetEBVrf
}
else if (model == "rrblup") {
    print("Model chosen: RRBLUP")
    fileTrain <- "rrblup_sc.R"
    fileRetrain <- "rrblup_sc_retrain.R"
    getEBV <- GetEBVrrblup
} else {
    print("Model not recognized. Defaulting to RRBLUP")
    fileTrain <- "rrblup_sc.R"
    fileRetrain <- "rrblup_sc_retrain.R"
    getEBV <- GetEBVrrblup
}


