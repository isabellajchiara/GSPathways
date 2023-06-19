# Defines some variables and functions based on what is set 
# in the model chosen

# - rf: Random Forest
# - rrblup: RRBLUP

# ADD NEW MODELS HERE
if (model == "rf"){
    print("Model chosen: Random Forest")
    fileTrain <- "RF_RD.R"
    fileRetrain <- "RF_RD_retrain.R"
    getEBV <- GetEBVrf
}
if (model == "rrblup") {
    print("Model chosen: RRBLUP")
    fileTrain <- "rrblup_sc.R"
    fileRetrain <- "rrblup_sc_retrain.R"
    getEBV <- GetEBVrrblup
}


