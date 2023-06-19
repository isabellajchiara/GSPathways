methods <- list()

methods$rf <- list()
methods$rf$fileTrain <- "RF_RD.R"
methods$rf$fileRetrain <- "RF_RD_retrain.R"
methods$rf$getEBV <- GetEBVrf

methods$rrblup <- list()
methods$rrblup$fileTrain <- "rrblup_sc.R"
methods$rrblup$fileRetrain <- "rrblup_sc_retrain.R"
methods$rrblup$getEBV <- GetEBVrrblup
