suppressMessages(library(argparse))
suppressMessages(library(tictoc))
suppressMessages(library(doParallel))
source("ParameterSettings.R")
source("InterfaceLibrary.R")
source("FunctionsLibrary.R")

DATA_DIR <- "data"
MODEL_DIR <- "models"

args <- parseArgs()
if (args$noInteraction == FALSE)
  interactive_menu()

# validateArgs(args)

## define variables ##

nModels = 7
# nReps = args$nReps
nGen = 10
nVar = 9
# nCores = args$nCores

# model = args$model
# nCycles = args$nCycles
# trainGen = args$trainGen

## Create model definitions
source("DefineModelVariables.R")

## establish empty matrices to hold outputs for Selfing and Recombination Population ##

## Run repeat loop to run reps ##
if (args$nCores == 1 || hasParallelVersion) { # Run reps serially
  cli_alert_info("Importing simulation libraries...")

  res <- lapply(1:args$nReps, function(rep){
    source("SimplifiedBreedingCyclePipeline.R") ##Source the SCript for the SCenario you would like to run##
    ret
  })
} else { # Run reps in parallel
  cli_alert_info("Running {args$nReps} reps in {args$nCores} cores...")

  # Create parallel cluster and export variables
  cl <- makeCluster(args$nCores)
  clusterExport(cl, c("args", "loadModelLibs", "DATA_DIR", "MODEL_DIR"))

  res <- parLapply(cl, 1:args$nReps, function(rep){
    source("SimplifiedBreedingCyclePipeline.R") ##Source the SCript for the SCenario you would like to run##
    ret
  })

  # Stop parallel cluster
  stopCluster(cl)
}
res <- bindSimResults(res)

cli_alert_success("Simulation finished!")
cli_text()


##create results directory and enter it##
dirName <- args$outputDir
dir.create(file.path(dirName))

workingDir <- getwd()
setwd(file.path(dirName))

##create all output files##
Allgeneticvalues <- list()
for (cycle in 1:args$nCycles){
  cli_alert_info("Writing output files for cycle {cycle}...")
  Allgeneticvalues[[cycle]] <- getAllGeneticValues(res$geneticvalues[[cycle]], 10, 2)
  res$correlations[[cycle]] <- getCorrelations(res$correlations[[cycle]])
  res$variances[[cycle]] <- getVariances(res$variances[[cycle]])

  write.csv(Allgeneticvalues[[cycle]], paste("1C", cycle, "_", args$model, "_rd_gvs_snp_yield.csv", sep=""))
  write.csv(res$correlations[[cycle]], paste("1C", cycle, "_", args$model,"_rd_cors_snp_yield.csv", sep=""))
  write.csv(res$variances[[cycle]], paste("1C", cycle, "_", args$model,"_rd_vars_snp_yield.csv", sep=""))
  saveRDS(res$alleles[[cycle]], file=paste("1C", cycle, "_", args$model,"_rd_alleles_snp_yield.rds", sep=""))
  saveRDS(res$bv_ebv[[cycle]], file=paste("1C", cycle, "_", args$model,"_rd_bvebv_snp_yield.rds", sep=""))
}

cli_text()
cli_alert_success("Results saved!")
setwd(workingDir) # Go back to previous directory
