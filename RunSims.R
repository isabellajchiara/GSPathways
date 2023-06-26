library(argparse)
library(tictoc)
source("FunctionsLibrary.R")

DATA_DIR <- "data"
MODEL_DIR <- "models"

args <- parseArgs()
validateArgs(args)

## define variables ##

nModels = 7
nReps = args$nReps
nGen = 10
nVar = 9
nCores = args$nCores

model = args$model
nCycles = args$nCycles
trainGen = args$trainGen

## Create model definitions
source("DefineModelVariables.R")

## establish empty matrices to hold outputs for Selfing and Recombination Population ##

geneticvalues <- list()
correlations <- list()
variances <- list()
alleles <- list()
bv_ebv <- list()

for (cycle in paste("C", 1:nCycles, sep="")){
  geneticvalues[[cycle]] <- matrix(nrow=nGen, ncol=nReps)
  correlations[[cycle]] <- matrix(nrow=nModels, ncol=nReps)
  variances[[cycle]] <- matrix(nrow=nVar,ncol=nReps)
  alleles[[cycle]] <- vector("list", length = nReps)
  bv_ebv[[cycle]] <- vector("list", length = nReps)
}

## Run repeat loop to run reps ##

for (rep in 1:nReps){
  source("SimplifiedBreedingCyclePipeline.R") ##Source the SCript for the SCenario you would like to run##
}

##create results directory and enter it##
dirName <- args$outputDir
if (is.null(args$outputDir))
  dirName <- getDirName(model)
dir.create(file.path(dirName))

workingDir <- getwd()
setwd(file.path(dirName))

##create all output files##
Allgeneticvalues <- list()
for (cycle in paste("C", 1:nCycles, sep="")){
  Allgeneticvalues[[cycle]] <- getAllGeneticValues(geneticvalues[[cycle]], 10, 2)
  correlations[[cycle]] <- getCorrelations(correlations[[cycle]])
  variances[[cycle]] <- getVariances(variances[[cycle]])

  write.csv(Allgeneticvalues[[cycle]], paste("1", cycle, "_", model, "_rd_gvs_snp_yield.csv", sep=""))
  write.csv(correlations[[cycle]], paste("1", cycle, "_", model,"_rd_cors_snp_yield.csv", sep=""))
  write.csv(variances[[cycle]], paste("1", cycle, "_", model,"_rd_vars_snp_yield.csv", sep=""))
  saveRDS(alleles[[cycle]], file=paste("1", cycle, "_", model,"_rd_alleles_snp_yield.rds", sep=""))
  saveRDS(bv_ebv[[cycle]], file=paste("1", cycle, "_", model,"_rd_bvebv_snp_yield.rds", sep=""))
}

setwd(workingDir) # Go back to previous directory
