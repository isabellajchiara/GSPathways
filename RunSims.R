suppressMessages(library(argparse))
suppressMessages(library(tictoc))
suppressMessages(library(reticulate))
suppressMessages(library(doParallel))
suppressMessages(library(AlphaSimR))
source("ParameterSettings.R")
source("InterfaceLibrary.R")
source("FunctionsLibrary.R")

# For the ANN model
args <- parseArgs()

activeLog <- args$nCores == 1 || hasParallelVersion

Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1) 

DATA_DIR <- "data"
MODEL_DIR <- "models"


if (args$noInteraction == FALSE)
  interactive_menu()


## define variables ##
genZero <- list()


allTrainingDataGeno = list()
allTrainingDataPheno = list()
nModels = 7
nGen = 10
nVar = 9

#load data and establish founder pop at outset so we do not reload data with every rep

genMap <- readRDS(file.path(DATA_DIR, "genMap4616SNPs.rds")) # can load other genZeroMaps 
haplotypes <- readRDS(file.path(DATA_DIR, "doubleHaplo4616.rds")) # can load other genZerootype data, must match genZeroMap


founderPop = newMapPop(genMap, 
                       haplotypes, 
                       inbred = FALSE, 
                       ploidy = 2L)
defineTraitAEG(4,3.6,0.14) # nQtl per chr, mean,heritability

## CYCLE 0 TO BUILD INITIAL TRAINING POP = Phenotypic Selection

# create base pop randomly cross 200 parents 
Base = newPop(founderPop)
Base = setPheno(Base)

newParents <- selectNewParents(Base,5,"pheno")
genZero$F1 = randCross(newParents, 200, nProgeny=3)

if (activeLog)
  cli_text("generated F1s...")

# self and bulk genZero$F1 to form genZero$F2 
genZero$F2 = self(genZero$F1, nProgeny = 20)
genZero$F2 = setPheno(genZero$F2)

if (activeLog)
  cli_text("generated F2s...")

# select top individuals from each family to form genZero$F2. Bulk and self to form genZero$F3
genZero$F3 = TopWithinFam(genZero$F2,5,200,"pheno") # top 10 F2 fam, 100 ind per fam using pheno
genZero$F3 = setPheno(genZero$F3)

if (activeLog)
  cli_text("generated F3s...")

# select top individuals within genZero$F3 families to form genZero$F4 
genZero$F4 = TopWithinFam(genZero$F3,5,30,"pheno") # top 5 F3 fam, 50 ind per fam using pheno
genZero$F4 = setPheno(genZero$F4)

if (activeLog)
  cli_text("generated F4s...")

# select top families from genZero$F4 to form genZero$F5 
genZero$F5 = TopFamily(genZero$F4,4,"pheno") #select top 4 F4 families 
genZero$F5 = setPheno(genZero$F5)


if (activeLog)
  cli_text("generated F5s...")


# collect data for next cycle's model training
trainingGenotypes = list()
trainingPhenotypes = list()

x = 1
for (generation in 1:length(genZero)){
  set = genZero[[generation]]
  M = pullSegSiteGeno(set)
  y = pheno(set)
  trainingGenotypes[[x]] = M
  trainingPhenotypes[[x]] = y
  x = x+1
}
allTrainingDataGeno[[1]] = trainingGenotypes
allTrainingDataPheno[[1]] = trainingPhenotypes


## Run repeat loop to run reps ##
if (args$nCores == 1 || modelParallelism) { # Run reps serially
  cli_alert_info("Importing simulation libraries...")

  res <- lapply(1:args$nReps, function(rep){
    cli_alert_info("Simulating rep {rep}/{args$nReps}...")
    source("SimplifiedBreedingCyclePipeline.R") ##Source the SCript for the SCenario you would like to run##
    cli_text("Rep {rep} finished.")

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

tic()
##create all output files##
Allgeneticvalues <- list()
for (cycle in 1:args$nCycles){
  cli_alert_info("Writing output files for cycle {cycle}...")
  Allgeneticvalues[[cycle]] <- getAllGeneticValues(res$geneticvalues[[cycle]])
  res$correlations[[cycle]] <- getCorrelations(res$correlations[[cycle]])
  res$variances[[cycle]] <- getVariances(res$variances[[cycle]])

  write.csv(Allgeneticvalues[[cycle]], paste("C", cycle, "_", args$model,"_trainAt",args$trainGen,"_trainWith",args$trainingData,"_",args$parentSelections, "Parents_gvs_snp_yield.csv", sep=""))
  write.csv(res$correlations[[cycle]], paste("C", cycle, "_",args$model,"_trainAt",args$trainGen,"_trainWith",args$trainingData,"_",args$parentSelections,"Parents_cors_snp_yield.csv", sep=""))
  write.csv(res$variances[[cycle]], paste("C", cycle, "_",args$model,"_trainAt",args$trainGen,"_trainWith",args$trainingData,"_",args$parentSelections,"Parents_vars_snp_yield.csv", sep=""))
  saveRDS(res$pheno[[cycle]], paste("C", cycle, "_",args$model,"_trainAt",args$trainGen,"_trainWith",args$trainingData,"_",args$parentSelections,"Parents_pheno_snp_yield.rds", sep=""))
  saveRDS(res$alleles[[cycle]], file=paste("C", cycle, "_",args$model,"_trainAt",args$trainGen,"_trainWith",args$trainingData,"_",args$parentSelections,"Parents_alleles_snp_yield.rds", sep=""))
}


cli_text("Time taken to write results:")
toc()
cli_text()
cli_alert_success("Results saved!")
setwd(workingDir) # Go back to previous directory
