
## define variables ##

nModels = 7
nReps = 1
nGen = 10
nVar = 9

nCycles = 1

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

##create data frames and label##
Allgeneticvalues <- list()

for (cycle in paste("C", 1:nCycles, sep="")){
  Allgeneticvalues[[cycle]] <- getAllGeneticValues(geneticvalues[[cycle]], 10, 2)
  correlations[[cycle]] <- getCorrelations(correlations[[cycle]])
  variances[[cycle]] <- getVariances(variances[[cycle]])

  write.csv(Allgeneticvalues[[cycle]], paste("1", cycle, "_rrblup_rd_gvs_snp_yield.csv", sep=""))
  write.csv(correlations[[cycle]], paste("1", cycle, "_rrblup_rd_cors_snp_yield.csv", sep=""))
  write.csv(variances[[cycle]], paste("1", cycle, "_rrblup_rd_vars_snp_yield.csv", sep=""))
  saveRDS(alleles[[cycle]], file=paste("1", cycle, "rrblup_rd_alleles_snp_yield.rds", sep=""))
  saveRDS(bv_ebv[[cycle]], file=paste("1", cycle, "rrblup_rd_bvebv_snp_yield.rds", sep=""))
}

