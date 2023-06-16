
## define variables ##

nModels = 7
nReps = 1
nGen = 10
nVar = 9

## establish empty matrices to hold outputs for Selfing and Recombination Population ##

geneticvaluesC1 <- matrix(nrow=nGen, ncol=nReps)
geneticvaluesC2 <- matrix(nrow=nGen, ncol=nReps)
geneticvaluesC3 <- matrix(nrow=nGen, ncol=nReps)

correlationsC1 <- matrix(nrow=nModels, ncol=nReps)
correlationsC2 <- matrix(nrow=nModels, ncol=nReps)
correlationsC3 <- matrix(nrow=nModels, ncol=nReps)

variancesC1 <- matrix(nrow=nVar,ncol=nReps)
variancesC2 <- matrix(nrow=nVar,ncol=nReps)
variancesC3 <- matrix(nrow=nVar,ncol=nReps)

allelesC1 <- vector("list", length = nReps)
allelesC2 <- vector("list", length = nReps)
allelesC3 <- vector("list", length = nReps)

bv_ebvC1 <- vector("list", length = nReps)
bv_ebvC2 <- vector("list", length = nReps)
bv_ebvC3 <- vector("list", length = nReps)



## Run repeat loop to run reps ##

for (rep in 1:nReps){
  source("SimplifiedBreedingCyclePipeline.R") ##Source the SCript for the SCenario you would like to run##
  
  geneticvaluesC1[,rep] <- gvMat
  # geneticvaluesC2[,rep] <- gvMatC2
  # geneticvaluesC3[,rep] <- gvMatC3
  
  correlationsC1[,rep] <- corMat
  # correlationsC2[,rep] <- corMatC2
  # correlationsC3[,rep] <- corMatC3
  
  variancesC1[,rep] <- varMat
  # variancesC2[,rep] <- varMatC2
  # variancesC3[,rep] <- varMatC3
  
  allelesC1[[rep]] <- allelesMat
  # allelesC2[[rep]] <- allelesMatC2
  # allelesC3[[rep]] <- allelesMatC3
  
  bv_ebvC1[[rep]] <- bv_ebv
  # bv_ebvC2[[rep]] <- bv_ebvC2
  # bv_ebvC3[[rep]] <- bv_ebvC3
  
}

##create data frames and label##
AllgeneticvaluesC1 <- getAllGeneticValues(geneticvaluesC1, 10, 2)
# AllGeneticValuesC2 <- getAllGeneticValues(geneticvaluesC2, 10, 3)
# AllGeneticValuesC3 <- getAllGeneticValues(geneticvaluesC3, 10, 3)

correlationsC1 <- getCorrelations(correlationsC1)
# correlationsC2 <- getCorrelations(correlationsC2)
# correlationsC3 <- getCorrelations(correlationsC3)

variancesC1 <- getVariances(variancesC1)
# variancesC2 <- getVariances(variancesC2)
# variancesC3 <- getVariances(variancesC3)


##write files
write.csv(AllgeneticvaluesC1, "1C1_rrblup_rd_gvs_snp_yield.csv")
# write.csv(AllgeneticvaluesC2, "1C2_rrblup_rd_gvs_snp_yield.csv")
# write.csv(AllgeneticvaluesC3, "1C3_rrblup_rd_gvs_snp_yield.csv")

write.csv(correlationsC1, "1C1_rrblup_rd_cors_snp_yield.csv")
# write.csv(correlationsC2, "1C2_rrblup_rd_cors_snp_yield.csv")
# write.csv(correlationsC3, "1C3_rrblup_rd_cors_snp_yield.csv")

write.csv(variancesC1, "1C1_rrblup_rd_vars_snp_yield.csv")
# write.csv(variancesC2, "1C2_rrblup_rd_vars_snp_yield.csv")
# write.csv(variancesC3, "1C3_rrblup_rd_vars_snp_yield.csv")

saveRDS(allelesC1, file="1C1rrblup_rd_alleles_snp_yield.rds")
# saveRDS(allelesC2, file="1C2rrblup_rd_alleles_snp_yield.rds")
# saveRDS(allelesC3, file="1C3rrblup_rd_alleles_snp_yield.rds")

saveRDS(bv_ebvC1, file="1C1rrblup_rd_bvebv_snp_yield.rds")
# saveRDS(bv_ebvC2, file="1C2rrblup_rd_bvebv_snp_yield.rds")
# saveRDS(bv_ebvC3, file="1C3rrblup_rd_bvebv_snp_yield.rds")
  
