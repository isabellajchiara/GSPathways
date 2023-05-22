## define variables ##

library(purrr)

nModels = 7
nReps = 25
nGen = 10
nVar = 10

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

i = 1
repeat{
  source("1CycleOne_rrblup.R") ##Source the SCript for the SCenario you would like to run##
  
  
  geneticvaluesC1[,i] <- gvMatC1
  geneticvaluesC2[,i] <- gvMatC2
  geneticvaluesC3[,i] <- gvMatC3
  
  correlationsC1[,i] <- corMatC1
  correlationsC2[,i] <- corMatC2
  correlationsC3[,i] <- corMatC3
  
  variancesC1[,i] <- varMatC1
  variancesC2[,i] <- varMatC2
  variancesC3[,i] <- varMatC3
  
  allelesC1[[i]] <- allelesMatC1
  allelesC2[[i]] <- allelesMatC2
  allelesC3[[i]] <- allelesMatC3
  
  
  bv_ebvC1[[i]] <- bv_ebvC1
  bv_ebvC2[[i]] <- bv_ebvC2
  bv_ebvC3[[i]] <- bv_ebvC3
  
  i <- i + 1
  
  if (i > nReps){ ##break at number of desired reps##
    break
  }
  
  
  
  ##create data frames and label##
  geneticvaluesC1 <- as.data.frame(geneticvaluesC1)
  colnames(geneticvaluesC1) <- c(1:nReps)
  rownames(geneticvaluesC1) <- c("PrevCycPYT","NewParents","F1","F2","F3","F4","F5","PYT","AYT","Variety")
  
  geneticvaluesC2 <- as.data.frame(geneticvaluesC2)
  colnames(geneticvaluesC2) <- c(1:nReps)
  rownames(geneticvaluesC2) <- c("PrevCycPYT","NewParents","F1","F2","F3","F4","F5","PYT","AYT","Variety")
  
  geneticvaluesC3 <- as.data.frame(geneticvaluesC3)
  colnames(geneticvaluesC3) <- c(1:nReps)
  rownames(geneticvaluesC3) <- c("PrevCycPYT","NewParents","F1","F2","F3","F4","F5","PYT","AYT","Variety")
  
  correlationsC1 <- as.data.frame(correlationsC1)
  colnames(correlationsC1) <- c(1:nReps)
  rownames(correlationsC1) <- c("PrevCycPYT", "F2","F3","F4","F5","PYT","AYT")
  
  correlationsC2 <- as.data.frame(correlationsC2)
  colnames(correlationsC2) <- c(1:nReps)
  rownames(correlationsC2) <- c("PrevCycPYT", "F2","F3","F4","F5","PYT","AYT")
  
  correlationsC3 <- as.data.frame(correlationsC3)
  colnames(correlationsC3) <- c(1:nReps)
  rownames(correlationsC3) <- c("PrevCycPYT", "F2","F3","F4","F5","PYT","AYT")
  
  variancesC1 <- as.data.frame(variancesC1)
  colnames(variancesC1) <- c(1:nReps)
  rownames(variancesC1) <- c("PrevCycPYT", "newParents","F1","F2", "F3","F4", "F5", "PYT","AYT",
                           "Variety")
  
  variancesC2 <- as.data.frame(variancesC2)
  colnames(variancesC2) <- c(1:nReps)
  rownames(variancesC2) <- c("PrevCycPYT", "newParents","F1","F2", "F3","F4", "F5", "PYT","AYT",
                             "Variety")
  
  variancesC3 <- as.data.frame(variancesC3)
  colnames(variancesC3) <- c(1:nReps)
  rownames(variancesC3) <- c("PrevCycPYT", "newParents","F1","F2", "F3","F4", "F5", "PYT","AYT",
                             "Variety")
  
  
  
  
  ##write files
  write.csv(geneticvaluesC1, "1C1_rrblup_rd_gvs_snp_yield.csv")
  write.csv(geneticvaluesC2, "1C2_rrblup_rd_gvs_snp_yield.csv")
  write.csv(geneticvaluesC3, "1C3_rrblup_rd_gvs_snp_yield.csv")
  
  
  write.csv(correlationsC1, "1C1_rrblup_rd_cors_snp_yield.csv")
  write.csv(correlationsC2, "1C2_rrblup_rd_cors_snp_yield.csv")
  write.csv(correlationsC3, "1C3_rrblup_rd_cors_snp_yield.csv")
  
  write.csv(variancesC1, "1C1_rrblup_rd_vars_snp_yield.csv")
  write.csv(variancesC2, "1C2_rrblup_rd_vars_snp_yield.csv")
  write.csv(variancesC3, "1C3_rrblup_rd_vars_snp_yield.csv")
  
  saveRDS(allelesC1, file="1C1rrblup_rd_alleles_snp_yield.rds")
  saveRDS(allelesC2, file="1C2rrblup_rd_alleles_snp_yield.rds")
  saveRDS(allelesC3, file="1C3rrblup_rd_alleles_snp_yield.rds")
  
  saveRDS(bv_ebvC1, file="1C1rrblup_rd_bvebv_snp_yield.rds")
  saveRDS(bv_ebvC2, file="1C2rrblup_rd_bvebv_snp_yield.rds")
  saveRDS(bv_ebvC3, file="1C3rrblup_rd_bvebv_snp_yield.rds")
  
  
  
}
  
 
