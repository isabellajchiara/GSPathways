nGVs = 8
nCors = 6
nVars= 7
nReps = 2


## establish empty matrices to hold outputs 

C1GV = matrix(nrow=nGVs, ncol=nReps)
C2GV = matrix(nrow=nGVs, ncol=nReps)
C3GV = matrix(nrow=nGVs, ncol=nReps)

C1COR = matrix(nrow=nCors, ncol=nReps)
C2COR = matrix(nrow=nCors, ncol=nReps)
C3COR = matrix(nrow=nCors, ncol=nReps)

C1VAR = matrix(nrow=nVars, ncol=nReps)
C2VAR = matrix(nrow=nVars, ncol=nReps)
C3VAR = matrix(nrow=nVars, ncol=nReps)


## Run repeat loop to run reps ##

i = 1
repeat{
  source("1MultipleCyclesF2.R") ##Source the SCript for the SCenario you would like to run##
  
  C1GV[,i] <- cycle1gvs
  C2GV[,i] <- cycle2gvs
  C3GV[,i] <- cycle3gvs
  
  C1COR[,i] <- cycle1cors
  C2COR[,i] <- cycle2cors
  C3COR[,i] <- cycle3cors
  
  C1VAR[,i] <- cycle1vars
  C2VAR[,i] <- cycle2vars
  C3VAR[,i] <- cycle3vars
  
  i <- i + 1
  
  if (i > nReps){ ##break at number of desired reps##
    break
  }
  
  
  
  ##create data frames and label##
  cycle1gvs <- as.data.frame(cycle1gvs)
  C1GVmean <- rowMeans(cycle1gvs)
  
  cycle2gvs <- as.data.frame(cycle2gvs)
  C2GVmean <- rowMeans(cycle2gvs)
  
  cycle3gvs <- as.data.frame(cycle3gvs)
  C3GVmean <- rowMeans(cycle3gvs)
  
  AllCyclesGV <- cbind(C1GVmean, C2GVmean, C3GVmean )
  write.csv(AllCyclesGV,"MCF2_GVS_rrblup_RD_SNP_yield.csv")
  
  cycle1cors <- as.data.frame(cycle1cors)
  C1CORmean <- rowMeans(cycle1cors)
  cycle2cors <- as.data.frame(cycle2cors)
  C2CORmean <- rowMeans(cycle2cors)
  cycle3cors <- as.data.frame(cycle3cors)
  C3CORmean <- rowMeans(cycle3cors)
  
  AllCyclesCOR <- cbind(C1CORmean, C2CORmean, C3CORmean )
  write.csv(AllCyclesCOR,"MCF2_CORS_rrblup_RD_SNP_yield.csv")
  
  
  cycle1vars <- as.data.frame(cycle1vars)
  C1VARmean <- rowMeans(cycle1vars)
  cycle2vars <- as.data.frame(cycle2vars)
  C2VARmean <- rowMeans(cycle2vars)
  cycle3vars <- as.data.frame(cycle3vars)
  C3VARmean <- rowMeans(cycle3vars)
  
  AllCyclesVAR <- cbind(C1VARmean, C2VARmean, C3VARmean )
  write.csv(AllCyclesVAR,"MCF2_VARS_rrblup_RD_SNP_yield.csv")
  
  
  
  
  
}

################################################################################
