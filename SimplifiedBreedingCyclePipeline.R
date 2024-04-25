## PEDIGREE BREEDING METHOD USING GEBVs TO SELECT

#load requirements 
suppressMessages(library(AlphaSimR))
suppressMessages(library(argparse))
suppressMessages(library(tictoc))
suppressMessages(library(reticulate))
suppressMessages(library(doParallel))
suppressMessages(library(AlphaSimR))
source("ParameterSettings.R")
source("InterfaceLibrary.R")
source("FunctionsLibrary.R")
source("ModelVariables.R")
loadModelLibs()


activeLog <- args$nCores == 1 || hasParallelVersion

if (activeLog)
  cli_text("Generating parent population...")
## define variables ##

ret <- list( 
  geneticvalues = list(),
  correlations = list(),
  variances = list(),
  alleles = list(),
  pheno = list()
)

# create results matrices
gvMat <- matrix(nrow=nGen, ncol=1)
corMat <- matrix(nrow=nModels, ncol=1)
varMat <- matrix(nrow=nVar, ncol=1)
allelesMat <- NULL
valuesMat <- NULL


gen <- list()

if (activeLog)
  cli_text("Created outputs")

# establish simulation parameters

if (activeLog)
  cli_text("Generated F5s...")


if (activeLog)
  cli_text("Initial cycle complete")

if (activeLog)
  cli_text("Generating parent population...")



# INITAL TRAINING POP IS BUILT, START GS CYCLES

for (cycle in 1:args$nCycles){
  if (activeLog)
    cli_text("Running cycle {cycle}/{args$nCycles}...")
  
      if (cycle == 1) {
        
          trainModel()
          genZero[[args$parentSelections]]@ebv = getEBV(genZero[[args$parentSelections]]) # set EBVs for parent pool
          
          updateResults(1, genZero[[args$parentSelections]], "ParentPool")
        
          #select parents based on EBV
          newParents <- selectNewParents(genZero[[args$parentSelections]], 5, "ebv") 

          if (grepl("rrblup",args$model)==TRUE){
            corMat[1,] = cor(bv(genZero[[args$parentSelections]]), ebv(genZero[[args$parentSelections]])) #determine model performance
          }
    
          if (grepl("ann",args$model)==TRUE){
            corMat[1,] = cor(pheno(genZero[[args$parentSelections]]), ebv(genZero[[args$parentSelections]])) 
          } #determine model performance

          #collect data on new parents
          updateResults(2, newParents, "NP")
  
      } else { #for cycles 2+ we don't need to train the model, we already have it from last cycle

            #select new parents
            newParents <- selectNewParents(gen[[args$parentSelections]], 5, "ebv") 
        
            if (grepl("rrblup",args$model)==TRUE){
            corMat[1,] = cor(bv(gen$F2), ebv(gen$F2)) #determine model performance
            }
  
            if (grepl("ann",args$model)==TRUE){
            corMat[1,] = cor(pheno(gen$F2), ebv(gen$F2)) 
            } #determine model performance

            #collect data on new parents
            updateResults(2, newParents, "NP")
            }

  ## 200 random crosses of new parents
  gen$F1 = randCross(newParents, 200,nProgeny=3)

  updateResults(3, gen$F1, "F1")
  
  
  ## self and bulk gen$F1 to form gen$F2 ##
  if (args$trainGen == "F2")
    trainModel()
  
  gen$F2 = self(gen$F1, nProgeny = 20) 
  gen$F2@ebv = getEBV(gen$F2)   ## set EBV 

  updateResults(4, gen$F2, "F2")
  if (grepl("rrblup",args$model)==TRUE){
    corMat[2,] = cor(bv(gen$F2), ebv(gen$F2)) #determine model performance
  }
  
  if (grepl("ann",args$model)==TRUE){
    corMat[2,] = cor(pheno(gen$F2), ebv(gen$F2)) 
  } #determine model performance
  
  ## select top individuals from gen$F2 bulk to form gen$F3 
  if (args$trainGen == "F3")
    trainModel()
  
  gen$F3 = TopWithinFam(gen$F2, 5,200 , "ebv")
  gen$F3 = setPheno(gen$F3)
  gen$F3@ebv = getEBV(gen$F3) #set EBVS

  updateResults(5, gen$F3, "F3")
    

    if (grepl("rrblup",args$model)==TRUE){
    corMat[3,] = cor(bv(gen$F3), ebv(gen$F3)) #determine model performance
  }
  
  if (grepl("ann",args$model)==TRUE){
    corMat[2,] = cor(pheno(gen$F3), ebv(gen$F3)) 
  } #determine model performance
  
  ## select top within familiy from gen$F3 to form gen$F4 
  if (args$trainGen == "F4")
    trainModel("F4", gen$F4)
  
  gen$F4 = TopWithinFam(gen$F3, 5, 30, "ebv")
  gen$F4 = setPheno(gen$F4)
  gen$F4@ebv = getEBV(gen$F4) #set EBV

  updateResults(6, gen$F4, "F4")
  
   if (grepl("rrblup",args$model)==TRUE){
    corMat[4,] = cor(bv(gen$F4), ebv(gen$F4)) #determine model performance
  }
  
  if (grepl("ann",args$model)==TRUE){
    corMat[4,] = cor(pheno(gen$F4), ebv(gen$F4)) 
  } #determine model performance

  ## select top families from gen$F4 to form gen$F5 ##
  if (args$trainGen == "F5")
    trainModel()
  
  gen$F5 = TopFamily(gen$F4,4,"ebv")
  gen$F5 = setPheno(gen$F5)
  gen$F5@ebv = getEBV(gen$F5) #set EBV

  updateResults(7, gen$F5, "F5")

    if (grepl("rrblup",args$model)==TRUE){
    corMat[5,] = cor(bv(gen$F5), ebv(gen$F5)) #determine model performance
  }
  
  if (grepl("ann",args$model)==TRUE){
    corMat[5,] = cor(pheno(gen$F5), ebv(gen$F5)) 
  } #determine model performance

  ## select top gen$F5 families for preliminary yield trial ##
  gen$PYT = TopFamily(gen$F5,2,"ebv")
  gen$PYT = setPheno(gen$PYT, reps=2)
  gen$PYT@ebv = getEBV(gen$PYT) #set EBV

  updateResults(8, gen$PYT, "PYT")
  
    if (grepl("rrblup",args$model)==TRUE){
    corMat[6,] = cor(bv(gen$PYT), ebv(gen$PYT)) #determine model performance
  }
  
  if (grepl("ann",args$model)==TRUE){
    corMat[6,] = cor(pheno(gen$PYT), ebv(gen$PYT)) 
  } #determine model performance

  ## select top families from gen$PYT for gen$AYT ##
  gen$AYT = TopFamily(gen$PYT, 1, "ebv")
  gen$AYT = setPheno(gen$AYT, reps=5)
  gen$AYT@ebv = getEBV(gen$AYT)

  updateResults(9, gen$AYT, "AYT")
    
    if (grepl("rrblup",args$model)==TRUE){
    corMat[7,] = cor(bv(gen$AYT), ebv(gen$AYT)) #determine model performance
  }
  
  if (grepl("ann",args$model)==TRUE){
    corMat[7,] = cor(pheno(gen$AYT), ebv(gen$AYT)) 
  } #determine model performance
  
  ## select top plants to form variety ##
  VarietySel = selectInd(gen$AYT, 1, use="ebv")
  Variety = self(VarietySel)
  
  gvMat[10,] <- mean(gv(Variety))
  allelesMatVar <- getAllelesMat(Variety, "Variety")
  allelesMat <- rbind(allelesMat, allelesMatVar)
  
  curValuesMat = matrix(nrow=nInd(Variety),ncol=5)
  curValuesMat[,1] <- as.matrix(rep(paste0("Variety","C",cycle,sep=""), times=nInd(Variety)))
  curValuesMat[,2] <- pheno(Variety)
  curValuesMat[,3] <- gv(Variety)
  curValuesMat[,4] <- bv(Variety)
  curValuesMat[,5] <- NA
  colnames(curValuesMat) = c("gen","pheno","gv","tbv","ebv")
  valuesMat <<- rbind(valuesMat, curValuesMat)
  valuesMat = as.data.frame(valuesMat)

  if (activeLog)
    cli_alert_success("Finished cycle {cycle}/{args$nCycles}...")
  
  # collect data to be used for next cycle's training
  
  trainingGenotypes = list()
  trainingPhenotypes = list()
  x = 1
  for (generation in 1:length(gen)){
    set = gen[[generation]]
    M = pullSegSiteGeno(set)
    y = pheno(set)
    trainingGenotypes[[x]] = M
    trainingPhenotypes[[x]] = y
    x = x+1
  }
  
  allTrainingDataGeno[[cycle+1]] <- trainingGenotypes
  allTrainingDataPheno[[cycle+1]] <- trainingPhenotypes
  
  
  if (activeLog)
    cli_alert_success("Collected training data for next cycle")
    
  
  ret$geneticvalues[[cycle]] <- gvMat
  ret$correlations[[cycle]] <- corMat
  ret$variances[[cycle]] <- varMat
  ret$alleles[[cycle]] <- allelesMat
  ret$pheno[[cycle]] <- valuesMat
  
  
}
