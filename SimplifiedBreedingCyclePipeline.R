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

# data to be returned
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

allTrainingDataGeno = list()
allTrainingDataPheno = list()

      if (args$parentSelections == "F2"){# depending on parent selections, we will collect varying numbers of phenotypes
      nSim = 13366
      }

      if (args$parentSelections == "F5"){
      nSim = 1486
      }

valuesMat <- matrix(nrow=(nSim*args$nCycles),ncol=5) #gen, pheno, gv,tbv,ebv

gen <- list()

if (activeLog)
  cli_text("Created outputs")

# establish simulation parameters

defineTraitAEG(4,3.6,0.14) # nQtl per chr, mean,heritability

## CYCLE 0 TO BUILD INITIAL TRAINING POP = Phenotypic Selection

# create base pop randomly cross 200 parents 

Base = newPop(founderPop)
Base = setPheno(Base)

newParents <- selectNewParents(Base,5,"pheno")
gen$F1 = randCross(newParents, 200, nProgeny=3)

if (activeLog)
  cli_text("Generated F1s...")

# self and bulk gen$F1 to form gen$F2 

gen$F2 = self(gen$F1, nProgeny = 20)
gen$F2 = setPheno(gen$F2)

if (activeLog)
  cli_text("Generated F2s...")

# select top individuals from each family to form gen$F2. Bulk and self to form gen$F3

gen$F3 = TopWithinFam(gen$F2,5,200,"pheno") # top 10 F2 fam, 100 ind per fam using pheno
gen$F3 = setPheno(gen$F3)

if (activeLog)
  cli_text("Generated F3s...")

# select top individuals within gen$F3 families to form gen$F4 

gen$F4 = TopWithinFam(gen$F3,5,30,"pheno") # top 5 F3 fam, 50 ind per fam using pheno
gen$F4 = setPheno(gen$F4)

if (activeLog)
  cli_text("Generated F4s...")

# select top families from gen$F4 to form gen$F5 

gen$F5 = TopFamily(gen$F4,4,"pheno") #select top 4 F4 families 
gen$F5 = setPheno(gen$F5)

if (activeLog)
  cli_text("Generated F5s...")


if (args$parentSelections == "F2"){

gvMat[1,] <- mean(gv(gen$F2))
varMat[1,] <- varG(gen$F2)
phenos = pheno(gen$F2)
tbvs = bv(gen$F2)
gvs = gv(gen$F2)
valuesMat[1:nInd(gen$F2),1] = rep("ParentPool", times=nInd(gen$F2))
valuesMat[1:nInd(gen$F2),2] = phenos
valuesMat[1:nInd(gen$F2),2] = gvs
valuesMat[1:nInd(gen$F2),4] = tbvs
updateResults(1, gen$F2, "ParentPool")
}

if (args$parentSelections == "F5"){

gvMat[1,] <- mean(gv(gen$F5))
varMat[1,] <- varG(gen$F5)
phenos = pheno(gen$F5)
tbvs = bv(gen$F5)
gvs = gv(gen$F5)
valuesMat[1:nInd(gen$F5),1] = rep("ParentPool", times=nInd(gen$F2))
valuesMat[1:nInd(gen$F5),2] = phenos
valuesMat[1:nInd(gen$F5),2] = gvs
valuesMat[1:nInd(gen$F5),4] = tbvs
updateResults(1, gen$F5, "ParentPool")

}


# collect data for next cycle's model training

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
allTrainingDataGeno[[1]] = trainingGenotypes
allTrainingDataPheno[[1]] = trainingPhenotypes


if (activeLog)
  cli_text("Initial cycle complete")

# INITAL TRAINING POP IS BUILT, START GEBV CYCLES

for (cycle in 1:args$nCycles){
  if (activeLog)
    cli_text("Running cycle {cycle}/{args$nCycles}...")
  
      if (cycle == 1) {

        newParents <- selectNewParents(gen$F5, 5, "pheno")
        corMat[1,] = 0

        updatePheno(newParents,"NP")
        updateResults(2, newParents, "NP")

      } else {
        if (args$parentSelections == "F2"){
          newParents <- selectNewParents(gen$F2, 5, "ebv")
            if (grepl("rrblup",args$model)==TRUE){
            corMat[1,] = cor(bv(gen$F2), ebv(gen$F2)) #determine model performance
            }
  
            if (grepl("ann",args$model)==TRUE){
            corMat[1,] = cor(pheno(gen$F2), ebv(gen$F2)) 
            } #determine model performance
   
        }

        if (args$parentSelections == "F5"){
          newParents <- selectNewParents(gen$F5, 5, "ebv")
           if (grepl("rrblup",args$model)==TRUE){
           corMat[1,] = cor(bv(gen$F5), ebv(gen$F5)) #determine model performance
           }
  
           if (grepl("ann",args$model)==TRUE){
           corMat[1,] = cor(pheno(gen$F5), ebv(gen$F5)) 
           } #determine model performance
        }

        
        updatePheno()
        updateResults(2, newParents, "NP")
      }

   

  ## 200 random crosses of new parents
  
  gen$F1 = randCross(newParents, 200,nProgeny=3)

    updatePheno(gen$F1,"F1")  
    updateResults(3, gen$F1, "F1")
  
  
  ## self and bulk gen$F1 to form gen$F2 ##
  
  if (args$trainGen == "F2")
    trainModel()
  
  gen$F2 = self(gen$F1, nProgeny = 20) 
  gen$F2@ebv = getEBV(gen$F2)   ## set EBV 

  updateResults(4, gen$F2, "F2")
  updatePheno(gen$F2,"F2")


  
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
  updatePheno(gen$F3,"F3")
    

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
  updatePheno(gen$F4,"F4")
  


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
  updatePheno(gen$F5,"F5")
  


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
  updatePheno(gen$PYT,"PYT")
  

  
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
  updatePheno(gen$AYT,"AYT")
  
  ##set EBV using RRBLUP model##
  
  
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
  
  phenos = pheno(Variety)
  checkMat = as.data.frame(valuesMat)
  from = nrow(valuesMat) - sum(is.na(checkMat[,1])) +1
  to = from + nInd(Variety) -1
  phenos = pheno(Variety)
  gvs = gv(Variety)
  tbvs = bv(Variety)
  valuesMat[from:to,1] = rep("Variety", times=nInd(Variety))
  valuesMat[from:to,2] = phenos
  valuesMat[from:to,3] = gvs
  valuesMat[from:to,4] = tbvs
  valuesMat[from:to,5] = rep("NA", times=nInd(Variety))
    
  allelesMatVar <- getAllelesMat(Variety, "Variety")
  allelesMat <- rbind(allelesMat, allelesMatVar)
  
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
