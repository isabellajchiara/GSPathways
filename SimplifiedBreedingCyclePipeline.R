## PEDIGREE BREEDING METHOD USING GEBVs TO SELECT
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
  bv_ebv = list(),
  pheno = list()
)

# create results matrices

gvMat <- matrix(nrow=nGen, ncol=1)
corMat <- matrix(nrow=nModels, ncol=1)
varMat <- matrix(nrow=nVar, ncol=1)
allelesMat <- NULL

if (args$parentSelections == "F2"){
nSim = 13366
}

if (args$parentSelections == "F5"){
nSim = 1486
}

phenoMat <- matrix(nrow=(60 + nSim*args$nCycles),ncol=3)

allTrainingDataGeno = list()
allTrainingDataPheno = list()

gen <- list()

if (activeLog)
  cli_text("Created outputs")

# establish simulation parameters
set.seed(1206)
defineTraitAEG(40,3.6,0.25) # nQtl per chr, mean,heritability

## FIRST CYCLE TO BUILD INITIAL TRAINING POP

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

# select top families from gen$F5 for PYTs 

gen$PYT1 = TopFamily(gen$F5, 2,"pheno") #select top 3 F5 families
gen$PYT1 = setPheno(gen$PYT1, reps=2)

if (activeLog)
  cli_text("Generated PYTs...")


# train model for initial parent selections 



if (args$parentSelections == "F2"){

gvMat[1,] <- mean(gv(gen$F2))
varMat[1,] <- varG(gen$F2)
phenoData = pheno(gen$F2)
tvs = bv(gen$F2)
phenoMat[1:nInd(gen$F2),1] = phenoData
phenoMat[1:nInd(gen$F2),2] = tbvs
phenoMat[1:nInd(gen$F2),3] = rep("ParentPool", times=nInd(gen$F2))
updateResults(1, F2, "ParentPool")

}

if (args$parentSelections == "F5"){

gvMat[1,] <- mean(gv(gen$F5))
varMat[1,] <- varG(gen$F5)
phenoData = pheno(gen$F5)
tvs = bv(gen$52)
phenoMat[1:nInd(gen$F5),1] = phenoData
phenoMat[1:nInd(gen$F5),2] = tbvs
phenoMat[1:nInd(gen$F5),2] = rep("ParentPool", times=nInd(gen$F5))
updateResults(1, F5, "ParentPool")

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

newParents <- selectNewParents(gen[[args$parentSelections]], 5, "pheno")

  
  updateResults(2, newParents, "NP")
  phenoData = pheno(newParents)
  checkMat = as.data.frame(phenoMat)
  from = nrow(phenoMat) - sum(is.na(checkMat[,1])) +1
  to = from + nInd(newParents) -1
  phenoMat[from:to,1] = phenoData
  phenoMat[from:to,2] = rep("NP", times=nInd(newParents))

if (activeLog)
  cli_text("Initial cycle complete")

# INITAL TRAINING POP IS BUILT, START NEW CYCLE. WE WILL CALL THIS CYCLE 1 
for (cycle in 1:args$nCycles){
  if (activeLog)
    cli_text("Running cycle {cycle}/{args$nCycles}...")
  
  
  ## 200 random crosses of new parents
  
  gen$F1 = randCross(newParents, 200,nProgeny=3)
  updateResults(3, gen$F1, "F1")
  updatePheno(gen$F1,"F1")
  
  
  ## self and bulk gen$F1 to form gen$F2 ##
  
  if (args$trainGen == "F2")
    trainModel()
  
  gen$F2 = self(gen$F1, nProgeny = 20) 
  updateResults(4, gen$F2, "F2")
  updatePheno(gen$F2,"F2")
  
  
  ## set EBV using RRBLUP model
  
  gen$F2@ebv = getEBV(gen$F2)
  
  if (args$model =="rrblup_random"){
    corMat[2,] = cor(bv(gen$PYT), ebv(gen$PYT)) #determine model performance
  }
  
  if (args$model == "ann_random"){
    corMat[2,] = cor(pheno(gen$PYT), ebv(gen$PYT)) 
  }#determine model performance
  
  ## select top individuals from gen$F2 bulk to form gen$F3 
  if (args$trainGen == "F3")
    trainModel()
  
  gen$F3 = TopWithinFam(gen$F2, 5,200 , "ebv")
  gen$F3 = setPheno(gen$F3)
  updateResults(5, gen$F3, "F3")
  updatePheno(gen$F3,"F3")
  
  
  
  ## set EBV using BLUP model
  
  gen$F3@ebv = getEBV(gen$F3)
  if (args$model =="rrblup_random"){
    corMat[3,] = cor(bv(gen$PYT), ebv(gen$PYT)) #determine model performance
  }
  
  if (args$model == "ann_random"){
    corMat[3,] = cor(pheno(gen$PYT), ebv(gen$PYT))
  } #determine model performance
  
  ## select top within familiy from gen$F3 to form gen$F4 
  
  if (args$trainGen == "F4")
    trainModel("F4", gen$F4)
  
  gen$F4 = TopWithinFam(gen$F3, 5, 30, "ebv")
  gen$F4 = setPheno(gen$F4)
  updateResults(6, gen$F4, "F4")
  updatePheno(gen$F4,"F4")
  
  
  ##set EBV using BLUP model##
  gen$F4@ebv = getEBV(gen$F4)
  if (args$model =="rrblup_random"){
    corMat[4,] = cor(bv(gen$PYT), ebv(gen$PYT)) #determine model performance
  }
  
  if (args$model == "ann_random"){
    corMat[4,] = cor(pheno(gen$PYT), ebv(gen$PYT)) #determine model performance
  }
  ## select top families from gen$F4 to form gen$F5 ##
  if (args$trainGen == "F5")
    trainModel()
  
  gen$F5 = TopFamily(gen$F4,4,"ebv")
  gen$F5 = setPheno(gen$F5)
  updateResults(7, gen$F5, "F5")
  updatePheno(gen$F5,"F5")
  
  ##set EBV using RRBLUP model##
  EBV <- getEBV(gen$F5)
  gen$F5@ebv = EBV
  if (args$model =="rrblup_random"){
    corMat[5,] = cor(bv(gen$PYT), ebv(gen$PYT)) #determine model performance
  }
  
  if (args$model == "ann_random"){
    corMat[5,] = cor(pheno(gen$PYT), ebv(gen$PYT)) #determine model performance
  }
  ## select top gen$F5 families for preliminary yield trial ##
  gen$PYT = TopFamily(gen$F5,2,"ebv")
  gen$PYT = setPheno(gen$PYT, reps=2)
  updateResults(8, gen$PYT, "PYT")
  updatePheno(gen$PYT,"PYT")
  
  ##set EBV using RRBLUP model##
  EBV <- getEBV(gen$PYT)
  gen$PYT@ebv = EBV
  
  if (args$model =="rrblup_random"){
    corMat[6,] = cor(bv(gen$PYT), ebv(gen$PYT)) #determine model performance
  }
  
  if (args$model == "ann_random"){
    corMat[6,] = cor(pheno(gen$PYT), ebv(gen$PYT)) 
  }#determine model performance
  ## select top families from gen$PYT for gen$AYT ##
  
  gen$AYT = TopFamily(gen$PYT, 1, "ebv")
  gen$AYT = setPheno(gen$AYT, reps=5)
  updateResults(9, gen$AYT, "AYT")
  updatePheno(gen$AYT,"AYT")
  
  ##set EBV using RRBLUP model##
  EBV <- getEBV(gen$AYT)
  gen$AYT@ebv = EBV
  
  if (args$model =="rrblup_random"){
    corMat[7,] = cor(bv(gen$PYT), ebv(gen$PYT)) #determine model performance
  }
  
  if (args$model == "ann_random"){
    corMat[7,] = cor(pheno(gen$PYT), ebv(gen$PYT)) 
  }#determine model performance
  
  ## select top plants to form variety ##
  VarietySel = selectInd(gen$AYT, 1, use="ebv")
  Variety = self(VarietySel)
  gvMat[10,] <- mean(gv(Variety))
  
  phenoData = pheno(Variety)
  checkMat = as.data.frame(phenoMat)
  from = nrow(phenoMat) - sum(is.na(checkMat[,1])) +1
  to = from + nInd(Variety) -1
  phenoMat[from:to,1] <- phenoData
  phenoMat[from:to,2] = rep("Variety", times=nInd(Variety))
  
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
  
  # collect bvs and ebvs
  
  bvebv0 <- getBvEbv(gen$PYT1, "NP")
  bvebv1 <- getBvEbv(gen$F2, "F2")
  bvebv2 <- getBvEbv(gen$F3, "F3")
  bvebv3 <- getBvEbv(gen$F4, "F4")
  bvebv4 <- getBvEbv(gen$F5, "F5")
  bvebv5 <- getBvEbv(gen$PYT, "PYT")
  bvebv6 <- getBvEbv(gen$AYT, "AYT")
  
  bv_ebv_df <- as.data.frame(rbind(bvebv0,bvebv1,bvebv2,bvebv3,bvebv4,bvebv5,bvebv6))
  
  ret$geneticvalues[[cycle]] <- gvMat
  ret$correlations[[cycle]] <- corMat
  ret$variances[[cycle]] <- varMat
  ret$alleles[[cycle]] <- allelesMat
  ret$pheno[[cycle]] <- phenoMat
  
  
  ret$bv_ebv[[cycle]] <- bv_ebv_df
}