## PEDIGREE BREEDING METHOD USING GEBVs TO SELECT
suppressMessages(library(AlphaSimR))
source("ParameterSettings.R")
source("InterfaceLibrary.R")
source("FunctionsLibrary.R")
source("ModelVariables.R")
loadModelLibs()

# Log only if reps are being run serially
activeLog <- args$nCores == 1 || modelParallelism

if (activeLog)
  cli_text("Generating parent population...")

# Data to be returned
ret <- list( 
  geneticvalues = list(),
  correlations = list(),
  variances = list(),
  alleles = list(),
  bv_ebv = list()
)

#Create Results Matrices

gvMat <- matrix(nrow=10, ncol=1)
corMat <- matrix(nrow=7, ncol=1)
varMat <- matrix(nrow=9, ncol=1)
allelesMat <- NULL

gen <- list()

# establish simulation parameters


defineTraitAEG(10,8.8,0.25) # nQtl per chr, mean,heritability
#yield = (10,8.8,0.25)
#flowering = (5,35,0.8) and use defineTraitA

### FIRST CYCLE TO BUILD INITIAL TRAINING POP

## create base pop randomly cross 200 parents 

Base = newPop(founderPop)
Base = setPheno(Base)

newParents <- selectNewParents(Base,10,"pheno")
gen$F1 = randCross(newParents, 200, nProgeny=3)

## self and bulk gen$F1 to form gen$F2 ##

gen$F2 = self(gen$F1, nProgeny = 30)
gen$F2 = setPheno(gen$F2)

## select top individuals from each family to form gen$F2. Bulk and self to form gen$F3
gen$F3 = TopWithinFam(gen$F2,10,100,"pheno") # top 10 F2 fam, 100 ind per fam using pheno
gen$F3 = setPheno(gen$F3)

## select top individuals within gen$F3 families to form gen$F4 

gen$F4 = TopWithinFam(gen$F3,5,50,"pheno") # top 5 F3 fam, 50 ind per fam using pheno
gen$F4 = setPheno(gen$F4)

## select top families from gen$F4 to form gen$F5 

gen$F5 = TopFamily(gen$F4,4,"pheno") #select top 4 F4 families 
gen$F5 = setPheno(gen$F5)

## select top families from gen$F5 for PYTs 

gen$PYT = TopFamily(gen$F5, 3,"pheno") #select top 3 F5 families
gen$PYT = setPheno(gen$PYT, reps=2)

gvMat[1,] <- mean(gv(gen$PYT))
varMat[1,] <- varG(gen$PYT)

## use PYTs as training data and GS Prediction Model
trainModel()

# calculate EBVs of PYTs
EBV <- getEBV(gen$PYT) #get EBVs
gen$PYT@ebv = EBV #set EBVs
corMat[1,] = cor(bv(gen$PYT), ebv(gen$PYT)) #determine model performance

# INITAL TRAINING POP IS BUILT, START NEW CYCLE. WE WILL CALL THIS CYCLE 1 
for (cycle in 1:args$nCycles){
  if (activeLog)
    cli_text("Running cycle {cycle}/{args$nCycles}...")

  ## select new parents from previous cycle PYTs
  
  if (cycle == 1) {
    newParents <- selectNewParents(gen$PYT, 5, "ebv")
  } else {
    
    newParents <- selectNewParents(gen[[args$parentSelections]], 5, "ebv")
  }

  updateResults(2, newParents, "NP")
  ## 200 random crosses of new parents

  gen$F1 = randCross(newParents, 200)
  updateResults(3, gen$F1, "F1")
                              
  ## self and bulk gen$F1 to form gen$F2 ##
  
  if (args$trainGen == "F2")
    trainModel()
  
  gen$F2 = self(gen$F1, nProgeny = 20) 
  updateResults(4, gen$F2, "F2")
  
    
  ## set EBV using RRBLUP model

  gen$F2@ebv = getEBV(gen$F2)
  corMat[2,] = as.numeric(cor(bv(gen$F2), ebv(gen$F2)))

  ## select top individuals from gen$F2 bulk to form gen$F3 
  if (args$trainGen == "F3")
    trainModel()
  
  gen$F3 = TopWithinFam(gen$F2, 5,200 , "ebv")
  gen$F3 = setPheno(gen$F3)
  updateResults(5, gen$F3, "F3")

  

  ## set EBV using BLUP model

  gen$F3@ebv = getEBV(gen$F3)
  corMat[3,] = cor(bv(gen$F3),ebv(gen$F3))

  ## select top within familiy from gen$F3 to form gen$F4 
  
  if (args$trainGen == "F4")
    trainModel("F4", gen$F4)
  
  gen$F4 = TopWithinFam(gen$F3, 5, 30, "ebv")
  gen$F4 = setPheno(gen$F4)
  updateResults(6, gen$F4, "F4")


  ##set EBV using BLUP model##
  gen$F4@ebv = getEBV(gen$F4)
  corMat[4,] = cor(bv(gen$F4),ebv(gen$F4))

  ## select top families from gen$F4 to form gen$F5 ##
  if (args$trainGen == "F5")
    trainModel()
  
  gen$F5 = TopFamily(gen$F4,3,"ebv")
  gen$F5 = setPheno(gen$F5)
  updateResults(7, gen$F5, "F5")


  ##set EBV using RRBLUP model##
  EBV <- getEBV(gen$F5)
  gen$F5@ebv = EBV
  corMat[5,] = cor(bv(gen$F5),ebv(gen$F5))

  ## select top gen$F5 families for preliminary yield trial ##
  gen$PYT = TopFamily(gen$F4,2,"ebv")
  gen$PYT = setPheno(gen$PYT, reps=2)
  updateResults(8, gen$PYT, "PYT")

  ##set EBV using RRBLUP model##
  EBV <- getEBV(gen$PYT)
  gen$PYT@ebv = EBV
  corMat[6,] = cor(bv(gen$PYT),ebv(gen$PYT))

  ## select top families from gen$PYT for gen$AYT ##

  gen$AYT = TopFamily(gen$PYT, 1, "ebv")
  gen$AYT = setPheno(gen$AYT, reps=5)
  updateResults(9, gen$AYT, "AYT")

  ##set EBV using RRBLUP model##
  EBV <- getEBV(gen$AYT)
  gen$AYT@ebv = EBV
  corMat[7,] = cor(bv(gen$AYT),ebv(gen$AYT))

  ## select top plants to form variety ##
  VarietySel = selectInd(gen$AYT, 1, use="ebv")
  Variety = self(VarietySel)
  gvMat[10,] <- mean(gv(Variety))

  allelesMatVar <- getAllelesMat(Variety, "Variety")
  allelesMat <- rbind(allelesMat, allelesMatVar)

  ###collect bvs and ebvs###

  bvebv0 <- getBvEbv(newParents, "NP")
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
  ret$bv_ebv[[cycle]] <- bv_ebv_df
}

