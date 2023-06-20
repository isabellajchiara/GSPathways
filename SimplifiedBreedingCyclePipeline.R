## PEDIGREE BREEDING METHOD USING GEBVs TO SELECT
library(AlphaSimR)
library(rrBLUP)
library(caret)
library(ranger)
library(tidyverse)
library(e1071)
library(randomForest)
source("FunctionsLibrary.R")

stages <- list()

#Create Results Matrices

gvMat <- matrix(nrow=10, ncol=1)
corMat <- matrix(nrow=7, ncol=1)
varMat <- matrix(nrow=9, ncol=1)

# establish simulation parameters

genMap <- readRDS("genMapSNPs.RData") # can load other genMaps 
haplotypes <- readRDS("haplotypesSNPs.RData") # can load other genotype data, must match genMap

founderPop = newMapPop(genMap, 
                       haplotypes, 
                       inbred = FALSE, 
                       ploidy = 2L)

defineTraitAEG(10,8.8,0.25) # nQtl per chr, mean,heritability
#yield = (10,8.8,0.25)
#flowering = (5,35,0.8) and use defineTraitA

### FIRST CYCLE TO BUILD INITIAL TRAINING POP

## create base pop randomly cross 200 parents 

Base = newPop(founderPop)
Base = setPheno(Base)

newParents <- selectNewParents(Base,10,"pheno")
stages$F1 = randCross(newParents, 200, nProgeny=3)

## self and bulk stages$F1 to form stages$F2 ##

stages$F2 = self(stages$F1, nProgeny = 30)
stages$F2 = setPheno(stages$F2)

## select top individuals from each family to form stages$F2. Bulk and self to form stages$F3
stages$F3 = TopWithinFam(stages$F2,10,100,"pheno")
stages$F3 = setPheno(stages$F3)

## select top individuals within stages$F3 families to form stages$F4 

stages$F4 = TopWithinFam(stages$F3,5,50,"pheno")
stages$F4 = setPheno(stages$F4)

## select top families from stages$F4 to form stages$F5 

stages$F5 = TopFamily(stages$F4,4,"pheno")
stages$F5 = setPheno(stages$F5)

## select top families from stages$F5 for PYTs 

stages$PYT = TopFamily(stages$F5, 3,"pheno")
stages$PYT = setPheno(stages$PYT, reps=2)

gvMat[1,] <- mean(gv(stages$PYT))
varMat[1,] <- varG(stages$PYT)

## use PYTs as training data

TrainingGeno <- pullSegSiteGeno(stages$PYT)
TrainingPheno <- pheno(stages$PYT)

# source GS Prediction Model
source(fileTrain)

# calculate EBVs of PYTs
EBV <- getEBV(stages$PYT) #get EBVs
stages$PYT@ebv = EBV #set EBVs
corMat[1,] = cor(bv(stages$PYT), ebv(stages$PYT)) #determine model performance

# NEW CYCLE
for (cycle in 1:nCycles){
    ## select new parents from previous cycle PYTs
    
    if (cycle == 1) {
      newParents <- selectNewParents(stages$PYT, 5, "ebv")
    } else {
      newParents <- selectNewParents(stages$F2, 5, "ebv")
    }
  
    varMat[2,] = varG(newParents) #collect variance
    gvMat[2,] <- mean(gv(newParents)) #collect genetic values 
    allelesMatNP <- getAllelesMat(newParents, "NP") #collect genotypes

    ## 200 random crosses of new parents

    stages$F1 = randCross(newParents, 200)
                                
    varMat[3,] = varG(stages$F1)
    gvMat[3,] <- mean(gv(stages$F1))
    allelesMatF1 <- getAllelesMat(stages$F1, "F1")

    ## self and bulk stages$F1 to form stages$F2 ##

    stages$F2 = self(stages$F1, nProgeny = 30) 
    
    varMat[4,] = varG(stages$F2)
    gvMat[4,] <- mean(gv(stages$F2))
    allelesMatF2 <- getAllelesMat(stages$F2, "F2")

    if (trainStage == "F2")
      retrain()
      
    ## set EBV using RRBLUP model

    EBV <- getEBV(stages$F2)
    stages$F2@ebv = EBV
    corMat[2,] = as.numeric(cor(bv(stages$F2), ebv(stages$F2)))

    ## select top individuals from stages$F2 bulk to form stages$F3 

    stages$F3 = TopWithinFam(stages$F2, 10, 100, "ebv")
    stages$F3 = setPheno(stages$F3)
                                
    varMat[5,] = varG(stages$F3)
    gvMat[5,] <- mean(gv(stages$F3))
    allelesMatF3 <- getAllelesMat(stages$F3, "F3")

    if (trainStage == "F3")
      retrain()

    ## set EBV using BLUP model

    EBV <- getEBV(stages$F3)
    stages$F3@ebv = EBV
    corMat[3,] = cor(bv(stages$F3),ebv(stages$F3))

    ## select top within familiy from stages$F3 to form stages$F4 
    stages$F4 = TopWithinFam(stages$F3, 5, 50, "ebv")
    stages$F4 = setPheno(stages$F4)
                                
    varMat[6,] = varG(stages$F4)
    gvMat[6,] <- mean(gv(stages$F4))                            
    allelesMatF4 <- getAllelesMat(stages$F4, "F4")

    if (trainStage == "F4")
        retrain()

    ##set EBV using BLUP model##
    EBV <- getEBV(stages$F4)
    stages$F4@ebv = EBV
    corMat[4,] = cor(bv(stages$F4),ebv(stages$F4))

    ## select top families from stages$F4 to form stages$F5 ##

    stages$F5 = TopFamily(stages$F4,4,"ebv")
    stages$F5 = setPheno(stages$F5)

    varMat[7,]= varG(stages$F5)
    gvMat[7,] <- mean(gv(stages$F5))
    allelesMatF5 <- getAllelesMat(stages$F5, "F5")

    if (trainStage == "F5")
      retrain()

    ##set EBV using RRBLUP model##
    EBV <- getEBV(stages$F5)
    stages$F5@ebv = EBV
    corMat[5,] = cor(bv(stages$F5),ebv(stages$F5))

    ## select top stages$F5 families for preliminary yield trial ##
    stages$PYT = TopFamily(stages$F5,3,"ebv")
    stages$PYT = setPheno(stages$PYT, reps=2)
                                                        
    varMat[8,] = varG(stages$PYT)
    gvMat[8,] <- mean(gv(stages$PYT))
    allelesMatPYT <- getAllelesMat(stages$PYT, "PYT")

    ##set EBV using RRBLUP model##
    EBV <- getEBV(stages$PYT)
    stages$PYT@ebv = EBV
    corMat[6,] = cor(bv(stages$PYT),ebv(stages$PYT))

    ## select top families from stages$PYT for stages$stages$AYT ##

    stages$AYT = TopFamily(stages$PYT, 1, "ebv")
    stages$AYT = setPheno(stages$AYT, reps=5)
                                
    varMat[9,] = varG(stages$AYT)
    gvMat[9,] <- mean(gv(stages$AYT))
    allelesMatAYT <- getAllelesMat(stages$AYT, "AYT")

    ##set EBV using RRBLUP model##
    EBV <- getEBV(stages$AYT)
    stages$AYT@ebv = EBV
    corMat[7,] = cor(bv(stages$AYT),ebv(stages$AYT))

    ## select top plants to form variety ##
    VarietySel = selectInd(stages$AYT, 1, use="ebv")
    Variety = self(VarietySel)
    gvMat[10,] <- mean(gv(Variety))

    allelesMatVar <- getAllelesMat(Variety, "Variety")

    allelesMat <- rbind(allelesMatNP, allelesMatF1, allelesMatF2, allelesMatF3, allelesMatF4, allelesMatF5, allelesMatPYT, allelesMatAYT, allelesMatVar)


    ###collect bvs and ebvs###

    bvebv0 <- getBvEbv(newParents, "NP")
    bvebv1 <- getBvEbv(stages$F2, "F2")
    bvebv2 <- getBvEbv(stages$F3, "F3")
    bvebv3 <- getBvEbv(stages$F4, "F4")
    bvebv4 <- getBvEbv(stages$F5, "F5")
    bvebv5 <- getBvEbv(stages$PYT, "PYT")
    bvebv6 <- getBvEbv(stages$AYT, "AYT")

    bv_ebv_df <- as.data.frame(rbind(bvebv0,bvebv1,bvebv2,bvebv3,bvebv4,bvebv5,bvebv6))

    geneticvalues[[cycle]][,rep] <- gvMat
    correlations[[cycle]][,rep] <- corMat
    variances[[cycle]][,rep] <- varMat
    alleles[[cycle]][[rep]] <- allelesMat
    bv_ebv[[cycle]][[rep]] <- bv_ebv_df
}
