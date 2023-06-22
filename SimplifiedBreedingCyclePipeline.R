## PEDIGREE BREEDING METHOD USING GEBVs TO SELECT
library(AlphaSimR)
library(rrBLUP)
library(caret)
library(ranger)
library(tidyverse)
library(e1071)
library(randomForest)
library(foreach)
library(import)
library(doParallel)

gens <- list()

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
gens$F1 = randCross(newParents, 200, nProgeny=3)

## self and bulk gens$F1 to form gens$F2 ##

gens$F2 = self(gens$F1, nProgeny = 30)
gens$F2 = setPheno(gens$F2)

## select top individuals from each family to form gens$F2. Bulk and self to form gens$F3
gens$F3 = TopWithinFam(gens$F2,10,100,"pheno")
gens$F3 = setPheno(gens$F3)

## select top individuals within gens$F3 families to form gens$F4 

gens$F4 = TopWithinFam(gens$F3,5,50,"pheno")
gens$F4 = setPheno(gens$F4)

## select top families from gens$F4 to form gens$F5 

gens$F5 = TopFamily(gens$F4,4,"pheno")
gens$F5 = setPheno(gens$F5)

## select top families from gens$F5 for PYTs 

gens$PYT = TopFamily(gens$F5, 3,"pheno")
gens$PYT = setPheno(gens$PYT, reps=2)

gvMat[1,] <- mean(gv(gens$PYT))
varMat[1,] <- varG(gens$PYT)

## use PYTs as training data and GS Prediction Model
trainModel("PYT")

# calculate EBVs of PYTs
EBV <- getEBV(gens$PYT) #get EBVs
gens$PYT@ebv = EBV #set EBVs
corMat[1,] = cor(bv(gens$PYT), ebv(gens$PYT)) #determine model performance

# NEW CYCLE
for (cycle in 1:nCycles){
    ## select new parents from previous cycle PYTs
    
    if (cycle == 1) {
      newParents <- selectNewParents(gens$PYT, 5, "ebv")
    } else {
      newParents <- selectNewParents(gens$F2, 5, "ebv")
    }
  
    varMat[2,] = varG(newParents) #collect variance
    gvMat[2,] <- mean(gv(newParents)) #collect genetic values 
    allelesMatNP <- getAllelesMat(newParents, "NP") #collect genotypes

    ## 200 random crosses of new parents

    gens$F1 = randCross(newParents, 200)
                                
    varMat[3,] = varG(gens$F1)
    gvMat[3,] <- mean(gv(gens$F1))
    allelesMatF1 <- getAllelesMat(gens$F1, "F1")

    ## self and bulk gens$F1 to form gens$F2 ##

    gens$F2 = self(gens$F1, nProgeny = 30) 
    
    varMat[4,] = varG(gens$F2)
    gvMat[4,] <- mean(gv(gens$F2))
    allelesMatF2 <- getAllelesMat(gens$F2, "F2")

    if (trainGen == "F2")
      trainModel(trainGen)
      
    ## set EBV using RRBLUP model

    EBV <- getEBV(gens$F2)
    gens$F2@ebv = EBV
    corMat[2,] = as.numeric(cor(bv(gens$F2), ebv(gens$F2)))

    ## select top individuals from gens$F2 bulk to form gens$F3 

    gens$F3 = TopWithinFam(gens$F2, 10, 100, "ebv")
    gens$F3 = setPheno(gens$F3)
                                
    varMat[5,] = varG(gens$F3)
    gvMat[5,] <- mean(gv(gens$F3))
    allelesMatF3 <- getAllelesMat(gens$F3, "F3")

    if (trainGen == "F3")
      trainModel(trainGen)

    ## set EBV using BLUP model

    EBV <- getEBV(gens$F3)
    gens$F3@ebv = EBV
    corMat[3,] = cor(bv(gens$F3),ebv(gens$F3))

    ## select top within familiy from gens$F3 to form gens$F4 
    gens$F4 = TopWithinFam(gens$F3, 5, 50, "ebv")
    gens$F4 = setPheno(gens$F4)
                                
    varMat[6,] = varG(gens$F4)
    gvMat[6,] <- mean(gv(gens$F4))                            
    allelesMatF4 <- getAllelesMat(gens$F4, "F4")

    if (trainGen == "F4")
      trainModel(trainGen)

    ##set EBV using BLUP model##
    EBV <- getEBV(gens$F4)
    gens$F4@ebv = EBV
    corMat[4,] = cor(bv(gens$F4),ebv(gens$F4))

    ## select top families from gens$F4 to form gens$F5 ##

    gens$F5 = TopFamily(gens$F4,4,"ebv")
    gens$F5 = setPheno(gens$F5)

    varMat[7,]= varG(gens$F5)
    gvMat[7,] <- mean(gv(gens$F5))
    allelesMatF5 <- getAllelesMat(gens$F5, "F5")

    if (trainGen == "F5")
      trainModel(trainGen)

    ##set EBV using RRBLUP model##
    EBV <- getEBV(gens$F5)
    gens$F5@ebv = EBV
    corMat[5,] = cor(bv(gens$F5),ebv(gens$F5))

    ## select top gens$F5 families for preliminary yield trial ##
    gens$PYT = TopFamily(gens$F5,3,"ebv")
    gens$PYT = setPheno(gens$PYT, reps=2)
                                                        
    varMat[8,] = varG(gens$PYT)
    gvMat[8,] <- mean(gv(gens$PYT))
    allelesMatPYT <- getAllelesMat(gens$PYT, "PYT")

    ##set EBV using RRBLUP model##
    EBV <- getEBV(gens$PYT)
    gens$PYT@ebv = EBV
    corMat[6,] = cor(bv(gens$PYT),ebv(gens$PYT))

    ## select top families from gens$PYT for gens$gens$AYT ##

    gens$AYT = TopFamily(gens$PYT, 1, "ebv")
    gens$AYT = setPheno(gens$AYT, reps=5)
                                
    varMat[9,] = varG(gens$AYT)
    gvMat[9,] <- mean(gv(gens$AYT))
    allelesMatAYT <- getAllelesMat(gens$AYT, "AYT")

    ##set EBV using RRBLUP model##
    EBV <- getEBV(gens$AYT)
    gens$AYT@ebv = EBV
    corMat[7,] = cor(bv(gens$AYT),ebv(gens$AYT))

    ## select top plants to form variety ##
    VarietySel = selectInd(gens$AYT, 1, use="ebv")
    Variety = self(VarietySel)
    gvMat[10,] <- mean(gv(Variety))

    allelesMatVar <- getAllelesMat(Variety, "Variety")

    allelesMat <- rbind(allelesMatNP, allelesMatF1, allelesMatF2, allelesMatF3, allelesMatF4, allelesMatF5, allelesMatPYT, allelesMatAYT, allelesMatVar)


    ###collect bvs and ebvs###

    bvebv0 <- getBvEbv(newParents, "NP")
    bvebv1 <- getBvEbv(gens$F2, "F2")
    bvebv2 <- getBvEbv(gens$F3, "F3")
    bvebv3 <- getBvEbv(gens$F4, "F4")
    bvebv4 <- getBvEbv(gens$F5, "F5")
    bvebv5 <- getBvEbv(gens$PYT, "PYT")
    bvebv6 <- getBvEbv(gens$AYT, "AYT")

    bv_ebv_df <- as.data.frame(rbind(bvebv0,bvebv1,bvebv2,bvebv3,bvebv4,bvebv5,bvebv6))

    geneticvalues[[cycle]][,rep] <- gvMat
    correlations[[cycle]][,rep] <- corMat
    variances[[cycle]][,rep] <- varMat
    alleles[[cycle]][[rep]] <- allelesMat
    bv_ebv[[cycle]][[rep]] <- bv_ebv_df
}
