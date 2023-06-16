## PEDIGREE BREEDING METHOD USING GEBVs TO SELECT
library(AlphaSimR)
library(rrBLUP)
source("FunctionsLibrary.R")

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
F1 = randCross(newParents, 200, nProgeny=3)

## self and bulk F1 to form F2 ##

F2 = self(F1, nProgeny = 30)
F2 = setPheno(F2)

## select top individuals from each family to form F2. Bulk and self to form F3
F3 = TopWithinFam(F2,10,100,"pheno")
F3 = setPheno(F3)

## select top individuals within F3 families to form F4 

F4 = TopWithinFam(F3,5,50,"pheno")
F4 = setPheno(F4)

## select top families from F4 to form F5 

F5 = TopFamily(F4,4,"pheno")
F5 = setPheno(F5)

## select top families from F5 for PYTs 

PYT = TopFamily(F5, 3,"pheno")
PYT = setPheno(PYT, reps=2)

gvMat[1,] <- mean(gv(PYT))
varMat[1,] <- varG(PYT)

## use PYTs as training data

TrainingGeno <- pullSegSiteGeno(PYT)
TrainingPheno <- pheno(PYT)

# source GS Prediction Model
source("rrblup_sc.R")

# calculate EBVs of PYTs
GetEBVrrblup(PYT) #get EBVs
PYT@ebv = EBV #set EBVs
corMat[1,] = cor(bv(PYT), ebv(PYT)) #determine model performance

# NEW CYCLE
for (cycle in 1:nCycles){
    ## select new parents from previous cycle PYTs
    newParents = selectInd(PYT, 10, use="ebv", top=TRUE)

    varMat[2,] = varG(newParents) #collect variance
    gvMat[2,] <- mean(gv(newParents)) #collect genetic values 
    allelesMatNP <- getAllelesMat(newParents, "NP") #collect genotypes

    ## 200 random crosses of new parents

    F1 = randCross(newParents, 200)
                                
    varMat[3,] = varG(F1)
    gvMat[3,] <- mean(gv(F1))
    allelesMatF1 <- getAllelesMat(F1, "F1")

    ## self and bulk F1 to form F2 ##

    F2 = self(F1, nProgeny = 30) 
                                
    varMat[4,] = varG(F2)
    gvMat[4,] <- mean(gv(F2))
    allelesMatF2 <- getAllelesMat(F2, "F2")

    ## set EBV using RRBLUP model

    GetEBVrrblup(F2)
    F2@ebv = EBV
    corMat[2,] = as.numeric(cor(bv(F2), ebv(F2)))

    ## select top individuals from F2 bulk to form F3 

    F3 = TopWithinFam(F2, 10, 100, "ebv")
    F3 = setPheno(F3)
                                
    varMat[5,] = varG(F3)
    gvMat[5,] <- mean(gv(F3))
    allelesMatF3 <- getAllelesMat(F3, "F3")

    ## set EBV using BLUP model

    GetEBVrrblup(F3)
    F3@ebv = EBV
    corMat[3,] = cor(bv(F3),ebv(F3))

    ## select top within familiy from F3 to form F4 
    F4 = TopWithinFam(F3, 5, 50, "ebv")
    F4 = setPheno(F4)
                                
    varMat[6,] = varG(F4)
    gvMat[6,] <- mean(gv(F4))                            
    allelesMatF4 <- getAllelesMat(F4, "F4")

    ##set EBV using BLUP model##
    GetEBVrrblup(F4)
    F4@ebv = EBV
    corMat[4,] = cor(bv(F4),ebv(F4))

    ## select top families from F4 to form F5 ##

    F5 = TopFamily(F4,4,"ebv")
    F5 = setPheno(F5)
                                
    varMat[7,]= varG(F5)
    gvMat[7,] <- mean(gv(F5))
    allelesMatF5 <- getAllelesMat(F5, "F5")

    #use F5 to retrain the model

    source("rrblup_sc_retrain.R")

    ##set EBV using RRBLUP model##
    GetEBVrrblup(F5)
    F5@ebv = EBV
    corMat[5,] = cor(bv(F5),ebv(F5))

    ## select top F5 families for preliminary yield trial ##
    PYT = TopFamily(F5,3,"ebv")
    PYT = setPheno(PYT, reps=2)
                                                        
    varMat[8,] = varG(PYT)
    gvMat[8,] <- mean(gv(PYT))
    allelesMatPYT <- getAllelesMat(PYT, "PYT")

    ##set EBV using RRBLUP model##
    GetEBVrrblup(PYT)
    PYT@ebv = EBV
    corMat[6,] = cor(bv(PYT),ebv(PYT))

    ## select top families from PYT for AYT ##

    AYT = TopFamily(PYT, 1, "ebv")
    AYT = setPheno(AYT, reps=5)
                                
    varMat[9,] = varG(AYT)
    gvMat[9,] <- mean(gv(AYT))
    allelesMatAYT <- getAllelesMat(AYT, "AYT")

    ##set EBV using RRBLUP model##
    GetEBVrrblup(AYT)
    AYT@ebv = EBV
    corMat[7,] = cor(bv(AYT),ebv(AYT))

    ## select top plants to form variety ##
    VarietySel = selectInd(AYT, 1, use="ebv")
    Variety = self(VarietySel)
    gvMat[10,] <- mean(gv(Variety))

    allelesMatVar <- getAllelesMat(Variety, "Variety")

    allelesMat <- rbind(allelesMatNP, allelesMatF1, allelesMatF2, allelesMatF3, allelesMatF4, allelesMatF5, allelesMatPYT, allelesMatAYT, allelesMatVar)


    ###collect bvs and ebvs###

    bvebv0 <- getBvEbv(newParents, "NP")
    bvebv1 <- getBvEbv(F2, "F2")
    bvebv2 <- getBvEbv(F3, "F3")
    bvebv3 <- getBvEbv(F4, "F4")
    bvebv4 <- getBvEbv(F5, "F5")
    bvebv5 <- getBvEbv(PYT, "PYT")
    bvebv6 <- getBvEbv(AYT, "AYT")

    bv_ebv_df <- as.data.frame(rbind(bvebv0,bvebv1,bvebv2,bvebv3,bvebv4,bvebv5,bvebv6))

    ###Select parents for next cycle

    newParents <- selectNewParents(F2,5,"ebv")

    geneticvalues[[cycle]][,rep] <- gvMat
    correlations[[cycle]][,rep] <- corMat
    variances[[cycle]][,rep] <- varMat
    alleles[[cycle]][[rep]] <- allelesMat
    bv_ebv[[cycle]][[rep]] <- bv_ebv_df
}