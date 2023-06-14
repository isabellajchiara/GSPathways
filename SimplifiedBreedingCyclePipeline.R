#Pedigree Breeding Method using GEBV Selections

#Create Results Matrices

gvMat <- matrix(nrow=10, ncol=1)
corMat <- matrix(nrow=7, ncol=1)
varMat <- matrix(nrow=9, ncol=1)

#establish simulation parameters
genMap <- readRDS("genMapSNPs.RData")
haplotypes <- readRDS("haplotypesSNPs.RData")

founderPop = newMapPop(genMap, 
                       haplotypes, 
                       inbred = FALSE, 
                       ploidy = 2L)

defineTraitAEG(10,8.8,0.25) 
#yield = (10,8.8,0.25)
#flowering = (5,35,0.8)

#build training pop

## randomly cross 200 parents 
Base = newPop(founderPop)
Base = setPheno(Base)

selectNewParentsPheno(Base,10)
F1 = randCross(NewParents, 200, nProgeny=3)

## self and bulk F1 to form F2 ##

F2 = self(F1, nProgeny = 30)
F2 = setPheno(F2)

## select top individuals from each family to form F2. Bulk and self to form F3
F3 = TopWithinFamPheno(F2,10,100)
F3 = setPheno(F3)

##select top individuals within F3 families to form F4 ##

F4 = TopWithinFamPheno(F3,5,50)
F4 = setPheno(F4)

## select top families from F4 to form F5 ##

F5 = TopFamilyPheno(F4,4)
F5 = setPheno(F5)

## select top families from F5 for PYTs ##

PYT = TopFamilyPheno(F5, 3) 
PYT = setPheno(PYT, reps=2)

gvMat[1,] <- mean(gv(PYT))
varMat[1,] <- varG(PYT)

## use PYTs as training data

TrainingGeno <- pullSegSiteGeno(PYT)
TrainingPheno <- pheno(PYT)

#source GS Prediction Model
source("rrblup_sc.R")


GetEBVrrblup(PYT)
PYT@ebv = EBV

corMat[1,] = cor(bv(PYT), ebv(PYT))

newParents = selectInd(PYT, 10, use="ebv", top=TRUE)
varMat[2,] = varG(newParents)
gvMat[2,] <- mean(gv(newParents))

allelesMatNP <- pullSegSiteHaplo(newParents)
Gen <- as.data.frame(rep("NP", times=nInd(newParents)))
colnames(Gen) <- "Gen"
allelesMatNP <- cbind(Gen, allelesMatNP)

#start new cycle

##start with 200 random crosses

F1 = randCross(newParents, 200) 
varMat[3,] = varG(F1)
gvMat[3,] <- mean(gv(F1))

allelesMatF1 <- pullSegSiteHaplo(F1)
Gen <- as.data.frame(rep("F1", times=nInd(F1)))
colnames(Gen) <- "Gen"
allelesMatF1 <- cbind(Gen, allelesMatF1)

## self and bulk F1 to form F2 ##

F2 = self(F1, nProgeny = 30) 
varMat[4,] = varG(F2)
gvMat[4,] <- mean(gv(F2))

allelesMatF2 <- pullSegSiteHaplo(F2)
Gen <- as.data.frame(rep("F2", times=nInd(F2)))
colnames(Gen) <- "Gen"
allelesMatF2 <- cbind(Gen, allelesMatF2)

##set EBV using BLUP model##
GetEBVrrblup(F2)
F2@ebv = EBV
corMat[2] = as.numeric(cor(bv(F2), ebv(F2)))

## select top individuals from F2 bulk  to form F3 ##

TopFamF2 = selectFam(F2, 10, use="pheno", top=TRUE) 
SelectionsF2 = selectWithinFam(TopFamF2, 100, use="ebv", top=TRUE)

F3 = self(SelectionsF2)
F3 = setPheno(F3)
varMat[5,] = varG(F3)
gvMat[5,] <- mean(gv(F3))

allelesMatF3 <- pullSegSiteHaplo(F3)
Gen <- as.data.frame(rep("F3", times=nInd(F3)))
colnames(Gen) <- "Gen"
allelesMatF3 <- cbind(Gen, allelesMatF3)



##set EBV using BLUP model##
GetEBVrrblup(F3)
F3@ebv = EBV
corMat[3,] = cor(bv(F3),ebv(F3))

##select top within familiy from F3 to form F4 ##

TopFamF3 = selectFam(F3,5,use="pheno", top=TRUE) 
SelectionsF3 = selectWithinFam(TopFamF3, 50, use="pheno", top=TRUE)

F4 = self(SelectionsF3)
F4 = setPheno(F4)
varMat[6,] = varG(F4)
gvMat[6,] <- mean(gv(F4))

allelesMatF4 <- pullSegSiteHaplo(F4)
Gen <- as.data.frame(rep("F4", times=nInd(F4)))
colnames(Gen) <- "Gen"
allelesMatF4 <- cbind(Gen, allelesMatF4)


##set EBV using BLUP model##
GetEBVrrblup(F4)
F4@ebv = EBV
corMat[4,] = cor(bv(F4),ebv(F4))

## select top families from F4 to form F5 ##

SelectionsF4 = selectFam(F4, 4, use="ebv")
F5 = self(SelectionsF4)
varMat[7,]= varG(F5)
gvMat[7,] <- mean(gv(F5))

allelesMatF5 <- pullSegSiteHaplo(F5)
Gen <- as.data.frame(rep("F5", times=nInd(F5)))
colnames(Gen) <- "Gen"
allelesMatF5 <- cbind(Gen, allelesMatF5)

#use F5 to retrain the model

source("rrblup_sc_retrain.R")

#continue pipeline

##set EBV using BLUP model##
GetEBVrrblup(F5)
F5@ebv = EBV
corMat[5,] = cor(bv(F5),ebv(F5))

## select top F5 families for preliminary yield trial ##

SelectionsF5 = selectFam(F5, 3, use="ebv") 
PYT = self(SelectionsF5)
varMat[8,] = varG(PYT)
gvMat[8,] <- mean(gv(PYT))

allelesMatPYT <- pullSegSiteHaplo(PYT)
Gen <- as.data.frame(rep("PYT", times=nInd(PYT)))
colnames(Gen) <- "Gen"
allelesMatPYT <- cbind(Gen, allelesMatPYT)

##set EBV using BLUP model##
GetEBVrrblup(PYT)
PYT@ebv = EBV
corMat[6,] = cor(bv(PYT),ebv(PYT))

## select top families from PYT for AYT ##

SelectionsPYT = selectFam(PYT,  1, use="ebv", reps=5, top=TRUE) 
AYT = self(SelectionsPYT)
varMat[9,] = varG(AYT)
gvMat[9,] <- mean(gv(AYT))

allelesMatAYT <- pullSegSiteHaplo(AYT)
Gen <- as.data.frame(rep("AYT", times=nInd(AYT)))
colnames(Gen) <- "Gen"
allelesMatAYT <- cbind(Gen, allelesMatAYT)

##set EBV using BLUP model##
GetEBVrrblup(AYT)
AYT@ebv = EBV
corMat[7,] = cor(bv(AYT),ebv(AYT))

## select top plants to form variety ##
VarietySel = selectInd(AYT, 1, use="ebv")
Variety = self(VarietySel)
gvMat[10,] <- mean(gv(Variety))

allelesMatVar <- pullSegSiteHaplo(Variety)
Gen <- as.data.frame(rep("Variety", times=nInd(Variety)))
colnames(Gen) <- "Gen"
allelesMatVar <- cbind(Gen, allelesMatVar)

allelesMat <- rbind(allelesMatNP, allelesMatF1, allelesMatF2, allelesMatF3, allelesMatF4, allelesMatF5, allelesMatPYT, allelesMatAYT, allelesMatVar)


###collect bvs and ebvs###

bvebv0 <- cbind(bv(newParents), ebv(newParents))
Gen <- as.data.frame(rep("NP", times=nInd(newParents)))
bvebv0 <- cbind(Gen, bvebv0)
colnames(bvebv0) <- c("Gen","bv","ebv")

bvebv1 <- cbind(bv(F2), ebv(F2))
Gen <- as.data.frame(rep("F2", times=nInd(F2)))
bvebv1 <- cbind(Gen, bvebv1)
colnames(bvebv1) <- c("Gen","bv","ebv")

bvebv2 <- cbind(bv(F3), ebv(F3))
Gen <- as.data.frame(rep("F3", times=nInd(F3)))
bvebv2 <- cbind(Gen, bvebv2)
colnames(bvebv2) <- c("Gen","bv","ebv")

bvebv3 <- cbind(bv(F4), ebv(F4))
Gen <- as.data.frame(rep("F4", times=nInd(F4)))
bvebv3 <- cbind(Gen, bvebv3)
colnames(bvebv3) <- c("Gen","bv","ebv")

bvebv4 <- cbind(bv(F5), ebv(F5))
Gen <- as.data.frame(rep("F5", times=nInd(F5)))
bvebv4 <- cbind(Gen, bvebv4)
colnames(bvebv4) <- c("Gen","bv","ebv")

bvebv5 <- cbind(bv(PYT), ebv(PYT))
Gen <- as.data.frame(rep("PYT", times=nInd(PYT)))
bvebv5 <- cbind(Gen, bvebv5)
colnames(bvebv5) <- c("Gen","bv","ebv")

bvebv6 <- cbind(bv(AYT), ebv(AYT))
Gen <- as.data.frame(rep("AYT", times=nInd(AYT)))
bvebv6 <- cbind(Gen, bvebv6)
colnames(bvebv6) <- c("Gen","bv","ebv")


bv_ebv <- as.data.frame(rbind(bvebv0,bvebv1,bvebv2,bvebv3,bvebv4,bvebv5,bvebv6))


###Select parents for next cycle

selectNewParentsEBV(F2,5)
