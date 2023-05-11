#load packages
library(AlphaSimR)
library(rrBLUP)

#establish simulation parameters
genMap <- readRDS("genMapSNPs.RData")
haplotypes <- readRDS("haplotypesSNPs.RData")

founderPop = newMapPop(genMap, 
                       haplotypes, 
                       inbred = FALSE, 
                       ploidy = 2L)

SP <- SimParam$new(founderPop)
SP$addTraitAEG(10, mean=900)
SP$setVarE(h2=0.25)

#build training pop

## randomly cross 200 parents 
Parents = newPop(founderPop)
F1 = randCross(Parents, 200)

## self and bulk F1 to form F2 ##

F2 = self(F1, nProgeny = 10)
F2 = setPheno(F2)

## select top 100 individuals from F2 bulk and self to form F3
F3Sel = selectInd(F2, 100, use="pheno", top=TRUE) 
F3 = self(F3Sel, nProgeny=30)
F3 = setPheno(F3)

##select top individuals within F3 families to form F4 ##

F4Sel = selectWithinFam(F3, 10, use="pheno", top=TRUE) 
F4 = self(F4Sel)
F4 = setPheno(F4)

## select top families from F4 to form F5 ##

F5Sel = selectFam(F4, 30, use="pheno", top=TRUE)
F5 = self(F5Sel)
F5 = setPheno(F5)


## select top families from F5 for PYTs ##

PYTSel = selectFam(F5, 16, use="pheno", top=TRUE) 
PYT = self(PYTSel, nProgeny = 2)
PYT = setPheno(PYT, reps=2)

## self PYT for sufficiently large TP ##

TP = self(PYT, nProgeny = 4)

TrainingGeno <- pullSegSiteGeno(TP)
TrainingPheno <- pheno(TP)


#load in GS prediction model and use it to select parents for next cycle

#source GS Prediction Model
source("rrBLUP_RD.R")

M_PYT <-pullSegSiteGeno(PYT)
G_PYT = M_PYT-1
EBVPYT <- G_PYT %*% markerEffects

PYT@ebv <- as.matrix(EBVPYT)

newParents = selectInd(PYT, 10, use="ebv", top=TRUE)

F1 = randCross(newParents, 100) 
var1 <-varG(F1)

## self and bulk F1 to form F2 ##

F2 = self(F1, nProgeny = 5) ##nProgeny = number of progeny per cross## 
var2 <-varG(F2)

##set EBV using BLUP model##
M_F2 <-pullSegSiteGeno(F2)
G_F2 = M_F2-1
EBVF2 <- G_F2 %*% markerEffects

F2@ebv <- as.matrix(EBVF2)

cor2 = cor(gv(F2), ebv(F2))

source("SelectParentsF2.R")

## select top families to form F3 ##

F3Sel = selectWithinFam(F2, 6, use="ebv", top=TRUE) 
F3 = self(F3Sel, nProgeny = 6)

##set EBV using BLUP model##
M_F3 <-pullSegSiteGeno(F3)
G_F3 = M_F3-1
EBVF3 <- G_F3 %*% markerEffects

F3@ebv <- as.matrix(EBVF3)

cor3 = cor(gv(F3),ebv(F3))
var3 = varG(F3)


##select top families from F3 to form F4 ##

F4Sel = selectWithinFam(F3, 10, use="ebv", top=TRUE) 
F4 = self(F4Sel)

##set EBV using BLUP model##
M_F4 <-pullSegSiteGeno(F4)
G_F4 = M_F4-1
EBVF4 <- G_F4 %*% markerEffects

F4@ebv <- as.matrix(EBVF4)

cor4 = cor(gv(F4),ebv(F4))
var4 = varG(F4)

## select top families from F4 to form F5 ##

F5Sel = selectFam(F4, 10, use="ebv")
F5 = self(F5Sel, nProgeny=3)

source("RF_RD_RetrainF2F5.R")

#continue pipeline

##set EBV using BLUP model##
M = as.data.frame(pullSegSiteGeno(F5))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
F5ebv <- as.numeric(predict(rf_fit2, M))
F5@ebv <- as.matrix(F5ebv)

var5 = varG(F5)
cor5 = cor(gv(F5),ebv(F5))


## select top F5 families for preliminary yield trial ##

PYTSel = selectFam(F5, 4, use="ebv") 
PYT = self(PYTSel)
var6 = varG(PYT)

##set EBV using BLUP model##
M = as.data.frame(pullSegSiteGeno(PYT))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
PYTebv <- as.numeric(predict(rf_fit2, M))
PYT@ebv <- as.matrix(PYTebv)

cor6 = cor(gv(PYT),ebv(PYT))

## select top families from PYT for AYT ##

AYTSel = selectFam(PYT,  3, use="ebv", reps=5, top=TRUE) 
AYT = self(AYTSel)
var7 = varG(AYT)

##set EBV using BLUP model##
M = as.data.frame(pullSegSiteGeno(AYT))
colnames(M) <- paste("ID",2:(ncol(M)+1),sep="")
AYTebv <- as.numeric(predict(rf_fit2, M))
AYT@ebv <- as.matrix(AYTebv)

cor7 = cor(gv(AYT),ebv(AYT))

## select top plants to form variety ##
VarietySel = selectInd(AYT, 1, use="ebv")
Variety = self(VarietySel)
var8 = varG(Variety)

## pull genetic value for each generation and write results ##


F1gv1 <- as.vector(mean(gv(F1)))
F2gv1 <- as.vector(mean(gv(F2)))
F3gv1 <- as.vector(mean(gv(F3)))
F4gv1 <- as.vector(mean(gv(F4)))
F5gv1 <- as.vector(mean(gv(F5)))
PYTgv1 <- as.vector(mean(gv(PYT)))
AYTgv1 <- as.vector(mean(gv(AYT)))
Varietygv1 <- as.vector(mean(gv(Variety)))

cycle1gv <- matrix(nrow=8, ncol=1)
cycle1gv[1,] <- F1gv1
cycle1gv[2,] <- F2gv1
cycle1gv[3,] <- F3gv1
cycle1gv[4,] <- F4gv1
cycle1gv[5,] <- F5gv1
cycle1gv[6,] <- PYTgv1
cycle1gv[7,] <- AYTgv1
cycle1gv[8,] <- Varietygv1

cycle1gvs <- as.matrix(cycle1gv)




###list correlations to view model performacne ##
corMat <- matrix(nrow=6, ncol=1)
corMat[1,] <- cor2
corMat[2,] <- cor3
corMat[3,] <- cor4
corMat[4,] <- cor5
corMat[5,] <- cor6
corMat[6,] <- cor7

cycle1cors <- as.matrix(corMat)


varMat <- matrix(nrow=7, ncol=1)
varMat[1,] <- var1
varMat[2,] <- var2
varMat[3,] <- var3
varMat[4,] <- var4
varMat[5,] <- var5
varMat[6,] <- var6
varMat[7,] <- var7


cycle1vars <- as.matrix(varMat)


source("1CycleTwoF2.R")

